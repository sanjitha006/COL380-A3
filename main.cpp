/*
 * MPI Parallel Branch-and-Bound — Budgeted Maximum Clique Problem
 *
 * Architecture: DECENTRALIZED WORK STEALING
 * ─────────────
 * • No Master Rank: All ranks perform computation.
 * • Work Stealing: Idle ranks send TAG_STEAL_REQ to a randomly chosen peer.
 * • Termination Detection: Dijkstra's Dual-Pass Token Ring algorithm.
 * • Strict Sequential Bounding: Retains unmodified, strict bounding heuristics
 *   to comply with assignment constraints.
 *
 * Sequential Optimisations (data-structure level, no algorithmic change):
 * ─────────────────────────────────────────────────────────────────────────
 * [OPT-1] BITSET ADJACENCY
 *   adj[] is stored as packed uint64_t words (NWORDS words per vertex) instead
 *   of a uint8_t flat array.  One 64-byte cache line now covers 512 neighbour
 *   bits vs 64 for uint8_t — 8× better cache density.  fastEdge() replaces
 *   edge() and uses a bit-shift instead of an integer multiply.
 *
 * [OPT-2] BITSET colorBound
 *   Each colour class is represented as a compact bitset (NWORDS words).
 *   The conflict check — "is v adjacent to any vertex already in colour c?" —
 *   is now a word-level AND of adjRow(v) and the colour-class bitset, running
 *   in O(NWORDS) instead of O(i) (a full scan of prior vertices).
 *   This turns the overall bound from O(n²) to O(n × C × NWORDS) where
 *   C ≤ n is the number of colour classes and NWORDS = ⌈N/64⌉.
 *
 * [OPT-3] STATIC COLOUR BUFFER
 *   colorBound() is called once per node expansion.  Instead of allocating a
 *   fresh vector<uint64_t> on every call, a single process-wide buffer
 *   (g_color_buf) is grown on demand and reused, eliminating heap pressure.
 *
 * [OPT-4] ZERO-ALLOCATION KNAPSACK BOUND (pre-existing)
 *   global_eff_order + in_cand[] lookup avoid re-sorting candidates inside
 *   knapsackBound.  Unchanged from the original.
 *
 * Compile : mpicxx -O3 -std=c++17 -o clique_mpi main.cpp
 * Run     : mpirun -n <P> ./clique_mpi input.txt output.txt
 */

#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
//  Message tags
// ═══════════════════════════════════════════════════════════════════════════
enum Tag {
    TAG_STEAL_REQ  = 1,   // Idle rank asking for work
    TAG_WORK       = 2,   // Busy rank sending half its stack
    TAG_NACK       = 3,   // Busy rank has no work to spare
    TAG_NEW_BEST   = 4,   // Someone found a new P_max
    TAG_TOKEN      = 5,   // Dijkstra termination token
    TAG_TERMINATE  = 6    // System is finished
};

// Colors for Dijkstra's Token algorithm
enum Color { WHITE = 0, BLACK = 1 };

// ═══════════════════════════════════════════════════════════════════════════
//  Graph & Global Lookups
// ═══════════════════════════════════════════════════════════════════════════
static int           N, E, B;
static vector<int>   profit, cost;

// ── [OPT-1] BITSET ADJACENCY ───────────────────────────────────────────────
// Packed as N rows × NWORDS uint64_t words (NWORDS = ⌈N/64⌉).
// Adjacency bit for edge (u,v) sits at adj_bits[u*NWORDS + (v>>6)], bit (v&63).
// One 64-byte cache line holds 8 words = 512 neighbour bits (vs 64 for uint8_t).
static int             NWORDS = 0;
static vector<uint64_t> adj_bits;

inline const uint64_t* adjRow(int u) {
    return adj_bits.data() + (size_t)u * NWORDS;
}
inline bool fastEdge(int u, int v) {
    return (adjRow(u)[v >> 6] >> (v & 63)) & 1ULL;
}
// ──────────────────────────────────────────────────────────────────────────

static int         g_pmax = 0;
static vector<int> g_best;

// --- [OPT-4] Global Efficiency Lookup Arrays (pre-existing) ---
static vector<int>     global_eff_order;
static vector<uint8_t> in_cand;

// ── [OPT-3] STATIC COLOUR BUFFER ──────────────────────────────────────────
// Reused across all colorBound() calls to avoid repeated heap allocation.
// Grows on demand, never shrinks (monotone).
static vector<uint64_t> g_color_buf;
// ──────────────────────────────────────────────────────────────────────────

// ═══════════════════════════════════════════════════════════════════════════
//  Subproblem + serialisation
// ═══════════════════════════════════════════════════════════════════════════
struct Sub {
    vector<int> cand;
    vector<int> clique;
    int P = 0;
    int W = 0;
};

static vector<int> packSub(const Sub& s) {
    vector<int> b;
    b.reserve(4 + s.cand.size() + s.clique.size());
    b.push_back((int)s.cand.size());
    b.insert(b.end(), s.cand.begin(), s.cand.end());
    b.push_back((int)s.clique.size());
    b.insert(b.end(), s.clique.begin(), s.clique.end());
    b.push_back(s.P);
    b.push_back(s.W);
    return b;
}

static Sub unpackSub(const int* d) {
    Sub s; int i = 0;
    int cs = d[i++]; s.cand.assign(d + i, d + i + cs);   i += cs;
    int cl = d[i++]; s.clique.assign(d + i, d + i + cl); i += cl;
    s.P = d[i++];
    s.W = d[i++];
    return s;
}

static vector<int> packBatch(const vector<Sub>& subs) {
    vector<int> b;
    b.push_back((int)subs.size());
    for (const Sub& s : subs) {
        auto p = packSub(s);
        b.push_back((int)p.size());
        b.insert(b.end(), p.begin(), p.end());
    }
    return b;
}

static vector<Sub> unpackBatch(const int* d) {
    vector<Sub> out;
    int idx = 0, n = d[idx++];
    for (int i = 0; i < n; i++) {
        int sz = d[idx++];
        out.push_back(unpackSub(d + idx));
        idx += sz;
    }
    return out;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Bounding functions
// ═══════════════════════════════════════════════════════════════════════════

// ── [OPT-2] BITSET colorBound ─────────────────────────────────────────────
// Colour-class conflict detection is now O(NWORDS) per class per vertex:
//   conflict = any_set_bit( adjRow(v)  &  colour_class_bitset )
// Previously this required scanning all j < i (O(i) per class), giving O(n²)
// overall.  Now the overall work is O(n × C × NWORDS) where C ≤ n.
//
// The [OPT-3] static buffer g_color_buf is used here to avoid per-call
// heap allocation.  Its layout is: g_color_buf[c * NWORDS .. (c+1)*NWORDS-1]
// = bitset for colour class c.
//
// Algorithmic logic is IDENTICAL to the original: we assign each vertex to
// the first feasible colour class and add profit[v] to the upper bound only
// when a new colour class is opened.
static int colorBound(const vector<int>& cand) {
    int n = (int)cand.size();
    if (!n) return 0;

    int num_colors = 0;
    int ub = 0;

    for (int i = 0; i < n; i++) {
        int v       = cand[i];
        const uint64_t* vrow = adjRow(v);   // [OPT-1] bitset row for v
        bool placed = false;

        for (int c = 0; c < num_colors; c++) {
            // [OPT-2] O(NWORDS) conflict check via bitwise AND
            const uint64_t* cs = g_color_buf.data() + c * NWORDS;
            bool conflict = false;
            for (int w = 0; w < NWORDS; w++) {
                if (vrow[w] & cs[w]) { conflict = true; break; }
            }
            if (!conflict) {
                // Place v in colour class c — set its bit
                g_color_buf[c * NWORDS + (v >> 6)] |= 1ULL << (v & 63);
                placed = true;
                break;
            }
        }

        if (!placed) {
            // Open a new colour class for v
            // [OPT-3] grow static buffer on demand, zero-initialise new slot
            size_t new_size = (size_t)(num_colors + 1) * NWORDS;
            if (g_color_buf.size() < new_size)
                g_color_buf.resize(new_size, 0ULL);
            // Zero-out the new slot (previous uses may have left stale bits)
            fill(g_color_buf.begin() + num_colors * NWORDS,
                 g_color_buf.begin() + new_size, 0ULL);
            g_color_buf[num_colors * NWORDS + (v >> 6)] |= 1ULL << (v & 63);
            num_colors++;
            ub += profit[v];
        }
    }
    return ub;
}

// ── [OPT-4] Zero-Allocation O(N) Knapsack Bound (pre-existing, unchanged) ─
static int knapsackBound(const vector<int>& cand, int rem) {
    if (rem <= 0 || cand.empty()) return 0;

    for (int v : cand) in_cand[v] = 1;

    int ub = 0, left = rem;
    int found = 0;
    int target = (int)cand.size();

    for (int v : global_eff_order) {
        if (!in_cand[v]) continue;

        if (cost[v] <= left) {
            ub += profit[v];
            left -= cost[v];
        } else {
            ub += (profit[v] * left) / cost[v];
            break;
        }

        if (++found == target) break;
    }

    for (int v : cand) in_cand[v] = 0;

    return ub;
}

// ═══════════════════════════════════════════════════════════════════════════
//  expand()
// ═══════════════════════════════════════════════════════════════════════════
static void expand(const Sub& sp, deque<Sub>& stk) {
    if (sp.P + colorBound(sp.cand) <= g_pmax) return;
    if (sp.P + knapsackBound(sp.cand, B - sp.W) <= g_pmax) return;

    int n = (int)sp.cand.size();
    for (int i = n - 1; i >= 0; i--) {
        int v = sp.cand[i];
        if (sp.W + cost[v] > B) continue;

        int np = sp.P + profit[v];
        if (np > g_pmax) {
            g_pmax = np;
            g_best = sp.clique;
            g_best.push_back(v);
        }
        if (sp.W + cost[v] == B) continue;

        Sub child;
        child.P = np;
        child.W = sp.W + cost[v];
        child.clique = sp.clique;
        child.clique.push_back(v);
        child.cand.reserve(n - i - 1);

        // [OPT-1] Use the bitset adjacency row for v — avoids integer multiply
        // and gives better cache behaviour for large N.  One cache line now
        // covers 512 bits of neighbour information.
        const uint64_t* vrow = adjRow(v);
        for (int j = i + 1; j < n; j++) {
            int u = sp.cand[j];
            if ((vrow[u >> 6] >> (u & 63)) & 1ULL)
                child.cand.push_back(u);
        }

        stk.push_back(move(child));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  rootDecompose() — used for initial Static Scatter
// ═══════════════════════════════════════════════════════════════════════════
static deque<Sub> rootDecompose(int target) {
    deque<Sub> q;
    Sub root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });
    q.push_back(root);

    while ((int)q.size() < target && !q.empty()) {
        auto it = max_element(q.begin(), q.end(), [](const Sub& a, const Sub& b){
            return a.cand.size() < b.cand.size();
        });
        Sub sp = *it; q.erase(it);
        if (sp.cand.empty()) { q.push_back(sp); break; }

        if (sp.P + colorBound(sp.cand) <= g_pmax) continue;
        if (sp.P + knapsackBound(sp.cand, B - sp.W) <= g_pmax) continue;

        int v = sp.cand[0];

        // Branch 1 : WITHOUT v
        Sub without = sp;
        without.cand.erase(without.cand.begin());
        if (!without.cand.empty()) q.push_back(without);

        // Branch 2 : WITH v
        if (sp.W + cost[v] <= B) {
            Sub with_v;
            with_v.P = sp.P + profit[v];
            with_v.W = sp.W + cost[v];
            with_v.clique = sp.clique;
            with_v.clique.push_back(v);
            // [OPT-1] bitset row for v used in candidate filtering
            const uint64_t* vrow = adjRow(v);
            for (size_t i = 1; i < sp.cand.size(); i++) {
                int u = sp.cand[i];
                if ((vrow[u >> 6] >> (u & 63)) & 1ULL)
                    with_v.cand.push_back(u);
            }
            if (with_v.P > g_pmax) { g_pmax = with_v.P; g_best = with_v.clique; }
            q.push_back(with_v);
        }
    }
    return q;
}

// ═══════════════════════════════════════════════════════════════════════════
//  Decentralized Compute Engine
// ═══════════════════════════════════════════════════════════════════════════
static void runDecentralized(int rank, int P, const char* out_filename) {
    vector<int> rbuf;
    deque<Sub> local;

    // --- SAFE ASYNC MEMORY: Persistent buffers and requests for MPI_Isend ---
    vector<vector<int>> init_sbufs(P);
    vector<MPI_Request> init_reqs(P, MPI_REQUEST_NULL);

    // --- State Variables ---
    bool is_idle = false;
    bool done = false;
    bool steal_pending = false;

    // --- Dijkstra Token Variables ---
    int  my_color = WHITE;
    bool has_token = (rank == 0);
    int  token_color = WHITE;

    int poll_counter = 0;
    const int POLL_FREQ = 1024;

    srand(time(NULL) + rank);

    // ── 1. THE DECOMPOSITION PHASE ─────────────────────────────────────────
    deque<Sub> initial_chunks;
    if (rank == 0) {
        const int CHUNK_MULTIPLIER = 4;
        initial_chunks = rootDecompose(P * CHUNK_MULTIPLIER);
    }

    MPI_Bcast(&g_pmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // ── 2. THE KICKSTART (Asynchronous Non-Blocking Scatter) ───────────────
    if (rank == 0) {
        if (!initial_chunks.empty()) {
            local.push_back(initial_chunks.front());
            initial_chunks.pop_front();
        }

        int target_rank = 1;
        while (!initial_chunks.empty() && target_rank < P) {
            vector<Sub> init_batch;
            init_batch.push_back(initial_chunks.front());
            initial_chunks.pop_front();

            init_sbufs[target_rank] = packBatch(init_batch);
            MPI_Isend(init_sbufs[target_rank].data(), (int)init_sbufs[target_rank].size(),
                      MPI_INT, target_rank, TAG_WORK, MPI_COMM_WORLD, &init_reqs[target_rank]);
            target_rank++;
        }

        while (target_rank < P) {
            init_sbufs[target_rank] = {0};
            MPI_Isend(init_sbufs[target_rank].data(), 1, MPI_INT, target_rank,
                      TAG_NACK, MPI_COMM_WORLD, &init_reqs[target_rank]);
            target_rank++;
        }

        while (!initial_chunks.empty()) {
            local.push_back(initial_chunks.front());
            initial_chunks.pop_front();
        }

        if (local.empty()) is_idle = true;
    }
    else {
        MPI_Status st;
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
        if ((size_t)cnt > rbuf.size()) rbuf.resize((size_t)cnt);

        MPI_Recv(rbuf.data(), cnt, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

        if (st.MPI_TAG == TAG_WORK) {
            auto subs = unpackBatch(rbuf.data());
            for (auto& s : subs) local.push_back(move(s));
        } else {
            is_idle = true;
        }
    }

    if (rank == 0) {
        MPI_Waitall(P - 1, init_reqs.data() + 1, MPI_STATUSES_IGNORE);
    }

    auto passToken = [&]() {
        if (has_token) {
            int out_token;
            if (!is_idle) {
                out_token = BLACK;
            } else {
                out_token = (token_color == BLACK || my_color == BLACK) ? BLACK : WHITE;
            }
            int next_rank = (rank + 1) % P;
            MPI_Send(&out_token, 1, MPI_INT, next_rank, TAG_TOKEN, MPI_COMM_WORLD);
            my_color = WHITE;
            has_token = false;
        }
    };

    while (!done) {
        // ── 3. Throttled Message Poll ──────────────────────────────────────
        if (is_idle || (++poll_counter % POLL_FREQ == 0)) {
            for (;;) {
                int flag = 0;
                MPI_Status st;
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
                if (!flag) break;

                int src = st.MPI_SOURCE;
                int tag = st.MPI_TAG;

                int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
                if ((size_t)cnt > rbuf.size()) rbuf.resize((size_t)cnt);

                MPI_Recv(rbuf.data(), cnt, MPI_INT, src, tag, MPI_COMM_WORLD, &st);

                switch (tag) {
                    case TAG_NEW_BEST: {
                        int np = rbuf[0];
                        if (np > g_pmax) g_pmax = np;
                        break;
                    }

                    case TAG_STEAL_REQ: {
                        if (local.size() >= 2) {
                            int half = (int)local.size() / 2;
                            vector<Sub> donated;
                            for (int i = 0; i < half; i++) {
                                donated.push_back(local.front());
                                local.pop_front();
                            }
                            auto sbuf = packBatch(donated);
                            MPI_Send(sbuf.data(), (int)sbuf.size(), MPI_INT, src, TAG_WORK, MPI_COMM_WORLD);
                            my_color = BLACK;
                        } else {
                            int d = 0;
                            MPI_Send(&d, 1, MPI_INT, src, TAG_NACK, MPI_COMM_WORLD);
                        }
                        break;
                    }

                    case TAG_WORK: {
                        auto subs = unpackBatch(rbuf.data());
                        for (auto& s : subs) local.push_back(move(s));
                        steal_pending = false;
                        is_idle = false;
                        break;
                    }

                    case TAG_NACK: {
                        steal_pending = false;
                        break;
                    }

                    case TAG_TOKEN: {
                        if (rank == 0) {
                            int in_color = rbuf[0];
                            if (in_color == WHITE && is_idle) {
                                int d = 0;
                                for (int r = 1; r < P; r++) {
                                    MPI_Send(&d, 1, MPI_INT, r, TAG_TERMINATE, MPI_COMM_WORLD);
                                }
                                done = true;
                            } else {
                                has_token = true;
                                token_color = WHITE;
                            }
                        } else {
                            has_token = true;
                            token_color = rbuf[0];
                        }
                        break;
                    }

                    case TAG_TERMINATE: {
                        if (steal_pending) {
                            MPI_Status drain_st;
                            int dflag = 0;
                            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &dflag, &drain_st);
                            if (dflag) {
                                int dcnt; MPI_Get_count(&drain_st, MPI_INT, &dcnt);
                                if ((size_t)dcnt > rbuf.size()) rbuf.resize((size_t)dcnt);
                                MPI_Recv(rbuf.data(), dcnt, MPI_INT,
                                         drain_st.MPI_SOURCE, drain_st.MPI_TAG,
                                         MPI_COMM_WORLD, &drain_st);
                            }
                            steal_pending = false;
                        }
                        done = true;
                        break;
                    }
                }
            }
        }

        if (done) break;

        // ── 4. Token Ring (Pass immediately to prevent stalling) ───────────
        if (local.empty()) {
            is_idle = true;
        }
        passToken();

        // ── 5. Work Stealing Logic (If Idle) ───────────────────────────────
        if (is_idle) {
            if (!steal_pending && P > 1) {
                int target;
                do { target = rand() % P; } while (target == rank);
                int d = 0;
                MPI_Send(&d, 1, MPI_INT, target, TAG_STEAL_REQ, MPI_COMM_WORLD);
                steal_pending = true;
            }
            continue;
        }

        // ── 6. Compute (If Active) ─────────────────────────────────────────
        Sub sp = local.back();
        local.pop_back();

        int prev = g_pmax;
        expand(sp, local);

        if (g_pmax > prev) {
            for (int r = 0; r < P; r++) {
                if (r != rank) {
                    MPI_Send(&g_pmax, 1, MPI_INT, r, TAG_NEW_BEST, MPI_COMM_WORLD);
                }
            }
        }
    }

    // ── 7. Final Aggregation (Safe Collective Version) ─────────────────────
    int actual_local_profit = 0;
    for (int v : g_best) actual_local_profit += profit[v];

    struct { int val; int rank; } in_val, out_val;
    in_val.val = actual_local_profit;
    in_val.rank = rank;

    MPI_Allreduce(&in_val, &out_val, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

    int global_max_profit = out_val.val;
    int winner_rank = out_val.rank;

    if (rank == 0) {
        if (winner_rank != 0) {
            MPI_Status st;
            MPI_Probe(winner_rank, 999, MPI_COMM_WORLD, &st);
            int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
            if ((size_t)cnt > rbuf.size()) rbuf.resize((size_t)cnt);

            MPI_Recv(rbuf.data(), cnt, MPI_INT, winner_rank, 999, MPI_COMM_WORLD, &st);
            int csz = rbuf[1];
            g_best.clear();
            for (int i = 0; i < csz; i++) g_best.push_back(rbuf[2 + i]);
        }

        sort(g_best.begin(), g_best.end());
        ofstream fout(out_filename);
        if (!fout.is_open()) {
            cerr << "Error: Could not open output file " << out_filename << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fout << global_max_profit << "\n";
        for (int i = 0; i < (int)g_best.size(); i++) {
            if (i) fout << " ";
            fout << g_best[i];
        }
        fout << "\n";
        fout.close();
    }
    else {
        if (rank == winner_rank) {
            vector<int> b;
            b.push_back(actual_local_profit);
            b.push_back((int)g_best.size());
            for (int v : g_best) b.push_back(v);
            MPI_Send(b.data(), (int)b.size(), MPI_INT, 0, 999, MPI_COMM_WORLD);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Sequential fallback (P == 1)
// ═══════════════════════════════════════════════════════════════════════════
static void runSequential(const char* out_filename) {
    deque<Sub> stk;
    Sub root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });
    stk.push_back(root);

    while (!stk.empty()) {
        Sub sp = stk.back(); stk.pop_back();
        expand(sp, stk);
    }

    sort(g_best.begin(), g_best.end());
    ofstream fout(out_filename);
    if (!fout.is_open()) {
        cerr << "Error: Could not open output file " << out_filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fout << g_pmax << "\n";
    for (int i = 0; i < (int)g_best.size(); i++) {
        if (i) fout << " ";
        fout << g_best[i];
    }
    fout << "\n";
    fout.close();
}

// ═══════════════════════════════════════════════════════════════════════════
//  main
// ═══════════════════════════════════════════════════════════════════════════
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (argc != 3) {
        if (rank == 0) cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        ifstream fin(argv[1]);
        if (!fin.is_open()) {
            cerr << "Error: Could not open input file " << argv[1] << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fin >> N >> E >> B;
        profit.resize(N); cost.resize(N);
        for (int i = 0; i < N; i++) fin >> profit[i] >> cost[i];

        // [OPT-1] Build bitset adjacency instead of uint8_t flat array
        NWORDS = (N + 63) / 64;
        adj_bits.assign((size_t)N * NWORDS, 0ULL);
        for (int i = 0; i < E; i++) {
            int u, v; fin >> u >> v;
            adj_bits[(size_t)u * NWORDS + (v >> 6)] |= 1ULL << (v & 63);
            adj_bits[(size_t)v * NWORDS + (u >> 6)] |= 1ULL << (u & 63);
        }
        fin.close();
    }

    {
        int neb[3] = { N, E, B };
        MPI_Bcast(neb, 3, MPI_INT, 0, MPI_COMM_WORLD);
        N = neb[0]; E = neb[1]; B = neb[2];
    }

    // [OPT-1] Broadcast NWORDS then the bitset adjacency matrix
    NWORDS = (N + 63) / 64;
    profit.resize(N); cost.resize(N);
    MPI_Bcast(profit.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cost.data(),   N, MPI_INT, 0, MPI_COMM_WORLD);
    adj_bits.resize((size_t)N * NWORDS, 0ULL);
    // Broadcast as bytes — MPI_UINT64_T is in MPI-3 but casting to byte is portable
    MPI_Bcast(adj_bits.data(), (int)(adj_bits.size() * sizeof(uint64_t)),
              MPI_BYTE, 0, MPI_COMM_WORLD);

    // [OPT-4] Setup global knapsack efficiency order (pre-existing)
    in_cand.assign(N, 0);
    global_eff_order.resize(N);
    iota(global_eff_order.begin(), global_eff_order.end(), 0);
    sort(global_eff_order.begin(), global_eff_order.end(), [&](int a, int b) {
        long long numA = profit[a], denA = cost[a];
        long long numB = profit[b], denB = cost[b];
        if (numA * denB != numB * denA) return numA * denB > numB * denA;
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });

    // [OPT-3] Pre-allocate the static colour buffer with a reasonable
    // starting capacity (avoids first-call reallocation for small graphs).
    g_color_buf.reserve((size_t)N * NWORDS);

    // Initial greedy pass for a strong starting P_max
    for (int v = 0; v < N; v++) {
        if (cost[v] <= B && profit[v] > g_pmax) {
            g_pmax = profit[v];
            g_best = { v };
        }
    }

    if (P == 1) {
        runSequential(argv[2]);
    } else {
        runDecentralized(rank, P, argv[2]);
    }

    MPI_Finalize();
    return 0;
}