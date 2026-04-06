/*
 * MPI Parallel Branch-and-Bound — Budgeted Maximum Clique Problem
 *
 * Architecture
 * ─────────────
 *  Rank 0  : Master coordinator
 *              • reads input, broadcasts graph
 *              • root-decomposes search space → initial work queue
 *              • dispatches subproblems to workers dynamically
 *              • collects TAG_NEW_BEST → rebroadcasts tighter P_max
 *              • detects global termination, prints final answer
 *
 *  Ranks 1…P-1 : Workers
 *              • maintain a local LIFO stack of subproblems
 *              • run the exact BnB kernel (colour bound + knapsack bound)
 *              • non-blocking message poll at every iteration:
 *                  – TAG_NEW_BEST  : update local P_max (tighter pruning)
 *                  – TAG_STEAL_REQ : donate half local stack if large enough
 *                  – TAG_STEAL_RESP/NACK : resolve pending steal
 *                  – TAG_WORK / TAG_NO_MORE : master response
 *              • when local stack empties:
 *                  1. one work-steal attempt from a peer
 *                  2. if steal fails → request work from master
 *                  3. if master has none → send TAG_RESULT, terminate
 *
 * Compile : mpicxx -O2 -std=c++17 -o clique_mpi clique_mpi.cpp
 * Run     : mpirun -n <P> ./clique_mpi  input.txt output.txt
 */

#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
//  Message tags
// ═══════════════════════════════════════════════════════════════════════════
enum Tag {
    TAG_WORK       = 1,   // master → worker  : here is a subproblem
    TAG_WORK_REQ   = 2,   // worker → master  : give me more work
    TAG_NO_MORE    = 3,   // master → worker  : no work left, terminate
    TAG_NEW_BEST   = 4,   // worker ↔ master  : updated P_max
    TAG_RESULT     = 5,   // worker → master  : final local best
    TAG_STEAL_REQ  = 6,   // worker → worker  : can I steal half your stack?
    TAG_STEAL_RESP = 7,   // worker → worker  : stolen subproblems
    TAG_STEAL_NACK = 8,   // worker → worker  : nothing to steal
};

// ═══════════════════════════════════════════════════════════════════════════
//  Graph — replicated on every rank after MPI_Bcast
// ═══════════════════════════════════════════════════════════════════════════
static int           N, E;
static long long     B;
static vector<int>   profit, cost;
static vector<uint8_t> adj;           // flat N×N adjacency matrix (1 byte/cell)
inline bool edge(int u, int v) { return adj[u * N + v]; }

// ═══════════════════════════════════════════════════════════════════════════
//  Per-rank best solution
// ═══════════════════════════════════════════════════════════════════════════
static long long   g_pmax = 0;
static vector<int> g_best;

// ═══════════════════════════════════════════════════════════════════════════
//  Subproblem + serialisation
//  Packed as a flat vector<int>; long longs split into (lo, hi) int pairs.
// ═══════════════════════════════════════════════════════════════════════════
struct Sub {
    vector<int> cand;     // candidate vertices (all pairwise adjacent to clique)
    vector<int> clique;   // vertices already chosen
    long long   P = 0;    // profit of current partial clique
    long long   W = 0;    // cost   of current partial clique
};

static void pushLL(vector<int>& b, long long v) {
    b.push_back((int)((uint64_t)v & 0xFFFFFFFFu));
    b.push_back((int)((uint64_t)v >> 32));
}
static long long pullLL(const int* b, int& i) {
    uint64_t lo = (uint32_t)b[i++];
    uint64_t hi = (uint32_t)b[i++];
    return (long long)(lo | (hi << 32));
}

static vector<int> packSub(const Sub& s) {
    vector<int> b;
    b.reserve(4 + s.cand.size() + s.clique.size());
    b.push_back((int)s.cand.size());
    b.insert(b.end(), s.cand.begin(), s.cand.end());
    b.push_back((int)s.clique.size());
    b.insert(b.end(), s.clique.begin(), s.clique.end());
    pushLL(b, s.P);
    pushLL(b, s.W);
    return b;
}
static Sub unpackSub(const int* d) {
    Sub s; int i = 0;
    int cs = d[i++]; s.cand.assign(d + i, d + i + cs);   i += cs;
    int cl = d[i++]; s.clique.assign(d + i, d + i + cl); i += cl;
    s.P = pullLL(d, i);
    s.W = pullLL(d, i);
    return s;
}

// Pack a batch of subproblems: [count, sz0, data0..., sz1, data1..., ...]
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
//  P_max helpers
// ═══════════════════════════════════════════════════════════════════════════
static void sendPmax(int dest, long long p) {
    int b[2] = { (int)((uint64_t)p & 0xFFFFFFFFu), (int)((uint64_t)p >> 32) };
    MPI_Send(b, 2, MPI_INT, dest, TAG_NEW_BEST, MPI_COMM_WORLD);
}
static long long extractPmax(const int* b) {
    return (long long)((uint64_t)(uint32_t)b[0] | ((uint64_t)(uint32_t)b[1] << 32));
}

// ═══════════════════════════════════════════════════════════════════════════
//  Bounding functions   (identical to sequential — per-rank, no MPI)
// ═══════════════════════════════════════════════════════════════════════════

// Structural bound: greedy graph colouring.
// Upper bound = sum of max profit per colour class (at most 1 vertex per class).
static long long colorBound(const vector<int>& cand) {
    int n = (int)cand.size();
    if (!n) return 0;
    vector<int>      col(n, -1);
    vector<long long> cbest;
    for (int i = 0; i < n; i++) {
        int v = cand[i];
        bool placed = false;
        for (int c = 0; c < (int)cbest.size(); c++) {
            bool conf = false;
            for (int j = 0; j < i; j++)
                if (col[j] == c && edge(v, cand[j])) { conf = true; break; }
            if (!conf) {
                col[i] = c;
                cbest[c] = max(cbest[c], (long long)profit[v]);
                placed = true; break;
            }
        }
        if (!placed) { col[i] = (int)cbest.size(); cbest.push_back(profit[v]); }
    }
    long long ub = 0;
    for (long long x : cbest) ub += x;
    return ub;
}

// Resource bound: fractional knapsack on remaining budget.
static long long knapsackBound(const vector<int>& cand, long long rem) {
    if (rem <= 0 || cand.empty()) return 0;
    vector<int> ord(cand.size()); iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](int a, int b) {
        return (long long)profit[cand[a]] * cost[cand[b]] >
               (long long)profit[cand[b]] * cost[cand[a]];
    });
    long long ub = 0, left = rem;
    for (int idx : ord) {
        int v = cand[idx];
        if (cost[v] <= left) { ub += profit[v]; left -= cost[v]; }
        else { ub += ((long long)profit[v] * left + cost[v] - 1) / cost[v]; break; }
    }
    return ub;
}

// ═══════════════════════════════════════════════════════════════════════════
//  expand() — process one BnB node, push children onto local_stack.
//
//  Sorting trick: sort cand ascending by profit so the highest-profit child
//  is pushed LAST and therefore sits on TOP of the LIFO stack — the best
//  branch is always explored first, tightening P_max early.
//
//  C_next construction: for vertex rem[i], C_next = rem[i+1..] ∩ Nbrs(v).
//  This mirrors the sequential algorithm's shrinking-cand technique and
//  ensures each combination is explored exactly once.
//
//  Returns true if a new local best was discovered.
// ═══════════════════════════════════════════════════════════════════════════
static bool expand(const Sub& sp, stack<Sub>& stk) {
    // ── Prune with both bounds ────────────────────────────────────────────
    if (sp.P + colorBound(sp.cand)               <= g_pmax) return false;
    if (sp.P + knapsackBound(sp.cand, B - sp.W)  <= g_pmax) return false;

    // Sort ascending by profit → best child pushed last → top of stack
    vector<int> rem = sp.cand;
    sort(rem.begin(), rem.end(), [](int a, int b){ return profit[a] < profit[b]; });

    bool found = false;
    for (int i = 0; i < (int)rem.size(); i++) {
        int v = rem[i];
        if (sp.W + cost[v] > B) continue;

        long long np = sp.P + profit[v];
        if (np > g_pmax) {
            g_pmax = np;
            g_best = sp.clique;
            g_best.push_back(v);
            found = true;
        }

        // Build C_next = rem[i+1..] ∩ Neighbours(v)
        Sub child;
        child.P = np;
        child.W = sp.W + cost[v];
        child.clique = sp.clique;
        child.clique.push_back(v);
        for (int j = i + 1; j < (int)rem.size(); j++)
            if (edge(v, rem[j])) child.cand.push_back(rem[j]);
        stk.push(move(child));
    }
    return found;
}

// ═══════════════════════════════════════════════════════════════════════════
//  rootDecompose() — rank 0 only.
//
//  BFS-style expansion of the root node until we have ≥ target subproblems
//  in the queue, or the graph is exhausted.  Each expansion step picks the
//  subproblem with the largest candidate set (most work) and splits it into:
//    • "with v" : clique ∪ {v},  cand ← cand ∩ Nbrs(v)
//    • "without v" : clique unchanged, cand ← cand \ {v}
//  where v is the highest-profit remaining candidate (best first).
//
//  Both bounds are applied during decomposition to prune unpromising nodes
//  immediately, so workers never receive obviously dead subproblems.
// ═══════════════════════════════════════════════════════════════════════════
static deque<Sub> rootDecompose(int target) {
    deque<Sub> q;
    Sub root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    // Sort descending by profit so the first split peels off the best vertex
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        return profit[a] > profit[b];
    });
    q.push_back(root);

    while ((int)q.size() < target && !q.empty()) {
        // Pick the largest candidate set to expand next
        auto it = max_element(q.begin(), q.end(), [](const Sub& a, const Sub& b){
            return a.cand.size() < b.cand.size();
        });
        Sub sp = *it; q.erase(it);
        if (sp.cand.empty()) { q.push_back(sp); break; }   // nothing to split

        // Apply both bounds; prune if hopeless
        if (sp.P + colorBound(sp.cand)              <= g_pmax) continue;
        if (sp.P + knapsackBound(sp.cand, B - sp.W) <= g_pmax) continue;

        // Split on the highest-profit candidate (last in ascending sort = front
        // in descending sort; our cand is already sorted descending from init)
        vector<int> rem = sp.cand;
        sort(rem.begin(), rem.end(), [](int a, int b){ return profit[a] < profit[b]; });
        int v = rem.back(); rem.pop_back();   // peel highest-profit vertex

        // ── Branch 1 : WITHOUT v (explore rest of candidates) ────────────
        if (!rem.empty()) {
            Sub without = sp;
            without.cand = rem;
            q.push_back(without);
        }

        // ── Branch 2 : WITH v (add v to clique, restrict candidates) ─────
        if (sp.W + cost[v] <= B) {
            Sub with_v;
            with_v.P = sp.P + profit[v];
            with_v.W = sp.W + cost[v];
            with_v.clique = sp.clique;
            with_v.clique.push_back(v);
            for (int u : rem) if (edge(v, u)) with_v.cand.push_back(u);
            // Update g_pmax/g_best with this partial clique
            if (with_v.P > g_pmax) { g_pmax = with_v.P; g_best = with_v.clique; }
            q.push_back(with_v);
        }
    }

    return q;
}

// ═══════════════════════════════════════════════════════════════════════════
//  MASTER — rank 0
// ═══════════════════════════════════════════════════════════════════════════
static void runMaster(int P, const char* out_filename) {
    const int W = P - 1;                    // number of worker ranks

    // ── Generate initial work pool ────────────────────────────────────────
    deque<Sub> Q = rootDecompose(W * 4);

    long long   global_best_p = g_pmax;
    vector<int> global_best_c = g_best;

    vector<bool> is_idle(P, false);
    vector<bool> terminated(P, false);
    terminated[0] = true;                   // master itself
    int idle_count = 0;
    int done_count = 0;

    // ── Initial dispatch ──────────────────────────────────────────────────
    for (int r = 1; r <= W; r++) {
        if (!Q.empty()) {
            auto b = packSub(Q.front()); Q.pop_front();
            MPI_Send(b.data(), (int)b.size(), MPI_INT, r, TAG_WORK, MPI_COMM_WORLD);
        } else {
            int d = 0;
            MPI_Send(&d, 1, MPI_INT, r, TAG_NO_MORE, MPI_COMM_WORLD);
            terminated[r] = true;
            done_count++;
        }
    }

    // ── Coordination loop ─────────────────────────────────────────────────
    const int RBUF_SZ = 1 << 24;
    vector<int> rbuf(RBUF_SZ);

    while (done_count < W) {
        MPI_Status st;
        MPI_Recv(rbuf.data(), RBUF_SZ, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        int src = st.MPI_SOURCE;
        int tag = st.MPI_TAG;
        int cnt; MPI_Get_count(&st, MPI_INT, &cnt);

        if (tag == TAG_WORK_REQ) {
            if (!Q.empty()) {
                auto b = packSub(Q.front()); Q.pop_front();
                MPI_Send(b.data(), (int)b.size(), MPI_INT, src, TAG_WORK, MPI_COMM_WORLD);
            } else {
                is_idle[src] = true;
                idle_count++;

                int active = W - done_count;
                if (idle_count == active) {
                    for (int r = 1; r <= W; r++) {
                        if (is_idle[r] && !terminated[r]) {
                            int d = 0;
                            MPI_Send(&d, 1, MPI_INT, r, TAG_NO_MORE, MPI_COMM_WORLD);
                        }
                    }
                }
            }
        }
        else if (tag == TAG_NEW_BEST) {
            long long np = extractPmax(rbuf.data());
            if (np > global_best_p) {
                global_best_p = np;
                for (int r = 1; r <= W; r++)
                    if (!terminated[r] && !is_idle[r] && r != src)
                        sendPmax(r, global_best_p);
            }
        }
        else if (tag == TAG_RESULT) {
            long long lp = extractPmax(rbuf.data());
            if (lp > global_best_p) {
                global_best_p = lp;
                global_best_c.clear();
                int csz = rbuf[2];
                for (int i = 0; i < csz; i++) global_best_c.push_back(rbuf[3 + i]);
            }
            terminated[src] = true;
            done_count++;
        }
    }

    // ── Write answer to output file ───────────────────────────────────────
    sort(global_best_c.begin(), global_best_c.end());
    
    ofstream fout(out_filename);
    if (!fout.is_open()) {
        cerr << "Error: Could not open output file " << out_filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    fout << global_best_p << "\n";
    for (int i = 0; i < (int)global_best_c.size(); i++) {
        if (i) fout << " ";
        fout << global_best_c[i];
    }
    fout << "\n";
    fout.close();
}
// ═══════════════════════════════════════════════════════════════════════════
//  WORKER — ranks 1 … P-1
// ═══════════════════════════════════════════════════════════════════════════
static void runWorker(int rank, int P) {
    const int W = P - 1;
    const int RBUF_SZ = 1 << 24;
    vector<int> rbuf(RBUF_SZ);

    stack<Sub> local;
    bool wait_master   = false;  // waiting for TAG_WORK / TAG_NO_MORE from master
    bool steal_pending = false;  // waiting for TAG_STEAL_RESP / TAG_STEAL_NACK
    int  steal_target  = -1;     // which rank we sent STEAL_REQ to
    bool done          = false;

    // ── sendResult: report local best to master, then worker is finished ──
    auto sendResult = [&]() {
        vector<int> b;
        pushLL(b, g_pmax);
        b.push_back((int)g_best.size());
        for (int v : g_best) b.push_back(v);
        MPI_Send(b.data(), (int)b.size(), MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
    };

    // ── Receive initial assignment from master ────────────────────────────
    // Use MPI_ANY_TAG: master might send TAG_WORK or TAG_NO_MORE if queue empty.
    {
        MPI_Status st;
        MPI_Recv(rbuf.data(), RBUF_SZ, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        if (st.MPI_TAG == TAG_WORK)
            local.push(unpackSub(rbuf.data()));
        else
            done = true;          // TAG_NO_MORE — nothing to do
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  Main BnB loop
    // ═══════════════════════════════════════════════════════════════════════
    while (!done) {

        // ── Non-blocking message poll ─────────────────────────────────────
        // Drain ALL pending messages before processing one subproblem.
        // This ensures TAG_NEW_BEST updates arrive as early as possible,
        // maximising pruning effectiveness.
        for (;;) {
            int flag = 0;
            MPI_Status st;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
            if (!flag) break;

            int src = st.MPI_SOURCE;
            int tag = st.MPI_TAG;
            int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
            MPI_Recv(rbuf.data(), cnt, MPI_INT, src, tag, MPI_COMM_WORLD, &st);

            // ── Updated global best → sharpen local pruning bound ─────────
            if (tag == TAG_NEW_BEST) {
                long long np = extractPmax(rbuf.data());
                if (np > g_pmax) g_pmax = np;
            }

            // ── Master granted work we requested ─────────────────────────
            else if (tag == TAG_WORK && wait_master) {
                local.push(unpackSub(rbuf.data()));
                wait_master = false;
            }

            // ── Master says no more work → terminate ──────────────────────
            else if (tag == TAG_NO_MORE && wait_master) {
                done = true;
                wait_master = false;
            }

            // ── Another worker wants to steal half our local stack ─────────
            // We respond even while steal_pending (our stack might still have
            // items if the steal came in before our own stack emptied).
            else if (tag == TAG_STEAL_REQ) {
                const int STEAL_THRESHOLD = 4;
                if ((int)local.size() > STEAL_THRESHOLD && !wait_master) {
                    int half = (int)local.size() / 2;
                    vector<Sub> donated;
                    donated.reserve(half);
                    for (int i = 0; i < half; i++) {
                        donated.push_back(local.top());
                        local.pop();
                    }
                    auto sbuf = packBatch(donated);
                    MPI_Send(sbuf.data(), (int)sbuf.size(), MPI_INT,
                             src, TAG_STEAL_RESP, MPI_COMM_WORLD);
                } else {
                    int d = 0;
                    MPI_Send(&d, 1, MPI_INT, src, TAG_STEAL_NACK, MPI_COMM_WORLD);
                }
            }

            // ── Victim gave us stolen subproblems ─────────────────────────
            // Accept regardless of steal_pending flag to handle any delayed
            // response that arrives while we've moved on.
            else if (tag == TAG_STEAL_RESP) {
                auto subs = unpackBatch(rbuf.data());
                for (auto& s : subs) local.push(move(s));
                if (steal_pending && src == steal_target)
                    steal_pending = false;
            }

            // ── Victim had nothing to give ────────────────────────────────
            // Fall back to asking master for work.
            else if (tag == TAG_STEAL_NACK && steal_pending && src == steal_target) {
                steal_pending = false;
                int d = 0;
                MPI_Send(&d, 1, MPI_INT, 0, TAG_WORK_REQ, MPI_COMM_WORLD);
                wait_master = true;
            }
        } // end poll loop

        if (done) break;
        // Waiting for a response — spin and keep polling rather than blocking.
        if (wait_master || steal_pending) continue;

        // ── Local stack is empty — try to get more work ───────────────────
        if (local.empty()) {
            if (W > 1) {
                // Ring steal: ask the next worker in rank order.
                // Workers are numbered 1..W; the ring wraps around.
                steal_target = (rank % W) + 1;          // next worker in ring
                if (steal_target == rank)                // edge case: W==1
                    steal_target = (rank == 1) ? 2 : 1;
                int d = rank;
                MPI_Send(&d, 1, MPI_INT, steal_target, TAG_STEAL_REQ, MPI_COMM_WORLD);
                steal_pending = true;
            } else {
                // Only one worker: go straight to master.
                int d = 0;
                MPI_Send(&d, 1, MPI_INT, 0, TAG_WORK_REQ, MPI_COMM_WORLD);
                wait_master = true;
            }
            continue;
        }

        // ── Expand one subproblem from the top of the local stack ─────────
        Sub sp = local.top(); local.pop();
        long long prev = g_pmax;
        expand(sp, local);

        // If we found a new best, tell master so it can rebroadcast.
        if (g_pmax > prev) sendPmax(0, g_pmax);
    }

    sendResult();
}

// ═══════════════════════════════════════════════════════════════════════════
//  Sequential fallback (P == 1)
// ═══════════════════════════════════════════════════════════════════════════
static void runSequential(const char* out_filename) {
    stack<Sub> stk;
    Sub root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    stk.push(root);
    
    while (!stk.empty()) {
        Sub sp = stk.top(); stk.pop();
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

    // Ensure proper argument count
    if (argc != 3) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        }
        MPI_Finalize();
        return 1;
    }

    // ── Rank 0 reads input from file ──────────────────────────────────────
    if (rank == 0) {
        ifstream fin(argv[1]);
        if (!fin.is_open()) {
            cerr << "Error: Could not open input file " << argv[1] << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1); // Abort all ranks if Rank 0 can't open file
        }

        fin >> N >> E >> B;
        profit.resize(N); cost.resize(N);
        for (int i = 0; i < N; i++) fin >> profit[i] >> cost[i];
        
        adj.assign((size_t)N * N, 0);
        for (int i = 0; i < E; i++) {
            int u, v; fin >> u >> v;
            adj[u * N + v] = adj[v * N + u] = 1;
        }
        fin.close();
    }

    // ── Broadcast graph to all ranks ──────────────────────────────────────
    {
        long long neb[3] = { N, E, B };
        MPI_Bcast(neb, 3, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        N = (int)neb[0]; E = (int)neb[1]; B = neb[2];
    }
    profit.resize(N); cost.resize(N);
    MPI_Bcast(profit.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cost.data(),   N, MPI_INT, 0, MPI_COMM_WORLD);
    adj.resize((size_t)N * N);
    MPI_Bcast(adj.data(), (int)adj.size(), MPI_BYTE, 0, MPI_COMM_WORLD);

    // ── Seed g_pmax with the best single-vertex solution ─────────────────
    for (int v = 0; v < N; v++) {
        if (cost[v] <= B && profit[v] > g_pmax) {
            g_pmax = profit[v];
            g_best = { v };
        }
    }

    // ── Dispatch ──────────────────────────────────────────────────────────
    if (P == 1) {
        runSequential(argv[2]);
    } else if (rank == 0) {
        runMaster(P, argv[2]);
    } else {
        runWorker(rank, P);
    }

    MPI_Finalize();
    return 0;
}