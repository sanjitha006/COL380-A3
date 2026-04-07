/*
 * MPI Parallel Branch-and-Bound — Budgeted Maximum Clique Problem
 *
 * Architecture (Pure Master-Worker)
 * ─────────────
 * • Master (Rank 0): Manages a dynamic task queue of subproblems.
 * • Workers (Ranks 1..P-1): Pull work, run the STRICT sequential BnB 
 * engine (local sorting for knapsack as per pseudocode), and report 
 * back tight bounds to prune the search space.
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
    TAG_WORK       = 1,   // master → worker  : here is a subproblem
    TAG_WORK_REQ   = 2,   // worker → master  : give me more work
    TAG_NO_MORE    = 3,   // master → worker  : no work left, terminate
    TAG_NEW_BEST   = 4,   // worker ↔ master  : updated P_max
    TAG_RESULT     = 5,   // worker → master  : final local best
};

// ═══════════════════════════════════════════════════════════════════════════
//  Graph — replicated on every rank
// ═══════════════════════════════════════════════════════════════════════════
static int           N, E, B;
static vector<int>   profit, cost;
static vector<uint8_t> adj;           // flat N×N adjacency matrix
inline bool edge(int u, int v) { return adj[u * N + v]; }

// ═══════════════════════════════════════════════════════════════════════════
//  Per-rank best solution
// ═══════════════════════════════════════════════════════════════════════════
static int         g_pmax = 0;
static vector<int> g_best;

// ═══════════════════════════════════════════════════════════════════════════
//  Subproblem + serialisation
// ═══════════════════════════════════════════════════════════════════════════
struct Sub {
    vector<int> cand;     // strictly ordered DESCENDING by profit
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

// ═══════════════════════════════════════════════════════════════════════════
//  Bounding functions (Strict Sequential Implementation)
// ═══════════════════════════════════════════════════════════════════════════

// Structural bound: greedy graph colouring.
static int colorBound(const vector<int>& cand) {
    int n = (int)cand.size();
    if (!n) return 0;
    
    vector<int> col(n, -1);
    int ub = 0;
    int num_colors = 0;
    
    for (int i = 0; i < n; i++) {
        int v = cand[i];
        bool placed = false;
        for (int c = 0; c < num_colors; c++) {
            bool conflict = false;
            for (int j = 0; j < i; j++) {
                if (col[j] == c && edge(v, cand[j])) { 
                    conflict = true; break; 
                }
            }
            if (!conflict) {
                col[i] = c;
                placed = true; break;
            }
        }
        if (!placed) {
            col[i] = num_colors++;
            ub += profit[v]; // First item in new class is the max profit
        }
    }
    return ub;
}

// Resource bound: fractional knapsack on remaining budget.
// Reverted to local sorting to strictly follow the pseudocode without pre-processing.
static int knapsackBound(const vector<int>& cand, int rem) {
    if (rem <= 0 || cand.empty()) return 0;
    
    vector<int> ord(cand.size()); 
    iota(ord.begin(), ord.end(), 0);
    
    // Sort locally by efficiency (p/c) using strict cross-multiplication
    sort(ord.begin(), ord.end(), [&](int a, int b) {
        int numA = profit[cand[a]], denA = cost[cand[a]];
        int numB = profit[cand[b]], denB = cost[cand[b]];
        if (numA * denB != numB * denA) return numA * denB > numB * denA;
        return profit[cand[a]] > profit[cand[b]];
    });
    
    int ub = 0, left = rem;
    for (int idx : ord) {
        int v = cand[idx];
        if (cost[v] <= left) { 
            ub += profit[v]; 
            left -= cost[v]; 
        } else { 
            ub += (profit[v] * left) / cost[v]; // floor of the fraction
            break; 
        }
    }
    return ub;
}

// ═══════════════════════════════════════════════════════════════════════════
//  expand() — process one BnB node
// ═══════════════════════════════════════════════════════════════════════════
static void expand(const Sub& sp, stack<Sub>& stk) {
    // Strict control flow as per typical BnB: Check bounds -> prune if needed.
    if (sp.P + colorBound(sp.cand) <= g_pmax) return;
    if (sp.P + knapsackBound(sp.cand, B - sp.W) <= g_pmax) return;

    int n = (int)sp.cand.size();

    // Loop through candidates. Going backwards pops from the "back" of the 
    // candidate list, exactly mimicking 'while (!cand.empty()) { pop_back() }'
    for (int i = n - 1; i >= 0; i--) {
        int v = sp.cand[i];
        if (sp.W + cost[v] > B) continue;

        int np = sp.P + profit[v];
        if (np > g_pmax) {
            g_pmax = np;
            g_best = sp.clique;
            g_best.push_back(v);
        }

        Sub child;
        child.P = np;
        child.W = sp.W + cost[v];
        child.clique = sp.clique;
        child.clique.push_back(v);
        child.cand.reserve(n - i - 1);

        for (int j = i + 1; j < n; j++) {
            if (edge(v, sp.cand[j])) {
                child.cand.push_back(sp.cand[j]);
            }
        }
        stk.push(move(child));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  rootDecompose() — rank 0 only
// ═══════════════════════════════════════════════════════════════════════════
static deque<Sub> rootDecompose(int target) {
    deque<Sub> q;
    Sub root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    
    // Initial assignment constraint: "sort the candidate vertices in descending order of profit"
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
            for (size_t i = 1; i < sp.cand.size(); i++) {
                if (edge(v, sp.cand[i])) with_v.cand.push_back(sp.cand[i]);
            }
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
    const int W = P - 1;

    deque<Sub> Q = rootDecompose(W * 8); 

    int global_best_p = g_pmax;
    vector<int> global_best_c = g_best;

    vector<bool> is_idle(P, false);
    vector<bool> terminated(P, false);
    terminated[0] = true;
    int idle_count = 0;
    int done_count = 0;

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

    const int RBUF_SZ = 1 << 18; 
    vector<int> rbuf(RBUF_SZ);

    while (done_count < W) {
        MPI_Status st;
        MPI_Recv(rbuf.data(), RBUF_SZ, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        int src = st.MPI_SOURCE;
        int tag = st.MPI_TAG;

        if (tag == TAG_WORK_REQ) {
            if (!Q.empty()) {
                auto b = packSub(Q.front()); Q.pop_front();
                MPI_Send(b.data(), (int)b.size(), MPI_INT, src, TAG_WORK, MPI_COMM_WORLD);
            } else {
                is_idle[src] = true;
                idle_count++;
                if (idle_count == W - done_count) {
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
            int np = rbuf[0];
            if (np > global_best_p) {
                global_best_p = np;
                for (int r = 1; r <= W; r++)
                    if (!terminated[r] && !is_idle[r] && r != src)
                        MPI_Send(&global_best_p, 1, MPI_INT, r, TAG_NEW_BEST, MPI_COMM_WORLD);
            }
        }
        else if (tag == TAG_RESULT) {
            int lp = rbuf[0];
            if (lp > global_best_p) {
                global_best_p = lp;
                global_best_c.clear();
                int csz = rbuf[1];
                for (int i = 0; i < csz; i++) global_best_c.push_back(rbuf[2 + i]);
            }
            terminated[src] = true;
            done_count++;
        }
    }

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
    const int RBUF_SZ = 1 << 18;
    vector<int> rbuf(RBUF_SZ);

    stack<Sub> local;
    bool wait_master = false;
    bool done        = false;

    auto sendResult = [&]() {
        vector<int> b;
        b.push_back(g_pmax);
        b.push_back((int)g_best.size());
        for (int v : g_best) b.push_back(v);
        MPI_Send(b.data(), (int)b.size(), MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
    };

    {
        MPI_Status st;
        MPI_Recv(rbuf.data(), RBUF_SZ, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        if (st.MPI_TAG == TAG_WORK) local.push(unpackSub(rbuf.data()));
        else done = true;
    }

    while (!done) {
        for (;;) {
            int flag = 0;
            MPI_Status st;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
            if (!flag) break;

            int src = st.MPI_SOURCE;
            int tag = st.MPI_TAG;
            int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
            MPI_Recv(rbuf.data(), cnt, MPI_INT, src, tag, MPI_COMM_WORLD, &st);

            if (tag == TAG_NEW_BEST) {
                int np = rbuf[0];
                if (np > g_pmax) g_pmax = np;
            }
            else if (tag == TAG_WORK && wait_master) {
                local.push(unpackSub(rbuf.data()));
                wait_master = false;
            }
            else if (tag == TAG_NO_MORE && wait_master) {
                done = true;
                wait_master = false;
            }
        }

        if (done) break;
        if (wait_master) continue;

        if (local.empty()) {
            int d = 0;
            MPI_Send(&d, 1, MPI_INT, 0, TAG_WORK_REQ, MPI_COMM_WORLD);
            wait_master = true;
            continue;
        }

        Sub sp = local.top(); local.pop();
        int prev = g_pmax;
        expand(sp, local);

        if (g_pmax > prev) {
            MPI_Send(&g_pmax, 1, MPI_INT, 0, TAG_NEW_BEST, MPI_COMM_WORLD);
        }
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
    
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });
    
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
        
        adj.assign((size_t)N * N, 0);
        for (int i = 0; i < E; i++) {
            int u, v; fin >> u >> v;
            adj[u * N + v] = adj[v * N + u] = 1;
        }
        fin.close();
    }

    {
        int neb[3] = { N, E, B };
        MPI_Bcast(neb, 3, MPI_INT, 0, MPI_COMM_WORLD);
        N = neb[0]; E = neb[1]; B = neb[2];
    }
    profit.resize(N); cost.resize(N);
    MPI_Bcast(profit.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cost.data(),   N, MPI_INT, 0, MPI_COMM_WORLD);
    adj.resize((size_t)N * N);
    MPI_Bcast(adj.data(), (int)adj.size(), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Initial Greedy pass for a strong starting P_max
    for (int v = 0; v < N; v++) {
        if (cost[v] <= B && profit[v] > g_pmax) {
            g_pmax = profit[v];
            g_best = { v };
        }
    }

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