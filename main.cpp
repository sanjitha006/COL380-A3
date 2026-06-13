#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;

enum Msg {
    NEED_WORK  = 1,  
    TAKE_WORK= 2,  
    NO_WORK= 3,  
    BEST_P= 4,  
    TOKEN= 5,   // Dijkstra token
    FINISH= 6   
};

enum Color { WHITE = 0, BLACK = 1 };

static int N;
static int E;
static int B;
static vector<int>   profit, cost;


static int NWORDS = 0;
static vector<uint64_t> adj_bits;

inline const uint64_t* adjRow(int u) {
    return adj_bits.data() + (size_t)u * NWORDS;
}
inline bool fastEdge(int u, int v) {
    return (adjRow(u)[v >> 6] >> (v & 63)) & 1ULL;
}


static int max_p = 0;
static vector<int> best;


static vector<int> sorted;
static vector<uint8_t> in_cand;

static long long msgs_sent = 0;
static long long msgs_recv = 0;
static vector<uint64_t> g_color_buf;

struct Prob {
    vector<int> cand;
    vector<int> clique;
    int P = 0;
    int W = 0;
};

static vector<int> packProb(const Prob& s) {
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

static Prob unpackProb(const int* d) {
    Prob s; int i = 0;
    int cs = d[i++]; s.cand.assign(d + i, d + i + cs);   
    i += cs;
    int cl = d[i++]; s.clique.assign(d + i, d + i + cl);
    i += cl;
    s.P = d[i++];
    s.W = d[i++];
    return s;
}

static vector<int> packBatch(const vector<Prob>& subs) {
    vector<int> b;
    b.push_back((int)subs.size());
    for (const Prob& s : subs) {
        auto p = packProb(s);
        b.push_back((int)p.size());
        b.insert(b.end(), p.begin(), p.end());
    }
    return b;
}

static vector<Prob> unpackBatch(const int* d) {
    vector<Prob> out;
    int idx = 0, n = d[idx++];
    for (int i = 0; i < n; i++) {
        int sz = d[idx++];
        out.push_back(unpackProb(d + idx));
        idx += sz;
    }
    return out;
}


static int structbound(const vector<int>& cand) {
   
    if (cand.empty()) return 0;

    int ncolor = 0;
    int res = 0;

    for (int i = 0; i < cand.size(); i++) {
        int v=cand[i];
        const uint64_t* vrow = adjRow(v);   
        int fnd=0;
        
        for (int c = 0; c < ncolor; c++) {
            
            const uint64_t* cs = g_color_buf.data() + c * NWORDS;
            int p=0;
            for (int w = 0; w < NWORDS; w++) {
                if (vrow[w] & cs[w]) { p=1; break; }
            }
            if (!p) {
           
                g_color_buf[c * NWORDS + (v >> 6)] |= 1ULL << (v & 63);
                fnd = true;
                break;
            }
        }

        if (!fnd) {
            size_t new_size =(ncolor + 1) * NWORDS;
            if (g_color_buf.size() < new_size)
                g_color_buf.resize(new_size, 0ULL);
            fill(g_color_buf.begin() + ncolor * NWORDS,
                 g_color_buf.begin() + new_size, 0ULL);
            g_color_buf[ncolor * NWORDS + (v >> 6)] |= 1ULL << (v & 63);
            ncolor++;
            res += profit[v];
        }
    }
    return res;
}


static int knapsackBound(const vector<int>& cand, int used) {
    if (used>=B || cand.empty()) return 0;

    for (int v : cand) in_cand[v] = 1;
    int n=cand.size();
    int res = 0;
    int cnt=0;
    

    for (int v : sorted) {
        if (!in_cand[v]) continue;

        if ((cost[v]+used) <= B) {
            res += profit[v];
            used += cost[v];
        } else {
            res += (int)(((long long)profit[v] * (B - used)) / cost[v]);
            break;
        }
        cnt;
        if(cnt==n) break;

    }

    for (int v : cand) in_cand[v] = 0;

    return res;
}


static void expand(const Prob& sp, deque<Prob>& stk) {
    if (sp.P + structbound(sp.cand) <= max_p) return;
    if (sp.P + knapsackBound(sp.cand,sp.W) <= max_p) return;

    int n =sp.cand.size();
    for (int i = n - 1; i >= 0; i--) {
        int v = sp.cand[i];
        if (sp.W + cost[v] > B) continue;

        int np = sp.P + profit[v];
        if (np > max_p) {
            max_p = np;
            best = sp.clique;
            best.push_back(v);
        }
        if (sp.W + cost[v] == B) continue;

        Prob newprob;
        newprob.P = np;
        newprob.W = sp.W + cost[v];
        newprob.clique = sp.clique;
        newprob.clique.push_back(v);
        newprob.cand.reserve(n - i - 1);

     
        const uint64_t* vrow = adjRow(v);
        for (int j = i + 1; j < n; j++) {
            int u = sp.cand[j];
            if ((vrow[u >> 6] >> (u & 63)) & 1ULL)
                newprob.cand.push_back(u);
        }

        stk.push_back(move(newprob));
    }
}


static deque<Prob> initroot(int target) {
    deque<Prob>q;
    Prob root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a > b;
    });
    q.push_back(root);

    while (q.size()<target && !q.empty()) {
        // auto it = max_element(q.begin(), q.end(), [](const Prob& a, const Prob& b){
        //     return a.cand.size() < b.cand.size();
        // });
        
       // auto it=q.front();
        Prob sp = q.front(); q.pop_front();
        if (sp.cand.empty()) { q.push_back(sp); break; }

        if (sp.P + structbound(sp.cand) <= max_p) continue;
        if (sp.P + knapsackBound(sp.cand, sp.W) <= max_p) continue;

        int v = sp.cand[0];

        Prob without = sp;
        without.cand.erase(without.cand.begin());
        if (!without.cand.empty()) q.push_back(without);

        if (sp.W + cost[v] < B) {
            Prob with_v;
            with_v.P = sp.P + profit[v];
            with_v.W = sp.W + cost[v];
            with_v.clique = sp.clique;
            with_v.clique.push_back(v);
            const uint64_t* vrow = adjRow(v);
            for (size_t i = 1; i < sp.cand.size(); i++) {
                int u = sp.cand[i];
                if ((vrow[u >> 6] >> (u & 63)) & 1ULL)
                    with_v.cand.push_back(u);
            }
            if (with_v.P > max_p) { max_p = with_v.P; best = with_v.clique; }
            q.push_back(with_v);
        }
        else if((sp.W+cost[v])==B){
            if ((sp.P+profit[v]) > max_p) { max_p = sp.P+profit[v];sp.clique.push_back(v); best =sp.clique; }
        }
    }
    return q;
}

static void runDecentralized(int rank, int P, const char* out_filename) {
    vector<int> rbuf;
    deque<Prob> local;

    vector<vector<int>> init_sbufs(P);
    vector<MPI_Request> init_reqs(P, MPI_REQUEST_NULL);

    bool is_idle = false;
    bool done = false;
    bool steal_pending = false;
    bool steal_exhausted = false; // True if our steal token made a full circle

    //dijk token vars
    int  my_color = WHITE;
    bool has_token = (rank == 0);
    int  token_color = WHITE;

    int poll_counter = 0;
    const int POLL_FREQ = 1024;

    srand(time(NULL) + rank);

    deque<Prob> initial_chunks;
    if (rank == 0) {
        const int CHUNK_MULTIPLIER = 4;
        initial_chunks = initroot(P * CHUNK_MULTIPLIER);
    }

    MPI_Bcast(&max_p, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        vector<vector<Prob>> peer_batches(P);
        int curr_peer = 0;
        
        while (!initial_chunks.empty()) {
            peer_batches[curr_peer].push_back(initial_chunks.front());
            initial_chunks.pop_front();
            curr_peer = (curr_peer + 1) % P;
        }

        for (auto& s : peer_batches[0]) {
            local.push_back(move(s));
        }

        for (int r = 1; r < P; r++) {
            if (!peer_batches[r].empty()) {
                init_sbufs[r] = packBatch(peer_batches[r]);
                MPI_Isend(init_sbufs[r].data(), (int)init_sbufs[r].size(),
                          MPI_INT, r, TAKE_WORK, MPI_COMM_WORLD, &init_reqs[r]);
                msgs_sent++; 
            } else {
                init_sbufs[r] = {0};
                MPI_Isend(init_sbufs[r].data(), 1, MPI_INT,
                          r, NO_WORK, MPI_COMM_WORLD, &init_reqs[r]);
                msgs_sent++; 
            }
        }

        if (local.empty()) is_idle = true;
    }
    else {
        MPI_Status st;
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
        if ((size_t)cnt > rbuf.size()) rbuf.resize((size_t)cnt);

        MPI_Recv(rbuf.data(), cnt, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        msgs_recv++; 

        if (st.MPI_TAG == TAKE_WORK) {
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
            MPI_Send(&out_token, 1, MPI_INT, next_rank, TOKEN, MPI_COMM_WORLD);
            msgs_sent++; // Track token send
            my_color = WHITE;
            has_token = false;
        }
    };

    while (!done) {
     
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
                msgs_recv++; 

                switch (tag) {
                    case BEST_P: {
                        int np = rbuf[0];
                        if (np > max_p) max_p = np;
                        break;
                    }

                    case NEED_WORK: {
                        int origin_rank = rbuf[0];

                        if (local.size() >= 2) {
                            int half = (int)local.size() / 2;
                            vector<Prob> donated;
                            for (int i = 0; i < half; i++) {
                                donated.push_back(local.front());
                                local.pop_front();
                            }
                            auto sbuf = packBatch(donated);
                            MPI_Send(sbuf.data(), (int)sbuf.size(), MPI_INT, origin_rank, TAKE_WORK, MPI_COMM_WORLD);
                            msgs_sent++;
                            my_color = BLACK;
                        } else {
                         
                            if (origin_rank == rank) {
                                steal_pending = false;
                                steal_exhausted = true; // The whole ring is empty.
                            } else {
                                int next_node = (rank + 1) % P;
                                MPI_Send(&origin_rank, 1, MPI_INT, next_node, NEED_WORK, MPI_COMM_WORLD);
                                msgs_sent++; 
                            }
                        }
                        break;
                    }

                    case TAKE_WORK: {
                        auto subs = unpackBatch(rbuf.data());
                        for (auto& s : subs) local.push_back(move(s));
                        
                        steal_pending = false;
                        steal_exhausted = false; 
                        is_idle = false;
                        break;
                    }

                    case NO_WORK: break; 

                    case TOKEN: {
                        if (rank == 0) {
                            int in_color = rbuf[0];
                            if (in_color == WHITE && is_idle) {
                                int d = 0;
                                for (int r = 1; r < P; r++) {
                                    MPI_Send(&d, 1, MPI_INT, r, FINISH, MPI_COMM_WORLD);
                                    msgs_sent++; // Track FINISH send
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

                    case FINISH: {
                        done = true;
                        break;
                    }
                }
            }
        }

        if (done) break;

        if (local.empty()) {
            is_idle = true;
        }
        passToken();

        if (is_idle) {
         
            if (!steal_pending && !steal_exhausted && P > 1) {
                int origin_rank = rank;
                int next_node = (rank + 1) % P;
                MPI_Send(&origin_rank, 1, MPI_INT, next_node, NEED_WORK, MPI_COMM_WORLD);
                msgs_sent++; 
                steal_pending = true;
            }
            continue; 
        }

        Prob sp = local.back();
        local.pop_back();

        int prev = max_p;
        expand(sp, local);

        if (max_p > prev) {
            for (int r = 0; r < P; r++) {
                if (r != rank) {
                    MPI_Request req;
                    MPI_Isend(&max_p, 1, MPI_INT, r, BEST_P, MPI_COMM_WORLD, &req);
                    MPI_Request_free(&req); 
                    msgs_sent++; 
                }
            }
        }
    }

    while (true) {
        long long local_counts[2] = {msgs_sent, msgs_recv};
        long long global_counts[2] = {0, 0};
        MPI_Allreduce(local_counts, global_counts, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (global_counts[0] == global_counts[1]) {
            break; 
        }

        int flag = 0;
        MPI_Status st;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
        if (flag) {
            int cnt; MPI_Get_count(&st, MPI_INT, &cnt);
            if ((size_t)cnt > rbuf.size()) rbuf.resize((size_t)cnt);
            MPI_Recv(rbuf.data(), cnt, MPI_INT, st.MPI_SOURCE, st.MPI_TAG, MPI_COMM_WORLD, &st);
            msgs_recv++; 
            
            if (st.MPI_TAG == BEST_P && rbuf[0] > max_p) {
                max_p = rbuf[0];
            }
        }
    }

    int actual_local_profit = 0;
    for (int v : best) actual_local_profit += profit[v];

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
            best.clear();
            for (int i = 0; i < csz; i++) best.push_back(rbuf[2 + i]);
        }

        sort(best.begin(), best.end());
        ofstream fout(out_filename);
        
        fout << global_max_profit << "\n";
        for (int i = 0; i < (int)best.size(); i++) {
            if (i) fout << " ";
            fout << best[i];
        }
        fout << "\n";
        fout.close();
    }
    else {
        if (rank == winner_rank) {
            vector<int> b;
            b.push_back(actual_local_profit);
            b.push_back((int)best.size());
            for (int v : best) b.push_back(v);
            MPI_Send(b.data(), (int)b.size(), MPI_INT, 0, 999, MPI_COMM_WORLD);
        }
    }
}


static void runSequential(const char* out_filename) {
    deque<Prob> stk;
    Prob root;
    root.cand.resize(N);
    iota(root.cand.begin(), root.cand.end(), 0);
    sort(root.cand.begin(), root.cand.end(), [](int a, int b){
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });
    stk.push_back(root);

    while (!stk.empty()) {
        Prob sp = stk.back(); stk.pop_back();
        expand(sp, stk);
    }

    sort(best.begin(), best.end());
    ofstream fout(out_filename);
   
    fout << max_p << "\n";
    for (int i = 0; i < (int)best.size(); i++) {
        if (i) fout << " ";
        fout << best[i];
    }
    fout << "\n";
    fout.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank,P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    
    if (rank == 0) {
        ifstream fin(argv[1]);

        fin >> N >> E >> B;
        profit.resize(N); cost.resize(N);
        for (int i = 0; i < N; i++) fin >> profit[i] >> cost[i];

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

    NWORDS = (N + 63) / 64;
    profit.resize(N); cost.resize(N);
    MPI_Bcast(profit.data(), N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(cost.data(),   N, MPI_INT, 0, MPI_COMM_WORLD);
    adj_bits.resize((size_t)N * NWORDS, 0ULL);
    MPI_Bcast(adj_bits.data(), (int)(adj_bits.size() * sizeof(uint64_t)),
              MPI_BYTE, 0, MPI_COMM_WORLD);

    in_cand.assign(N, 0);
    sorted.resize(N);
    iota(sorted.begin(), sorted.end(), 0);
    sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        long long numA = profit[a], denA = cost[a];
        long long numB = profit[b], denB = cost[b];
        if (numA * denB != numB * denA) return numA * denB > numB * denA;
        if (profit[a] != profit[b]) return profit[a] > profit[b];
        return a < b;
    });

    g_color_buf.reserve((size_t)N * NWORDS);

    for (int v = 0; v < N; v++) {
        if (cost[v] <= B && profit[v] > max_p) {
            max_p = profit[v];
            best = { v };
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