/*
 * Sequential Branch and Bound for Budgeted Maximum Clique Problem
 * 
 * Algorithm:
 *   - Structural Bound: Greedy Graph Coloring (upper bound on clique size/profit)
 *   - Resource Bound:   Fractional Knapsack   (upper bound on achievable profit)
 *
 * Input Format:
 *   Line 1:  N E B          (vertices, edges, budget)
 *   Lines 2..N+1: p(v) c(v) for v = 0,1,...,N-1
 *   Lines N+2..N+E+1: u v   (undirected edge)
 *
 * Output Format:
 *   Line 1: Maximum profit found
 *   Line 2: Space-separated vertex IDs forming the optimal clique (ascending)
 */

#include <bits/stdc++.h>
using namespace std;

// ── Global graph data ──────────────────────────────────────────────────────
int N, E;
long long B;
vector<int>      profit, cost;
vector<vector<bool>> adj;   // adjacency matrix for O(1) neighbour checks

// ── Best solution found so far ─────────────────────────────────────────────
long long P_max = 0;
vector<int> best_clique;

// ══════════════════════════════════════════════════════════════════════════
//  BOUND 1 – Structural Bound via Greedy Graph Colouring
//  Partition C_cand into independent sets (colour classes).
//  Each class can contribute at most one vertex → take the best profit per class.
// ══════════════════════════════════════════════════════════════════════════
long long colorBound(const vector<int>& cand)
{
    // colour[i] = colour class index assigned to cand[i]
    int n = (int)cand.size();
    vector<int> colour(n, -1);
    vector<long long> class_best; // best profit per colour class

    for (int i = 0; i < n; ++i) {
        int v = cand[i];
        // Find the first colour class where v has no neighbour
        bool placed = false;
        for (int c = 0; c < (int)class_best.size(); ++c) {
            // Check if v is adjacent to any vertex already in class c
            bool conflict = false;
            for (int j = 0; j < i; ++j) {
                if (colour[j] == c && adj[v][cand[j]]) {
                    conflict = true;
                    break;
                }
            }
            if (!conflict) {
                colour[i] = c;
                class_best[c] = max(class_best[c], (long long)profit[v]);
                placed = true;
                break;
            }
        }
        if (!placed) {
            colour[i] = (int)class_best.size();
            class_best.push_back((long long)profit[v]);
        }
    }

    long long ub = 0;
    for (long long best : class_best) ub += best;
    return ub;
}

// ══════════════════════════════════════════════════════════════════════════
//  BOUND 2 – Resource Bound via Fractional Knapsack
//  Remaining budget = B - W_curr.
//  Sort cand by efficiency p(v)/c(v) descending, add greedily, allow fraction
//  for the last item. Returns upper bound on additional profit.
// ══════════════════════════════════════════════════════════════════════════
long long knapsackBound(const vector<int>& cand, long long rem_budget)
{
    if (rem_budget <= 0) return 0;

    // Sort by efficiency descending
    vector<int> order(cand.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b){
        // p[cand[a]]/c[cand[a]] > p[cand[b]]/c[cand[b]]
        // Use cross-multiplication to avoid floating point
        return (long long)profit[cand[a]] * cost[cand[b]] >
               (long long)profit[cand[b]] * cost[cand[a]];
    });

    long long ub = 0;
    long long budget_left = rem_budget;

    for (int idx : order) {
        int v = cand[idx];
        if (cost[v] <= budget_left) {
            ub += profit[v];
            budget_left -= cost[v];
        } else {
            // Fractional part
            // ub += profit[v] * budget_left / cost[v]  (floor is safe as upper bound needs ceiling)
            // Use ceiling to keep it an upper bound
            ub += ((long long)profit[v] * budget_left + cost[v] - 1) / cost[v];
            break;
        }
    }
    return ub;
}

// ══════════════════════════════════════════════════════════════════════════
//  Branch and Bound  (Algorithm 1 from the assignment)
//
//  cand    – candidate vertices (all pairwise adjacent to current clique)
//  P_curr  – total profit of current partial clique
//  W_curr  – total cost  of current partial clique
//  clique  – current partial clique (for solution reconstruction)
// ══════════════════════════════════════════════════════════════════════════
void findClique(vector<int> cand, long long P_curr, long long W_curr,
                vector<int>& clique)
{
    // ── Bound 1: Structural (Colouring) ──────────────────────────────────
    long long U_color = colorBound(cand);
    if (P_curr + U_color <= P_max) return;   // prune

    // ── Bound 2: Resource (Fractional Knapsack) ───────────────────────────
    long long U_knap = knapsackBound(cand, B - W_curr);
    if (P_curr + U_knap <= P_max) return;    // prune

    // ── Branching ─────────────────────────────────────────────────────────
    // Work on a copy so we can pop from the back
    while (!cand.empty()) {
        int v = cand.back();
        cand.pop_back();

        if (W_curr + cost[v] <= B) {
            // Update global best if this vertex alone (added to clique) is better
            long long new_profit = P_curr + profit[v];
            if (new_profit > P_max) {
                P_max = new_profit;
                best_clique = clique;
                best_clique.push_back(v);
            }

            // Build C_next = cand ∩ Neighbours(v)
            vector<int> C_next;
            C_next.reserve(cand.size());
            for (int u : cand) {
                if (adj[v][u]) C_next.push_back(u);
            }

            clique.push_back(v);
            findClique(C_next, P_curr + profit[v], W_curr + cost[v], clique);
            clique.pop_back();
        }
        // If budget exceeded, just skip (no recursion needed)
    }
}
int main(int argc, char *argv[])
{
    // Check if the correct number of arguments are provided
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    // Open the input file
    ifstream fin(argv[1]);
    if (!fin.is_open()) {
        cerr << "Error: Could not open input file " << argv[1] << "\n";
        return 1;
    }

    // Read from fin instead of cin
    fin >> N >> E >> B;

    profit.resize(N);
    cost.resize(N);
    for (int v = 0; v < N; ++v) {
        fin >> profit[v] >> cost[v];
    }

    adj.assign(N, vector<bool>(N, false));
    for (int i = 0; i < E; ++i) {
        int u, v;
        fin >> u >> v;
        adj[u][v] = adj[v][u] = true;
    }
    fin.close();

    // Initial candidate set: all vertices
    vector<int> cand(N);
    iota(cand.begin(), cand.end(), 0);
    sort(cand.begin(), cand.end(), [](int a, int b){
        return a < b; 
    });

    // Initialise P_max with the best single-vertex solution within budget
    for (int v = 0; v < N; ++v) {
        if (cost[v] <= B && profit[v] > P_max) {
            P_max = profit[v];
            best_clique = {v};
        }
    }

    vector<int> clique;
    findClique(cand, 0, 0, clique);

    // Open the output file
    ofstream fout(argv[2]);
    if (!fout.is_open()) {
        cerr << "Error: Could not open output file " << argv[2] << "\n";
        return 1;
    }

    // Write to fout instead of cout
    sort(best_clique.begin(), best_clique.end());
    fout << P_max << "\n";
    for (int i = 0; i < (int)best_clique.size(); ++i) {
        if (i) fout << " ";
        fout << best_clique[i];
    }
    fout << "\n";
    
    fout.close();

    return 0;
}