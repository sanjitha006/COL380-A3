import sys

def verify_clique(input_file, output_file):
    # ─── 1. PARSE INPUT FILE ─────────────────────────────────────────
    try:
        with open(input_file, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        print(f"[ERROR] Could not find input file: {input_file}")
        return

    # Read N, E, B
    try:
        N, E, B = map(int, lines[0].split())
    except ValueError:
        print("[ERROR] Line 1 of input must be: N E B")
        return

    # Read Profits and Costs
    profits = []
    costs = []
    for i in range(1, N + 1):
        p, c = map(int, lines[i].split())
        profits.append(p)
        costs.append(c)

    # Read Edges into an Adjacency Matrix (Set of Tuples for fast lookup)
    edges = set()
    for i in range(N + 1, N + 1 + E):
        u, v = map(int, lines[i].split())
        edges.add((u, v))
        edges.add((v, u)) # Undirected graph

    # ─── 2. PARSE OUTPUT FILE ────────────────────────────────────────
    try:
        with open(output_file, 'r') as f:
            out_lines = [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        print(f"[ERROR] Could not find output file: {output_file}")
        return

    if len(out_lines) < 1:
        print("[ERROR] Output file is empty.")
        return

    reported_profit = int(out_lines[0])
    
    # Handle the case where the clique might be empty (profit 0)
    clique_vertices = []
    if len(out_lines) > 1 and out_lines[1]:
        clique_vertices = list(map(int, out_lines[1].split()))

    # ─── 3. RUN VALIDATION CHECKS ────────────────────────────────────
    print(f"--- Verifying {output_file} against {input_file} ---")
    
    passed_all = True

    # Check A: Are all vertices valid IDs?
    invalid_vertices = [v for v in clique_vertices if v < 0 or v >= N]
    if invalid_vertices:
        print(f"[FAIL] Invalid vertex IDs found in solution: {invalid_vertices}")
        passed_all = False
    else:
        print("[PASS] All vertex IDs are within valid range (0 to N-1).")

    # Check B: Is it sorted in ascending order?
    is_sorted = all(clique_vertices[i] < clique_vertices[i+1] for i in range(len(clique_vertices)-1))
    if not is_sorted:
        print("[FAIL] The vertices in the output are NOT sorted in ascending order.")
        passed_all = False
    else:
        print("[PASS] Vertices are sorted in ascending order.")

    # Check C: Is it a valid clique? (Every pair must share an edge)
    is_clique = True
    for i in range(len(clique_vertices)):
        for j in range(i + 1, len(clique_vertices)):
            u = clique_vertices[i]
            v = clique_vertices[j]
            if (u, v) not in edges:
                print(f"[FAIL] Not a valid clique! Missing edge between {u} and {v}.")
                is_clique = False
                break
        if not is_clique:
            break
    
    if is_clique:
        print("[PASS] The solution forms a valid complete subgraph (Clique).")
    else:
        passed_all = False

    # Check D: Does it violate the budget constraint?
    actual_cost = sum(costs[v] for v in clique_vertices)
    if actual_cost > B:
        print(f"[FAIL] Budget exceeded! Max Budget: {B}, Actual Cost: {actual_cost}")
        passed_all = False
    else:
        print(f"[PASS] Budget is respected. (Used {actual_cost} / {B})")

    # Check E: Does the calculated profit match the reported profit?
    actual_profit = sum(profits[v] for v in clique_vertices)
    if actual_profit != reported_profit:
        print(f"[FAIL] Profit mismatch! Printed profit: {reported_profit}, Actual sum of vertices: {actual_profit}")
        passed_all = False
    else:
        print(f"[PASS] Profit matches the reported maximum ({actual_profit}).")

    # ─── 4. FINAL VERDICT ────────────────────────────────────────────
    print("-" * 50)
    if passed_all:
        print("✅ SUCCESS: The output is a valid budgeted clique!")
    else:
        print("❌ FAILED: The output violates one or more constraints.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python verify.py <input.txt> <output.txt>")
        sys.exit(1)
        
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    verify_clique(in_file, out_file)