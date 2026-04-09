import sys

def check(input_file, output_file, expected_profit):
    with open(input_file, 'r') as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]
    N, E, B = map(int, lines[0].split())
    profits, costs = [], []
    for i in range(1, N + 1):
        p, c = map(int, lines[i].split())
        profits.append(p), costs.append(c)
    
    edges = set()
    for i in range(N + 1, N + 1 + E):
        u, v = map(int, lines[i].split())
        edges.add((u, v)), edges.add((v, u))

    with open(output_file, 'r') as f:
        out_lines = [l.strip() for l in f.readlines() if l.strip()]
    if not out_lines:
        print("FAIL: Output file is empty")
        sys.exit(1)
        
    reported_profit = int(out_lines[0])
    clique = list(map(int, out_lines[1].split())) if len(out_lines) > 1 and out_lines[1] else []

    if any(v < 0 or v >= N for v in clique):
        print("FAIL: Invalid vertex IDs")
        sys.exit(1)
    if not all(clique[i] < clique[i+1] for i in range(len(clique)-1)):
        print("FAIL: Vertices not sorted in ascending order")
        sys.exit(1)
    for i in range(len(clique)):
        for j in range(i + 1, len(clique)):
            if (clique[i], clique[j]) not in edges:
                print("FAIL: Not a clique! Missing edge ({}, {})".format(clique[i], clique[j]))
                sys.exit(1)
                
    actual_cost = sum(costs[v] for v in clique)
    if actual_cost > B:
        print("FAIL: Budget exceeded ({} > {})".format(actual_cost, B))
        sys.exit(1)
        
    actual_profit = sum(profits[v] for v in clique)
    if actual_profit != reported_profit:
        print("FAIL: Profit math mismatch (Printed: {}, Actual: {})".format(reported_profit, actual_profit))
        sys.exit(1)
        
    if actual_profit != expected_profit:
        print("FAIL: Suboptimal. Expected Pmax={}, Got={}".format(expected_profit, actual_profit))
        sys.exit(1)

    print("PASS")
    sys.exit(0)

check(sys.argv[1], sys.argv[2], int(sys.argv[3]))
