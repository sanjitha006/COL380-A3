import random
import os

# To make the computation highly intensive, we use high budgets (1000)
# and adjust densities to push the Branch and Bound tree to its limits.
# Expected lines = N + Edges + 1
configs = [
    ("01_warmup",        400,  0.20,  800),  # ~16,000 edges  (16k lines)
    ("02_small_dense",   400,  0.40,  800),  # ~32,000 edges  (32k lines)
    ("03_med_sparse",    800,  0.10, 1000),  # ~32,000 edges  (33k lines)
    ("04_med_dense",     800,  0.20, 1000),  # ~64,000 edges  (65k lines)
    ("05_ta_clone",     1200, 0.103, 1000),  # ~74,000 edges  (75k lines) <- Matches your TA input
    ("06_large_1",      1200,  0.12, 1000),  # ~86,000 edges  (87k lines)
    ("07_large_2",      1200,  0.15, 1000),  # ~107,000 edges (108k lines)
    ("08_stress_1",     1200,  0.18, 1000),  # ~129,000 edges (130k lines)
    ("09_stress_2",     1200,  0.22, 1000)   # ~158,000 edges (159k lines)
]

for name, N, density, B in configs:
    filepath = f"testcases/{name}.txt"
    if os.path.exists(filepath):
        print(f"Skipping {name} (already exists)")
        continue
        
    print(f"Generating {name} (N={N}, Target Edges=~{int((N*(N-1)/2)*density)})...")
    
    # Pre-calculate to speed up python generation
    edges = []
    for i in range(N):
        for j in range(i+1, N):
            if random.random() < density:
                edges.append((i, j))
    
    with open(filepath, "w") as f:
        f.write(f"{N} {len(edges)} {B}\n")
        
        # Keep costs slightly higher so the budget prunes deep trees,
        # forcing the algorithm to do more complex backtracking.
        for i in range(N):
            f.write(f"{random.randint(10, 100)} {random.randint(10, 100)}\n")
            
        for u, v in edges:
            f.write(f"{u} {v}\n")
            
print("Done generating massive graphs.")
