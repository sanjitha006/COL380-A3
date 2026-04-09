import random
import os

random.seed(42)

configs = [
    ("01_tight_s",   300,  0.55,    1,  20,   5),
    ("02_tight_m",   400,  0.55,    1,  20,   5),
    ("03_dense_s",   300,  0.65,    1,  20,   6),
    ("04_dense_m",   400,  0.65,    1,  20,   6),
    ("05_ta_clone",  500,  0.60,    1,  20,   9),
    ("06_hard_500",  500,  0.65,    1,  20,   7),
    ("07_hard_600",  600,  0.60,    1,  20,   8),
    ("08_stress_a",  600,  0.65,    1,  20,   7),
    ("09_stress_b",  700,  0.60,    1,  20,   8),
    ("10_ultra",     500,  0.70,    1,  15,   7),
]

os.makedirs("testcases", exist_ok=True) # Works in Python 3.2+

for name, N, density, clo, chi, bmul in configs:
    filepath = "testcases/{}.txt".format(name)
    if os.path.exists(filepath):
        print("Skipping {} (already exists)".format(name))
        continue

    avg_cost = (clo + chi) / 2.0
    B = int(bmul * avg_cost)

    print("Generating {}: N={}, density={}, cost=[{},{}], B={} (~{} nodes fit)...".format(
        name, N, density, clo, chi, B, bmul))

    edges = []
    for i in range(N):
        for j in range(i + 1, N):
            if random.random() < density:
                edges.append((i, j))

    with open(filepath, "w") as f:
        f.write("{} {} {}\n".format(N, len(edges), B))
        for i in range(N):
            f.write("{} {}\n".format(random.randint(1, 100), random.randint(clo, chi)))
        for u, v in edges:
            f.write("{} {}\n".format(u, v))

    print("  -> {} edges written.".format(len(edges)))

print("Done generating hard testcases.")
