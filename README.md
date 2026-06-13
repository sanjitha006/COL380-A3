# Max-Weight-Clique-MPI

**Language & Method:** This project is implemented in C++ and parallelized exclusively using MPI, utilizing a recursive Branch and Bound algorithm with mandatory Structural (Graph Coloring) and Resource (Fractional Knapsack) bounding heuristics.

**Motive & Objective:** Given an undirected graph where vertices have specific profits and costs, the goal is to find a complete subgraph (clique) that maximizes total profit without exceeding a strict maximum global budget.

**Execution Constraints:** The program must compile via standard MPI commands (`mpic++ -O3 -std=c++17`) and execute efficiently, scaling up to 16 processes to yield the exact optimal clique.
