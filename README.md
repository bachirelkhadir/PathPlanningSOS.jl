# PathPlanningSOS.jl


# Example usage

```
using MosekTools
using SumOfSquares

## Model
n = 2 # dimension of the space
max_deg_uv = 2
num_pieces = 5
scale_init = 1
reg = 0
num_iterations=100
weight_lenght= .1
seed = 3
a = [-1, -1]
b = [.1, 0.7]
edge_size = 1.

solver = with_optimizer(Mosek.Optimizer, QUIET=true)


moving_disk = [
    (t, x) -> dist_squared(x, [1-2*t, 0]) - .5^2  
]
                


opt_trajectory_pieces = find_path_using_heuristic(n, moving_disk, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    reg, weight_lenght,
    num_iterations,
    scale_init, seed=seed)
```
