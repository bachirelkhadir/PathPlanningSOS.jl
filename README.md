# PathPlanningSOS.jl

# Installation

Open Julia and type-in the following commands to install the required packages.

```
using Pkg
Pkg.add("PathPlanningSOS")
Pkg.add("JuMP")
Pkg.add("MosekTools") # can be replaced another solver like CSDP or SCS
```


# Example usage

```
using PathPlanningSOS
using MosekTools
using JuMP

n = 2 # dimension of the space
max_deg_uv = 2 # degree of moment relaxation
num_pieces = 5 # number of linear pieces
num_iterations=20
weight_lenght= .1
random_seed = 3
a = [-1, -1]
b = [.1, 0.7]
edge_size = 1.

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)


moving_disk = [
    (t, x) -> sum( (x .- [1-2*t, 0]).^2 ) - .5^2
]

# compute optimal piece-wise linear trajectory
opt_trajectory = find_path_using_heuristic(n, moving_disk, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    weight_lenght,
    num_iterations,
    seed=random_seed)
```


Show the optimal trajectory.

```
# plot setup at time `t`
t = 0.
PathPlanningSOS.plot_at_time(0., edge_size, a, b,
                        moving_disk, opt_trajectory)
```
