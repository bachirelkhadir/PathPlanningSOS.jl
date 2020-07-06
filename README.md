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

## Model
n = 2 # dimension of the space
max_deg_uv = 2
num_pieces = 5
scale_init = 1
reg = 0
num_iterations=10
weight_lenght= .1
seed = 3
a = [-1, -1]
b = [.1, 0.7]
edge_size = 1.

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

# functions describing the dynamic obstacles
moving_disk = [
    (t, x) -> dist_squared(x, [1-2*t, 0]) - .5^2
]

# compute optimal path
opt_trajectory_pieces = find_path_using_heuristic(n, moving_disk, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    reg, weight_lenght,
    num_iterations,
    scale_init, seed=seed)


# Create an animation of the particule moving through the dynamic obstacles
fn = "animation.mp4"
PathPlanningSOS.make_animation(fn, opt_trajectory_pieces, moving_disk, edge_size, a, b, )

# If you are using jupyter notebook, the following will display the video in the browser
PathPlanningSOS.show_animation_in_notebook(fn)

```
