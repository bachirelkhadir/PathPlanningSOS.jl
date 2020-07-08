# PathPlanningSOS.jl


[![Path Planning using SOS - Example 1](https://raw.githubusercontent.com/bachirelkhadir/PathPlanningSOS.jl/master/doc/path_planning_animation.gif)](https://youtu.be/8VXckZWe-VQ)

[![Path Planning using SOS Example 2](https://raw.githubusercontent.com/bachirelkhadir/PathPlanningSOS.jl/master/doc/path_planning_animation_2.gif)](https://youtu.be/6ThKwE0B9yA)


The role of this package is to find a trajectory (i.e., a
continuously-differentiable
function) 
<img src="https://render.githubusercontent.com/render/math?math=x: [0, T] \rightarrow  \mathbb R^n" />
that
- starts on a point `a` and ends on a point `b`, i.e. 
- satisfies `m` obstacle constraints given by inequalities
   <img src="https://render.githubusercontent.com/render/math?math=g_i(t, x(t)) \ge 0 " />
   <img src="https://render.githubusercontent.com/render/math?math=\forall t \in [0, T]" />
   <img src="https://render.githubusercontent.com/render/math?math=\forall i \in \{1, \ldots, s\}" />


# Installation

Open Julia and type-in the following commands to install the required packages.

```
using Pkg
# download this package
Pkg.add(PackageSpec(url="https://github.com/bachirelkhadir/PathPlanningSOS.jl"))
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

# Obstacle contraints
# here there is a single obstacle: a disk of center [1-2*t, 0] and radius .5
moving_disk = [
    (t, x) -> sum( (x .- [1-2*t, 0]).^2 ) - .5^2
]

# compute optimal piece-wise linear trajectory
opt_trajectory = find_path_using_heuristic(n, moving_disk, edge_size, a, b,
    max_deg_uv, num_pieces, solver,
    weight_lenght,
    num_iterations,
    seed=random_seed)
# opt_trajectory is a function that takes t as input, 
# and gives the location of the particle at time t.
```
## Optional

Plot the obtained trajectory.

```
using PyPlot

# plot setup at time `t`
for (i,t)=enumerate(0:.1:1)
    PyPlot.figure()
    PathPlanningSOS.plot_at_time(t, edge_size, a, b,
                        moving_disk, opt_trajectory)
    PyPlot.title("t = $t")
    fn_index = lpad(i, 3, "0")
    PyPlot.savefig("path_planning_frame_$(fn_index).png")
end
```


Make a video animation (make sure ffmpeg is installed)
```
;ffmpeg -r 10 -i path_planning_frame_%03d.png -vcodec libx264 -pix_fmt yuv420p path_planning_animation.mp4 -y
```

Show video on jupyter notebook
```
using Base64
base64_video = base64encode(open("path_planning_animation.mp4"))
display("text/html", """<video controls src="data:video/x-m4v;base64,$base64_video">""")
```




