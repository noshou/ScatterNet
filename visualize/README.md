# visualize

Optional visual and statistical harnesses for the forward mod

These depend on GLMakie, which lives in test/Project.toml`, so run them with `--project=test`.


| file             | what it does                                                                                |
| ------------------ | --------------------------------------------------------------------------------------------- |
| `plastic_vis.jl` | scatter-plots the plastic sequence on the unit sphere                                       |
| `sasa_bench.jl`  | accuracy / `area_tol` skip-firing / speed sweeps over ~450 systems, plus summary figures |
| `sasa_vis.jl`    | four occlusion regimes + mesh-convergence panels, points coloured by exposed/occluded state |

```
julia --project=visualize -e 'include("visualize/sasa_vis.jl");   vis_sasa_cases()'
julia --project=visualize -e 'include("visualize/plastic_vis.jl"); vis_plastic_points(2000)'
```
