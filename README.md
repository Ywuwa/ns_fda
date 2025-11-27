## General info

This is a study pet-project, connected with 3D solution of Navier-Stokes equations. The domain is a simple cube.

## Python files

+ **visualize.py** - module with reading data and its single-frame visualization
+ **video_maker_velocity.py** - uses visualize.py, create video of grid-function velocity change
+ **video_maker_velocity_pair.py** - uses visualize.py, create video of grid-function velocity change in comparsion with original one change



## C++ files

+ src
  - inout (input-output functionality)
  - main
  - mesh_n_model (init mesh and model parameters, initial conditions, precise solution)
  - compute_flow (construct scheme, solve equations)
+ headers
  - settings (contains model parameters structure and includes some basic libraries)
+ config (contains model data aka domain X-Y-Z size, domain partition, time partition, Reynolds number etc)

## Other
+ .vscode (vscode settings)
+ makefile (it probably needs to be modified)

> Absolute path were set inside program files, change them in order to be able run program
> makefile is not need if u want to run program via VS Code on Windows
> Eigen is needed to run C++ program

