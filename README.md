\# General info

This is a study pet-project, connected with 3D solution of Navier-Stokes equations. The domain is a simple cube.



\## Python files



\+ \*\*visualize.py\*\* - module with reading data and its single-frame visualization

\+ \*\*video\_maker\_velocity.py\*\* - uses visualize.py, create video of grid-function velocity change

\+ \*\*video\_maker\_velocity\_pair.py\*\* - uses visualize.py, create video of grid-function velocity change in comparsion with original one change



\## C++ files

\+ src

&nbsp; - inout (input-output functionality)

&nbsp; - main

&nbsp; - mesh\_n\_model (init mesh and model parameters, initial conditions, precise solution)

&nbsp; - compute\_flow (construct scheme, solve equations)

\+ headers

&nbsp; - settings (contains model parameters structure and includes some basic libraries)

\+ config (contains model data aka domain X-Y-Z size, domain partition, time partition, Reynolds number etc)



\## Other

\+ .vscode (vscode settings)

\+ makefile (it probably needs to be modified)



> Absolute path were set inside program files, change them in order to be able run program\\

> makefile is not need if u want to run program via VS Code on Windows

> Eigen is needed to run C++ program

