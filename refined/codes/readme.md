
### <span style="color: darkred; text-decoration: underline;">FOSA Code Overview</span>

This document provides an overview of the **FOSA** code, written in **Fortran**. The code utilizes a **6th-order finite differencing** method for spatial derivatives and a **3rd-order implicit Runge-Kutta** method for time integration ([reference](https://doi.org/10.48550/arXiv.astro-ph/0109497)). Designed with a modular structure, each module performs a distinct task, facilitating both flexibility and clarity in code management.

#### <span style="color: darkred; text-decoration: underline;">Main Modules</span>

1. **`alpha_profile.f90`**  
   - Defines the alpha profile. This module is relevant only if there is a kinematic alpha effect (i.e., if $ R_\alpha \neq 0 $.)
   
2. **`equations.f90`**  
   - Contains all routines for the evolution equations. Detailed comments are included within each subroutine.
   
3. **`field_initial.f90`**  
   - Sets the initial conditions for the simulation.
   
4. **`grid.f90`**  
   - Constructs the grid for the simulation. Modify spatial resolution and time steps in this module.
   
5. **`parameters_constants.f90`**  
   - Sets all control parameters. Each trial name and corresponding parameter values are listed here.
   
6. **`run_file.f90`**  
   - Manages output generation and the execution of all routines.
   
7. **`seed.f90`**  
   - Generates a random Gaussian field as part of the initial conditions.
   
8. **`spatial_derivatives.f90`**  
   - Calculates spatial derivatives using 6th-order finite differencing. This module includes options for various ghost zone types and supports finite differencing of orders from 2 to 12. Avoid modifications unless necessary.
   
9. **`time_stepping.f90`**  
   - Implements the 3rd-order implicit Runge-Kutta time integration scheme, with an additional RK4 option. Avoid modifications unless necessary.
   
10. **`velocity_profile.f90`**  
    - Defines the velocity profile for the simulation.

#### <span style="color: darkred; text-decoration: underline;">Plotting and Running the Code</span>

- **`plots.py`** and **`compare_plots.py`**  
  - `plots.py` generates plots from simulation output.  
  - `compare_plots.py` compares outputs from multiple simulations.
  
- **`run.sh`**  
  - This script compiles and executes the code. Run with the command: `bash run.sh`.

#### <span style="color: darkred; text-decoration: underline;">Code Capabilities</span>

Currently, this code can be used to replicate results from the following works and notes:

1. **L.Chamandy et al. 2014 (toolbox paper)**  
   - Implements algebraic and dynamic quenching with a kinematic alpha effect. Vishniac flux is not included.
   
2. **Sharanya's Notes**  
   - Uses an older expression of Vishniac flux without a kinematic alpha effect.
   
3. **2017 Notes**  
   - Explores various stratification configurations for \( u \) and \( \rho \).

With appropriate adjustments to parameters and flux expressions, this code has the potential to reproduce any results within the **First Order Smoothing Approximation (FOSA)**. 

--- 
