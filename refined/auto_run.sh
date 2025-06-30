#!/bin/bash

cd run_files
trial=525

# for xi in $(seq 0.01 0.01 3.00); do
for xi in 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5 5.5 6.0; do
# for xi in 4.5 5 5.5 6.0; do
    echo "Running trial $trial with xi_0 = $xi"

    # Modify xi_0 in parameters_constants.f90
    sed -i "s/^\s*double precision :: xi_0 = .*/      double precision :: xi_0 = $xi \!control of the small scale magnetic field strength ON VISHNIAC FLUX EXPRESSION/" ../codes/parameters_constants.f90

    # Compile
    start=$(date +%s)
    gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/temporal_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90

    # Run
    ./run
    end=$(date +%s)
    echo "Execution time for trial $trial: $(expr $end - $start) seconds"

    # Run Python script
    python3 ../codes/plot.py --trial $trial

    # Increment trial counter
    trial=$((trial + 1))
done

cd ..
