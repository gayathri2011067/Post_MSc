cd run_files
start=`date +%s`
gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/spatial_derivatives.f90 ../codes/temporal_derivatives.f90 ../codes/eta_profile.f90 ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 ../codes/run_file.f90
./run
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
python3 ../codes/plot.py --trial 187 
cd ..
# python3 codes/compare_plots.py 600 "alpha_k only" 601 "alpha_k + new diffusion case 1" 602 "alpha_k + new diffusion case 2" 
# python3 codes/compare_plots.py 600 "alpha_k only" 603 "alpha_k + nvf case 1" 604 "alpha_k + nvf case 2" 
# python3 codes/compare_plots.py 600 "alpha_k only" 607 "alpha_k + nvf 0.01 + new diff case 1" 608 "alpha_k + nvf 0.05 + new diff case 1" 609 "0.1" 610 "0.2" 611 "0.3" 612 "0.4" 613 "0.5" 614 "0.6" 615 "0.7" 616 "0.8"


# echo "bash run.sh 0.5 1" | at now
# echo "bash run.sh 1.0 2" | at now + 1 minute
# echo "bash run.sh 1.5 3" | at now + 2 minutes
# echo "bash run.sh 2.0 4" | at now + 3 minutes
# echo "bash run.sh 2.5 5" | at now + 4 minutes
# echo "bash run.sh 3.0 6" | at now + 5 minutes
# echo "bash run.sh 3.5 7" | at now + 6 minutes
# echo "bash run.sh 4.0 8" | at now + 7 minutes


# start=`date +%s`

# xi=$1         # First command line argument
# trial=$2      # Second command line argument

# gfortran -o run ../codes/parameters_constants.f90 ../codes/grid.f90 ../codes/courant_condition.f90 \
# ../codes/spatial_derivatives.f90 ../codes/temporal_derivatives.f90 ../codes/eta_profile.f90 \
# ../codes/omega_profile.f90 ../codes/alpha_profile.f90 ../codes/velocity_profile.f90 \
# ../codes/field_initial.f90 ../codes/seed.f90 ../codes/equations.f90 ../codes/time_stepping.f90 \
# ../codes/run_file.f90

# ./run $xi

# end=`date +%s`
# echo Execution time was `expr $end - $start` seconds.

# python3 ../codes/plot.py --trial $trial

# cd ..
