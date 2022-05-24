# Set name of GROMACS executable (gmx/gmx_mpi)
GMX=gmx_mpi
# Set number of threads - adjust if needed
# if using SLURM this can be set automatically using the variable SLURM_CPUS_PER_TASK
NT=10
# Set MPI launcher (srun/mpirun) - you might need to specify the number of MPI tasks
MPI=srun

# 1) extract starting conformation from NVT runs [conf_start.gro].
#    Choose different values of X [in ps] to obtain different starting conformations for the 3 production runs.
#    In the paper, we chose X=5000,7500,10000
echo 0 | $GMX trjconv -f ../../1-EQUIL/traj.trr -s ../../1-EQUIL/topol.tpr -dump X -o conf_start.gro

# 2) creation of tpr file - run once before first production run
$GMX grompp -f ../../../DATA/gromacs-mdps-charmm/3-nvt-production.mdp -c conf_start.gro -n ../../0-TOPO/index.ndx -p ../../0-TOPO/topol.top

# 3) production in chunks of 24h (adjust based on available resources) - resubmit the line below to restart the simulation from the state (.cpt) file
$MPI $GMX mdrun -ntomp $NT -maxh 24.0 -cpi -noappend
