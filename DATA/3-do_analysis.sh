# Set name of GROMACS executable (gmx/gmx_mpi)
GMX=gmx_mpi
# Set MPI launcher (srun/mpirun) - you might need to specify the number of MPI tasks
MPI=srun
# Set FoldX executable
FOLDX=~/bin/FoldX/foldx

###
### preparing trajectory
###
# 1) first we need to cat together all the production trajectories: we take a total of 10000 frames
$MPI $GMX trjcat -f ../2-PRODUCTION/traj_comp.part*.xtc -o traj_all.xtc -dt 100
# 2) now we fix PBCs
echo 0 | $MPI $GMX trjconv -f traj-all.xtc -o traj-PBC.xtc -s ../../0-TOPO/topol_xtc.tpr -pbc nojump -e 1000000

# 3) remove temporary trajectory
rm traj_all.xtc

###
### now we do the analysis on the full trajectory
###
# 4) PLUMED post-processing: RMSD calculation
plumed driver --plumed ../../../DATA/plumed_driver.dat --mf_xtc ../2-PRODUCTION/traj-PBC.xtc

# 5) Hydrogen bonds analysis
# execute python script: requires python3 and MDAnalysis
python ../../../DATA/count-HB.py ../../0-TOPO/topol_xtc.tpr ../2-PRODUCTION/traj-PBC.xtc

# 6) Salt bridges analysis
# execute python script: requires python3 and MDAnalysis
python ../../../DATA/count-SB.py ../../1-EQUIL/PDBs/conf_emin.pdb ../2-PRODUCTION/traj-PBC.xtc 

# 7) Create individual PDBs from traj-PBC.xtc
# you might want to parallelize this using a job-array
mkdir PDBs; cd PDBs
# execute python script: requires python3 and MDTraj
python ../../../../DATA/get_PDB_frames.py ../../../1-EQUIL/PDBs/conf_emin.pdb ../../2-PRODUCTION/traj-PBC.xtc  0 10001
# done
cd ../

# 8) FoldX scoring
mkdir FOLDX; cd FOLDX
# you might want to parallelize this loop using a job-array
for j in `seq 0 10000`
do
 # padding for PDB filename
 jj=$(printf "%05d" $j)
 # run FoldX
 $FOLDX --command=AnalyseComplex --pdb=frame_${jj}.pdb --analyseComplexChains=A,B --complexWithDNA=false --output-file=out_$jj clean-mode --pdb-dir="../PDBs/" --output-dir="./"
done
# collect "Interaction Energy" of all frames in FoldX.dat
grep -A2 "Interaction Energy" Summary_out_* | grep PDB | awk '{print $6}' > FoldX.dat
# done
cd ../
