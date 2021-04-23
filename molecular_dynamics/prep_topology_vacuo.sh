## generate parm and top file 

##load the mutated pdb
MUTPDB=G6PD_ABEN_prevac.pdb
MUT=ABEN

echo "
source ~/amber18/dat/leap/cmd/oldff/leaprc.ff14SB       
source ~/amber18/dat/leap/cmd/leaprc.gaff                                                               
loadamberparams frcmod.ionsjc_tip3p
                                                       
loadAmberPrep /home/jorge/Desktop/G6PD/Structures/mol_dynamics/G6PD/gen_lib_nad/nad.prep
loadAmberparams /home/jorge/Desktop/G6PD/Structures/mol_dynamics/G6PD/gen_lib_nad/nad.frcmod 

gnad = loadpdb $MUTPDB
         

savepdb gnad G6PD_${MUT}_vacuotomin.pdb                                                       
saveamberparm gnad G6PD_${MUT}_vacuotomin.prmtop G6PD_${MUT}_vacuotomin.rst7
quit                                                                             
" >prepare_complex.in                                                           

tleap -f prepare_complex.in 


