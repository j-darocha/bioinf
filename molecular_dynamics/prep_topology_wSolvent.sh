##make pdb from vacminned rst and the prmtop 
MUT=MEDI

ambpdb -p G6PD_${MUT}_vacuotomin.prmtop -c G6PD_${MUT}_vacminned.rst > G6PD_${MUT}_vacminned.pdb
## generate parm and top file with solvent post in_vacuo_minimisation

MUTPDB=G6PD_${MUT}_vacminned.pdb


echo "
source ~/amber18/dat/leap/cmd/oldff/leaprc.ff14SB       
source ~/amber18/dat/leap/cmd/leaprc.gaff                                                               
loadamberparams frcmod.ionsjc_tip3p
                                                       
loadAmberPrep /home/jorge/Desktop/G6PD/Structures/mol_dynamics/G6PD/gen_lib_nad/nad.prep
loadAmberparams /home/jorge/Desktop/G6PD/Structures/mol_dynamics/G6PD/gen_lib_nad/nad.frcmod 

gnad = loadpdb $MUTPDB
         
#add ions and neutralizing the system
addions  gnad  Cl- 0 			
solvateoct gnad  TIP3PBOX 15.0

savepdb gnad G6PD_nad_startMD.pdb                                                       
saveamberparm gnad G6PD_${MUT}_startMD.prmtop G6PD_${MUT}_startMD.rst7
quit                                                                             
" >prepare_complex.in                                                           

tleap -f prepare_complex.in 


