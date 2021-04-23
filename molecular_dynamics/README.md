# Base for MD scripts configured for AMBER in an HPC env

## Molecular dynamics requires the preparation of a strucutre with force parameters for it and its ligand. 
If you have a pdb structure of interest, and the prep files required to add parameters, begin with the below

Steps for Mutations

1.) Modify the mutate_var_MUT.sh file with positions of variants to chain and AA

2.) Run process_structure.sh on PDB_to_mutate.pdb

3.) Then run the vacuo prep - prep_topology_vacuo.sh and send the prmtop, rst7 and relax.pbs files to chpc

4.) Run the relax.pbs, and return the rst file to base pc

5.) Run prep_topology_wSolvent.sh - it will make the files ready for the start of min, heating, then MD (same script)

6.) Send the startMD.prmtop and rst files, along with protein_MUT_md.pbs to chpc

7.) Once finished, with no issue, that will make the first 10ns - now send the md_MUT_35ns.pbs to chpc - and run it. 

8.) Repeat at 9hr intervals for 25ns steps of MD

9.) Combine the stripped trajectory runs with combineTraj.pbs - scp the agr to look at
