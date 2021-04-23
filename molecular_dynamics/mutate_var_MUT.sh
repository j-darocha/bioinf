### mutate Ser 188 to Phe (Mediterranean variant)
# Only remove the atoms of the 
# side chain and tleap will correct them

#This PDB has the PKA residues and WT variant made
cp 6E07_to_mutate.pdb temp.pdb

#Replace these to match variant in title
sed -i '1313,1314d' temp.pdb 
sed -i 's/SER A 188/PHE A 188/' temp.pdb

#cp temp.pdb temp1sdel.pdb

sed -i '5252,5253d' temp.pdb
sed -i 's/SER B 188/PHE B 188/' temp.pdb

mv temp.pdb G6PD_MED_prevac.pdb
