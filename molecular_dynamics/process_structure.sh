# download the structure PDB file 

wget https://files.rcsb.org/view/6E07.pdb

# add chains Q and T with NADP and rename to chain A and B
convpdb.pl -chain Q -setchain A -out generic 6E07.pdb |grep -v END >temp.pdb
convpdb.pl -chain T -setchain B -out generic 6E07.pdb |grep -v END >>temp.pdb
grep NAP  6E07.pdb | grep HETATM | grep " Q " |sed 's/ Q / A /'>>temp.pdb
echo TER >>temp.pdb
grep NAP  6E07.pdb | grep HETATM | grep " T "|sed 's/ T / B /' >>temp.pdb
echo TER >>temp.pdb

### mutate L459 to Arg(WT) 
# I only remove the atoms of the 
# side chain and tleap will correct for
# msiign atoms later
sed -i '3495,3497d' temp.pdb 
sed -i 's/LEU A 459/ARG A 459/' temp.pdb

sed -i '7436,7438d' temp.pdb
sed -i 's/LEU B 459/ARG B 459/' temp.pdb

mv temp.pdb 6E07_clean_wt.pdb 



























