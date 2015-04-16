
outDir="testresults"
mkdir -p $outDir
siteAss=AntiPatch/tools/siteAssignment.py
outResidues=$outDir/allSelectedResidues.dat
outPatchList=$outDir/allPatches.dat

# get antipatch
if [ ! -f antipatch ]; then make; fi;

treeMode="intExt"
normMode="allVals"
surfMode="surfcOnly"
surfacedata="testdata/3hmg"
pdb="testdata/3HMG.pdb"
onlyPyList="-k false -m false -y false"
inFile="testdata/oneValue.onlyInternal"
stop="1.1"

surfOption="false"
name="surfaceAndBurried"
echo "> $name"
./antipatch -n $outDir/$name -g $inFile -z $stop -r $name -u $surfOption -s $surfacedata -b $pdb -c "A"  -p 1 -k true
echo -e "\n\n--------------------------------------\n\t\t\t$name\n--------------------------------------" >> $outPatchList


