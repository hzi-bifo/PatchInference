# GraphCut-Algorithm

Antipatch

--- install --- To use antipatch your system must have installed make and g++. If you do not want to use the g++ compiler, that is found in your path, you have to change it in the makefile. Otherwise open a console, change to the root path (where you found this README too) and type 'make'.

This will produce a path named AntiPatch/bin, where all compiled src will lay in and the programm 'antipatch'.

The command 'make dist-clean' will delete all exept the program 'antipatch' and the src files.

--- run --- To run antipatch you have the choice to edit the config file or to type the parameter in your system.

- use config: ./antipatch config_surfaceOnly
- use paramter option: ./antipatch -i testdata/SEAS_H1_dnds.pfd -n testresults/surfaceAndBurried -g testdata/testfiles/oneValue.onlyInternal -s testdata/3hmg -f 0.187 -b testdata/3HMG.pdb -r oneValueInternal
please ask questions and reports faults to tkl15@helmholtz-hzi.de
