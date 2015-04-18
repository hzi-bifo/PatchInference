## PatchInference

This repository contains the source code source code for AntiPatch, a software for inference of antigenicity-altering patches of sites on a protein structure.
For reference, please cite "C. Kratsch, L. Mümken, L. Steinbrück, T.R. Klingen, A.C. McHardy; Determination of antigenicity-altering patches of sites on the major surface protein of human influenza A/H3N2 viruses (in review)".


The source code is located under AntiPatch/src/.

### How to install
To use antipatch your system must have installed make and g++. If you do not want to use the g++ compiler, that is found in your path, you have to change it in the makefile. Otherwise open a console, change to the root path (where you found this README too) and type 'make'.

This will produce a path named AntiPatch/bin, where all compiled src will lay in and the programm 'antipatch'.

The command 'make dist-clean' will delete all exept the program 'antipatch' and the src files.

### How to run
To run antipatch you have the choice to edit the config file or to type the parameter in your system.

- use config: ./antipatch config_surfaceOnly
- use parameter option: ./antipatch -g testdata/oneValue.onlyInternal -z 1.1 -u false -s testdata/3hmg -b testdata/3HMG.pdb -c "A"  -p 1 -k true

Alternatively, you can run the bash script produceResults.sh to infer patches for the reference hemagglutinin structure (PDB identifier 3HMG).
- just type: bash produceResults.sh

With the -h option you can request help to start the program.

### Questions 
Please ask questions and report faults to <tkl15@helmholtz-hzi.de>
