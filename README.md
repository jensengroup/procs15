# ProCS15
ProCS15 is a DFT-based chemical shift predictor for protein backbone atoms and CB as described in https://peerj.com/articles/1344/ .

It is written as a module for the open-source protein simulation framework PHAISTOS (www.phaistos.org).
The newest version of PHAISTOS can be downloaded with
'svn checkout svn://svn.code.sf.net/p/phaistos/code/trunk phaistos'

Then just clone this repository to the modules folder and build.
Data files required for the predictor is available at http://www.erda.dk/public/archives/YXJjaGl2ZS1pMHN3TXE=/procs/ProCSnumpyfiles.tar.bz2

## Prediction from pdb-files
An executable to predict chemical shieldings from a set of pdb-files is included. Note that it is the shieldings and not chemical shifts that are predicted, so they will have to be scaled to the experimental values.
Run `make procs15_predictor` in the build directory to compile the executable, which will afterwards be located in ./build/test/procs15\_predictor.

Running `procs15_predictor --data-folder <path to numpy data files> --pdb-files file1.pdb (file2.pdb ...)` will save the prediction to a file replacing the .pdb extension with .procs.
