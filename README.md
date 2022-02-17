# diag-NEMO-PISCES
Complete NEMO-PISCES outputs with a catalogue of diagnostics

To install, enter a terminal an type: 

    git clone git@github.com:acapet/diag-NEMO-PISCES.git
    cd diag-NEMO-PISCES
    conda env create -f environment.yml

This may take a little while.
Note that you will need to load a conda module from your cluster environment beforehand, which depends on your system.

To use

   conda activate PISCES-DIAG
   python diag.py --dir <DIRECTORY>

`<DIRECTORY>` should be a directory including outputs of NEMO-PISCES simulations.
The script will search for any `*ptrc*.nc` files and issue a new file `*diag*.nc` with the same name structure, nc attributes, and containing the newly computed 2D diagnsotics.

Note that the file and variable requirements depends on the list of requested diagnostics.

