# RDTCC
SIMULTANEOUS EEG-FMRI INVESTIGATION OF RHYTHM-DEPENDENT THALAMO-CORTICAL CIRCUITS ALTERATION IN SCHIZOPHRENIA

This readme document describes software codes.
The source codes are provided as service to the scientific community and may be used for any non-commercial purpose.  Users should use the codes at their own risks.


To run these codes, you need to have MATLAB installed on your computer (the current code has been tested on MATLAB 2022b and 2023b). 

You can download and deploy EEGLAB toolbox from https://sccn.ucsd.edu/eeglab/index.php. Please use latest version of EEGLAB.

Function called "ADTF" and (its nested functions) are used in these codes which is borrowed from eConnectome (https://www.nitrc.org/projects/econnectome/). 
We modified the computation flow of the code to reduce the time and compute resource requirements.
Please note that this is a personalized adjustment to the current calculation. If you want to be able to adjust all possible parameters, use the original code in eConnectome. 
We included this function and (its nested functions) in Codes folder, for your convenience. eConnectome is a MATLAB toolbox that can be installed easily and freely. For the purposes of running our codes, its installation is not needed.

