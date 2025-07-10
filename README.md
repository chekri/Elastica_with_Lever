# Instructions
This repository contains the numerical continuation codes used in the article. The numerical continuation is performed using AUTO-07p, a powerful and easy-to-install software tool for continuation and bifurcation problems in ordinary differential equations (ODEs). A copy of this software package along with the instructions is provided.

## AUTO-07p Configuration

AUTO-07p is an easy-to-use software package. To install this software package, unzip the provided "auto-07p-master.zip file" and navigate to its parent folder. Then configure the files by running './configure' in the terminal. After successful configuration, execute 'make all'. Then, add the path of the bin folder to your system path file by appending the line "path_to_auto/bin" to your bash file. The installation can be verified by running 'auto' in the terminal. Once installed, the codes can be run easily through auto scripts (auto.filename). To run these scripts, run auto auto.filename in the terminal window. These scripts make use of input files (`filename.f90`) and constants files (`c.constants`) where the parameters running the continuation can be tuned. The variables in c.constants files and script files can be adjusted as per the requirement. Once, the computation is complete, the solution files are written in a "s.filename" file. These files are analyzed and plotted using Python tools provided in the respective folders. Ensure that the necessary Python dependencies, including numpy and matplotlib, are installed.


## Instructions on using the contination codes

The repository includes two example folders, namely, Elastica and Kirchhoff Rod. 

Elastica – A simplified planar elastic rod model with lower dimensionality.

Kirchhoff Rod – A 14-dimensional model representing the full Kirchhoff rod equations. It can be employed for solving more complex problems

Despite the difference in complexity, both models yield the same physical results, serving as a validation of the simplified approach.

Navigate to any of the subfolders, say for Example by typing 

`cd Elastica/Example1_Rotating_Arm` 

in the terminal  and run the appropriate auto file by running the command

`auto auto.filename`. 

### Reference

Please cite or refer to the associated article when using or referencing this repository:
Dhanakoti, Siva Prasad Chakri. “Stability analysis through folds: An end-loaded elastic with a lever arm.” (2025). https://doi.org/10.48550/arXiv.2501.04729



