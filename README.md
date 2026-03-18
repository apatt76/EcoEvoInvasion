# EcoEvoInvasion

This repo has 3 main parts: a folder for code, a folder for simulation output, and a folder for plots.

The code folder contains 4 code files. First run RandomParSearch.R to generate 10000 random parameter draws and solve for long term population dynamics for each scenario. Then run MainAnalysis.R for all structured parameter searching and plotting. Both of those files reference the custom R functions that are included in SolverFunctions.R and the main analysis file also references PcaFunctions.R using source. The custom functions deal with iterating analyses, managing data outputs, and producing custom figures. 

The output folder contains outputs that are created from RandomParSearch.R and MainAnalysis.R. These take a long time to run, and so major outputs are saved here for later reuse, post-analysis, and plotting.

The plots folder contains plot outputs from MainAnalysis.R. A few figures are not shown here, like the conceptual diagram (drawn on computer without data), and the compound supplemental figures. Compound were assembled in post-processing, but instructions for re-generating each piece are found in the MainAnalysis.R.