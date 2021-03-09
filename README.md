# Harpseal
A population model for Canadian harp seals, estimated using a Bayesian hierarchical model

## Notes:
A summary of methods and results from the latest analysis can be found in the report "HSmodel_summary.html". This report is generated from the R markdown file "HSmodel_summary.Rmd", and is based on the most recent model fitting results (which are located in the rdata file "Results.rdata"). The Bayesian model itself (coded in STAN language) is contained in the file "HSmodfit.stan". 

As new data become available, the data files contained in the /data sub-folder should be updated accordingly, and the model re-fit to the new data. The steps are as follows:
1. Run "HS_modfit_shell.r" to fit the model - this takes from 1-3 hours, depending on how many Cores (parallel processors) are available. The results will be saved to an rdata file in the /Results sub-folder that has date_time in the file name. 
2. To inspect the results from a given rdata file, run the "HS_modfit_plot.R" script, which prompts you to select an rdata file from the Results sub-folder and creates a number of plots.
3. Once you are satisfied with model fitting results, the associated rdata file should be copied and pasted into the root directory with the name "Results.rdata" (thus over-writing the existing file, which you may wish to archive). 
4. Run the script "HS_estTAC.r" to estimate K and TAC values via Monte Carlo simulations. 
5. Finally, the R markdown file "HSmodel_summary.Rmd" can be re-knit to generate a new html report.

