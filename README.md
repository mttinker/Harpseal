# Harpseal
Population model for Canadian harp seals

## Notes:
1. A summary of the methods and results from the latest analysis can be found in the report "HSmodel_summary.html". 
2. This report is generated from the R markdown file "HSmodel_summary.Rmd", and is based on the most recent model fitting results (which are located in the rdata file "Results.rdata")
3. The Bayesian model itself (coded in STAN language) is contained in the file "HSmodfit.stan". The model is fit to data by executing the R script "HS_modfit_shell.R". 
4. As new data become available, the data files contained in the /data sub-folder should be updated accordingly, and the model re-fit to the new data by running "HS_modfit_shell.r". Model fitting takes from 1-3 hours, depending on how many Cores (parallel processors) are available. The results will be saved to an rdata file in the /Results sub-folder that has date/time in the file name. To inspect the results, run the "HS_modfit_plot.R" script, which will call up and then plot results and run simulations.
5. Once you are satisfied with model fitting results, the associated rdata file should be copied and pasted into the root directory with the name "Results.rdata", thus over-writing the existing file. The R markdown file "HSmodel_summary.Rmd" can then be re-knit to generate a new report.

