# Harpseal
A population model for Canadian harp seals, estimated using a Bayesian hierarchical model

## Notes:
A summary of methods and results from the latest analysis can be found in the latest html report "Report_HS_Results..." in the results sub-folder. The Bayesian model itself (coded in STAN language) is contained in the file "HSmodfit.stan", which is fit to data by running "HS_modfit_shell.r". 

As new data become available, the data files contained in the /data sub-folder should be updated accordingly, and the model re-fit to the new data. 
The steps are as follows:
1. Update all data spreadsheets

2. Run "HS_modfit_shell.r" to fit the model - this takes about an hour, depending on how many Cores (parallel processors) are available. 
  The results will be saved to an rdata file in the /Results sub-folder that has date_time in the file name 
  (as well as some codes describing age range and bias corrction value used for age strucuture). 

3. To inspect the results from a given rdata results file, run the "HS_modfit_plot.R" script. The script prompts you to select an "HS_Results..." rdata file, or uses the one already loaded in environment, and creates a number of plots.

4. Once you are satisfied with model fitting results, run the script "HS_estTAC.r" to estimate K and TAC values via Monte Carlo simulations. This script will first prompt you to select an "HS_Results..." rdata file, or will use the one already loaded in environment. Once finished you can save the TAC results as prompted.

5. Once a TAC analysis has been completed, you can generate an HTML markdown report. To do so, run the script "HS_Make_report.r", which will prompt you to select an "HS_Results..." rdata file (Note that if no matching TAC results file is found, the script will stop). The script then generates an html report (by knitting "HSmodel_summary,rmd"), which is saved in the results folder with a file name that begins with "Report_" and then the name of the results file used.

