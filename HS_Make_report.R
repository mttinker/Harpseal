rm(list = ls())
require(svDialogs)
require(rJava)
require(rChoiceDialogs)
require(rmarkdown)

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
rspnse = dlg_message(c("This script is used to generate an html report that summarizes results of  ",
                       "an age-structure population model fit to harp seal data. ",
                       "You can begin by selecting a Results file. ",
                       "Continue?"), "okcancel")$res
if(rspnse == "cancel"){
  stop_quietly()
}
Sys.sleep(.1)
file_list = list.files(path = "./Results",pattern = "HS_Results",full.names = F)
rslt_list = grep("TAC", file_list, value = TRUE, invert = TRUE)
rslt_list = grep("Report", rslt_list, value = TRUE, invert = TRUE)
rdata_file = rselect.list(rslt_list, preselect = NULL, multiple = FALSE,
                          title = "Select results file" ,
                          graphics = getOption("menu.graphics")) 
if(length(rdata_file)==0){
  dlg_message(c("No data file selected"), "ok")
  stop_quietly()
}
TAC_list = grep("TAC", file_list, value = TRUE, invert = F)
TAC_file = TAC_list[which(TAC_list==paste0("TAC_",rdata_file))]
if(length(TAC_file)==0){
  dlg_message(c("No TAC sum for this results file: run 'HS_estTAC.r' script first"), "ok")
  stop_quietly()
}else{
  file.copy(paste0("./Results/",rdata_file),"Results.rdata",overwrite = TRUE)
  file.copy(paste0("./Results/",TAC_file),"TAC_K_est.rdata",overwrite = TRUE)
}
resultsfilename = substr(rdata_file,1,nchar(rdata_file)-6)

YearT = paste0("20",substr(resultsfilename,12,13))
vers = read.csv("./data/Version.csv"); vers = vers$Version
title = paste0("Harp Seal Population Model, ", YearT)
subtitle = paste0("Model ", vers,", Results file: ", resultsfilename)
Daterun = Sys.Date()
render("./HSmodel_summary.Rmd",
       output_dir = "./Results",
       output_file = paste0("Report_",resultsfilename,".html"),
       params = list(rep_title = title, rep_subtitle = subtitle, rep_date = Daterun)) # 
dlg_message(c("The results can be viewed by opening the approproiate 'Report' html file in the Results folder"),
            "ok")
