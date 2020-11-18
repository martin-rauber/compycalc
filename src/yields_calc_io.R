###################################################################################
#yields calc io
###################################################################################

#select the fitting type values: 
fitting_type = "poly"                                     #poly, expo, or manual
manual.coef = c(7645, -0.371, 5.32E-4)                    #enter manual coefficients if fitting_type manual is selected

#would you like to export the resuls as an .csv file (TRUE or FALSE)?
csv_raw=T                                                 #all EC-yields and charring data from each filter
csv_stat=T                                                #statistics: min, max, mean, median etc.
csv_mean=T                                                #export only the mean values 
#enter a filename for the export files
result_filename = basename(getwd())  

####################################################################################
#script
####################################################################################

#load libraries & install packages if necessary-------------------------------------
load.lib = c("ggpubr", "pastecs", "dplyr", "ggplot2", "pracma", "signal","gridGraphics", "pdftools")   
install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

#clean up environment---------------------------------------------------------------
rm(list=setdiff(ls(), c("result_filename", "csv_raw", "csv_stat", "csv_mean","fitting_type","manual.coef", "r_scripts", "home_wd","F14C_EC_raw_data","F14C_OC_raw_data","charr_corr_slope" )))
if(!is.null(dev.list())) dev.off()

#load function----------------------------------------------------------------------

data_load_func = function(filename) {
  cooldown = read.csv("../src/cooldown_data.csv", sep = ",", header = T)
  dat = as.data.frame(read.csv(file = filename, sep = ",", skip = 28, header = T ))[,c(1:18)]
  tabla_complete <<- rbind(dat, cooldown)
  yield_calc = function(tabla_complete, fitting_type, manual.coef) {source("../src/yields_calc_ext.R")}
  yield_calc(tabla_complete)
  fig_temp <- ggarrange( plot_fit, plot_TwoSide, plot_Sunset, labels = c("A", "B", "C"),  nrow = 3, ncol = 1)
  fig_temp <- annotate_figure(fig_temp, bottom = text_grob(paste(basename(getwd())," ",filename, " ", Sys.info()[["user"]],Sys.time(), "  ",sep=" "), color = "black",  face = "italic", size = 10),)
  ggsave(filename = paste(filename,"_temp.pdf",sep=""), plot = fig_temp , width = 8.3, height = 11.7)
}

#load data, run calculation --------------------------------------------------------
filename = dir(".",pattern="^(.*)txt$")
df = NULL
for (i in filename){
  data_load_func(i)
  df = rbind(df, data.frame(tabla_resultados2$EC_yield,tabla_resultados2$charringS1,tabla_resultados2$charringS2,tabla_resultados2$charringS3))
}

#load all .pdf files and create one pdf file
file.list <- paste(getwd(), "/",list.files(getwd(), pattern = "*.pdf"), sep = "")
pdf_combine(file.list, output = paste(basename(getwd()),"_calc_summary_plots.pdf",sep=""))
#remove temporary pdf files
file.list.rem <- paste(getwd(), "/",list.files(getwd(), pattern = "*_temp.pdf"), sep = "")
file.remove(file.list.rem)

if(!is.null(dev.list())) dev.off()
#rm(list=setdiff(ls(), c("df","result_filename", "csv_raw", "csv_stat", "csv_mean", "r_scripts", "home_wd","F14C_EC_raw_data")))
df$filter = c(rep(result_filename, length(df$tabla_resultados2.EC_yield)))
colnames(df) = c("EC_yield", "charring_S1", "charring_S2", "charring_S3", "filter_name")
#extract specific data--------------------------------------------------------------
df_length = length(df$EC_yield)
df_charr = data.frame("charring value" = c(df$charring_S1,df$charring_S2,df$charring_S3), "Sunset step" = c(rep("S1",df_length), rep("S2",df_length), rep("S3",df_length)), stringsAsFactors = FALSE)
df_yield = data.frame("EC-yield" = df$EC_yield, "dummy"= rep("",df_length), stringsAsFactors = FALSE)

#calculate stats and generate plots-------------------------------------------------

#general stats
stat = stat.desc(df[, -5])
#mean of each EC-yield and each Sunset step
df_mean = stat[9,]

#plot: EC-yield
# theme_set(theme_classic(base_size = 13))
# plot_yield = ggboxplot(data=df_yield, x = "dummy", y = "EC.yield",  xlab="", ylab = "EC-yield")
# plot_yield

#plot: charring for each Sunset step
# plot_charr = ggboxplot(data=df_charr, x="Sunset.step", y="charring.value", xlab="Sunset step", ylab = "charring")+
#   geom_hline(yintercept = 0, color="red", lty= 5)
# plot_charr 
# 
# fig = ggarrange(plot_yield, plot_charr, labels = c("A", "B"),ncol = 2, nrow = 1)
# fig

#save result as csv-----------------------------------------------------------------
if ( csv_raw ) {
  write.csv(df, file =  paste(result_filename,"-raw-results.csv",sep=""), row.names = T)
}
if ( csv_stat) {
  write.csv(stat, file =  paste(result_filename,"-stats.csv",sep=""), row.names = T)
}
if ( csv_mean) {
  write.csv(df_mean, file =  paste(result_filename,"-mean-results.csv",sep=""), row.names = T)
} 

#end--------------------------------------------------------------------------------


