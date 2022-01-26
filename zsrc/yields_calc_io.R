###################################################################################
#yields calc io
###################################################################################
#README
#The following files need to be in the same folder and set to the working directory:
      #this file
      #yields_calc_ext.R
      #cooldown_data.csv
      #all Sunset txt files you want to analyse (but only these!)


#DISCLAIMER
#only fitting_type = "poly" works, expo and manual does not

####################################################################################
#user: enter the necessary information
####################################################################################

#select the fitting type values: 
fitting_type = "poly"                                     #poly, expo, or manual
manual.coef = c(7645, -0.371, 5.32E-4)                    #enter manual coefficients if fitting_type manual is selected

#enter a filename for the export files
result_filename = basename(getwd())  

#set value for outlier removal (1.5*IQR recomended)
IQR <- 1.5

####################################################################################
#script
####################################################################################

#clean up environment---------------------------------------------------------------
rm(list=setdiff(ls(), c("result_filename","fitting_type","manual.coef", "r_scripts", "home_wd","F14C_EC_raw_data","F14C_OC_raw_data","F14C_EC_u","charr_corr_slope","IQR" )))
if(!is.null(dev.list())) dev.off()
#load function----------------------------------------------------------------------
data_load_func = function(filename) {
  cooldown = read.csv("../zsrc/cooldown_data.csv", sep = ",", header = T)
  dat = as.data.frame(read.csv(file = filename, sep = ",", skip = 28, header = T ))[,c(1:18)]
  tabla_complete <<- rbind(dat, cooldown)
  yield_calc = function(tabla_complete, fitting_type, manual.coef) {source("../zsrc/yields_calc_ext.R")}
  yield_calc(tabla_complete)
  
  fig_temp <- ggarrange( plot_TwoSide, plot_Sunset, labels = c("Sunset laser", "Sunset protocol"),  nrow = 2, ncol = 1)
  fig_temp <- annotate_figure(fig_temp, bottom = text_grob(paste(basename(getwd())," ",filename, " ", Sys.info()[["user"]],Sys.time(), "  ",sep=" "), color = "black",  face = "italic", size = 10),)
  ggsave(filename = paste(filename,"_temp.pdf",sep=""), plot = fig_temp , width = 8.3, height = 11.7)
  
} 
#
#load data, run calculation --------------------------------------------------------
filename <- dir(".",pattern="^(.*)txt$")
df_raw <- NULL
for (i in filename){
  data_load_func(i)
  df_raw <- rbind(df_raw, data.frame(tabla_resultados2$EC_yield_S3,tabla_resultados2$charring_S1,tabla_resultados2$charring_S2,tabla_resultados2$charring_S3))}

#rm(list=setdiff(ls(), c("df","result_filename", "csv_raw", "csv_stat", "csv_mean")))
df_raw$filter <- c(rep(result_filename, length(df_raw$tabla_resultados2.EC_yield)))
colnames(df_raw) <- c("EC_yield", "charring_S1", "charring_S2", "charring_S3", "filter_name")

#extract specific data--------------------------------------------------------------
df_length_raw <- length(df_raw$EC_yield)
df_charr_raw <- data.frame("charring value" = c(df_raw$charring_S1,df_raw$charring_S2,df_raw$charring_S3), "Sunset step" = c(rep("S1",df_length_raw), rep("S2",df_length_raw), rep("S3",df_length_raw)), stringsAsFactors = FALSE)
df_yield_raw <- data.frame("EC-yield" = df_raw$EC_yield, "dummy"= rep("",df_length_raw), stringsAsFactors = FALSE)

#plots with raw data----------------------------------------------------------------

#plot: EC-yield (raw data)
theme_set(theme_classic(base_size = 13))
plot_yield_raw <- ggboxplot(width=0.33, data=df_yield_raw, x = "dummy", y = "EC.yield",  xlab="", ylab = "EC-yield")
#plot: charring for each Sunset step (raw data)
plot_charr_raw <- ggboxplot(data=df_charr_raw, x="Sunset.step", y="charring.value", xlab="Sunset step", ylab = "charring")+
  geom_hline(yintercept = 0, color="red", lty= 5)

fig_EC_yield_charring_raw <- ggarrange(plot_yield_raw, plot_charr_raw, ncol = 2, nrow = 1)
fig_EC_yield_charring_raw <- fig_EC_yield_charring_raw + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
fig_EC_yield_charring_raw <- annotate_figure(fig_EC_yield_charring_raw, top = text_grob("Raw data", color = "darkred", face = "bold", size = 12))

#remove outliers-------------------------------------------------------------------
iqr_scores <- abs(scores(df_raw[,1:4],type="iqr"))
outlier_rows <- unique(which((iqr_scores>IQR) == TRUE, arr.ind=TRUE)[,1])
n_outlier_rows <- length(outlier_rows)
if (n_outlier_rows == 0) {
  df <- df_raw
} else {
  df <- df_raw[-c(outlier_rows),]
}

#extract specific data w/o outliers-------------------------------------------------
df_length <- length(df$EC_yield)
df_charr <- data.frame("charring value" = c(df$charring_S1,df$charring_S2,df$charring_S3), "Sunset step" = c(rep("S1",df_length), rep("S2",df_length), rep("S3",df_length)), stringsAsFactors = FALSE)
df_yield <- data.frame("EC-yield" = df$EC_yield, "dummy"= rep("",df_length), stringsAsFactors = FALSE)

#calculate stats and generate plots-------------------------------------------------

#general stats
stat <- stat.desc(df[, -5])
#mean of each EC-yield and each Sunset step
df_mean <- stat[9,]

#calculate stats and generate plots-------------------------------------------------

#EC-yield and charring plots after outlier removal

#plot: EC-yield
theme_set(theme_classic(base_size = 13))
plot_yield <- ggboxplot(width=0.33, data=df_yield, x = "dummy", y = "EC.yield",  xlab="", ylab = "EC-yield")
#plot: charring for each Sunset step
plot_charr <- ggboxplot(data=df_charr, x="Sunset.step", y="charring.value", xlab="Sunset step", ylab = "charring")+
  geom_hline(yintercept = 0, color="red", lty= 5)

fig_EC_yield_charring_wo_outliers <- ggarrange(plot_yield, plot_charr, ncol = 2, nrow = 1)
fig_EC_yield_charring_wo_outliers <- fig_EC_yield_charring_wo_outliers + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
fig_EC_yield_charring_wo_outliers <- annotate_figure(fig_EC_yield_charring_wo_outliers, top = text_grob("Outliers removed", color = "darkblue", face = "bold", size = 12))

#stats table
colnames(stat) <- c("EC-yield", "charring S1", "charring S2", "charring S3")
stat_table <- ggtexttable(format(stat[c(9,4:6,8,12,13),], digits = 3),  theme = ttheme("blank"))

#generate summary pdfs---------------------------------------------------------------

#summary pdf page 1
fig_summary_temp <- ggarrange( plot_fit, fig_EC_yield_charring_raw, fig_EC_yield_charring_wo_outliers, stat_table, nrow = 4, ncol = 1)
fig_summary_temp <- annotate_figure(fig_summary_temp, top = text_grob(paste("\n","Summary:","\n","slope, EC-yield, charring, and statistics",sep=""), color = "black" , face = "bold", size = 16),)
ggsave(filename = "aaa1_temp.pdf", plot = fig_summary_temp , width = 8.3, height = 11.7)

#summary QQ-Plots-------------------------------------------------------------------

#QQ-plots EC
plot_qq_EC_raw <- ggqqplot(df_raw$EC_yield, ylab = "EC-yield", title = "Raw data")
plot_qq_EC_corr <- ggqqplot(df$EC_yield, ylab = "EC-yield", title = "Outliers removed")
#QQ-plots charring S1
plot_qq_S1_raw <- ggqqplot(df_raw$charring_S1, ylab = "Charring S1")
plot_qq_S1_corr <- ggqqplot(df$charring_S1, ylab = "Charring S1")
#QQ-plots charring S2
plot_qq_S2_raw <- ggqqplot(df_raw$charring_S2, ylab = "Charring S2")
plot_qq_S2_corr <- ggqqplot(df$charring_S2, ylab = "Charring S2")
#QQ-plots charring S1
plot_qq_S3_raw <- ggqqplot(df_raw$charring_S3, ylab = "Charring S3")
plot_qq_S3_corr <- ggqqplot(df$charring_S3, ylab = "Charring S3")

#summary pdf page 2
fig_summary2_temp <- ggarrange( plot_qq_EC_raw, plot_qq_EC_corr, plot_qq_S1_raw, plot_qq_S1_corr,plot_qq_S2_raw, plot_qq_S2_corr, plot_qq_S3_raw, plot_qq_S3_corr, nrow = 4, ncol = 2, align = "hv")
fig_summary2_temp <- annotate_figure(fig_summary2_temp, top = text_grob(paste("\n","Outliers: QQ-Plots","\n",sep=""), color = "black" , face = "bold", size = 16),bottom = text_grob(paste("\n","Total number of filters: ",df_length_raw,"\n","Outliers removed: ",n_outlier_rows,"\n","Filters used for calculation: ", c(df_length_raw-n_outlier_rows),"\n", sep=""), color = "black" , face = "plain", size = 13),)
ggsave(filename = "aaa2_temp.pdf", plot = fig_summary2_temp , width = 8.3, height = 11.7)

#load all temp .pdf files and create one pdf file
file.list <- paste(getwd(), "/",list.files(getwd(), pattern = "*.pdf"), sep = "")
pdf_combine(file.list, output = paste(basename(getwd()),"_calc_summary_plots.pdf",sep=""))
#remove temporary pdf files
file.list.rem <- paste(getwd(), "/",list.files(getwd(), pattern = "*_temp.pdf"), sep = "")
file.remove(file.list.rem)

#save result as csv-----------------------------------------------------------------

  write.csv(df_raw, file =  paste(result_filename,"-raw-results.csv",sep=""), row.names = T)
  write.csv(df, file =  paste(result_filename,"-clean-results.csv",sep=""), row.names = T)
  write.csv(stat, file =  paste(result_filename,"-stats.csv",sep=""), row.names = T)
  write.csv(df_mean, file =  paste(result_filename,"-mean-results.csv",sep=""), row.names = T)

#end--------------------------------------------------------------------------------

####################################################################################
#changelog

#V1.0.6   Outlier removal IQR >1.5, output plots aaa1 and aaa2 pdf for: qqplots etc. 
#V1.0.5   minor changes
#V1.0.4   generate pdf summary of the main yields_calc_ext pots
#V1.0.3   keep plots from yields_calc_ext.R
#V1.0.2   corrections from Gary 02.04.2020
#V1.0.1   automatic script to load libraries and install package if necessary 
#V1.0.0   code based on dev. version mr20200330a


####################################################################################



