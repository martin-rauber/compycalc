####################################################################################
#COMPYCALC: COMprehensive Yield CALCulation 
####################################################################################
#Version information: 1.0.7
####################################################################################
#USER: follow the four instructions and run the script
####################################################################################

#set wd-----------------------------------------------------------------------------

# 1)  set the working directory (wd) for COMPYCALC: this folder must contain the compycalc R
#     file and the folder /src containing the scripts. The wd name will be used to name the 
#     result files. 

#setwd("")

#add data---------------------------------------------------------------------------

# 2)  Sunset data
#     Add your folder(s) with the individual Sunset measurement(s) to the wd folder. 

# 3)  enter F14C raw data
          ##import from csv
          F14C_raw_data = read.csv("MR01-143-EC-F14C-raw-data.csv", header = TRUE)
          F14C_raw_data = F14C_raw_data[,]
          #OR enter here manually
          #F14C_raw_data = c(0.8812040,0.5968259,0.6417999,0.6892420,0.5443450)

# 4)  run script
          
#OUTPUT          
#         - each folder with measurements will get three files: "last-digit-of-folder"-mean-results.csv, "last-digit-of-folder"-raw-results.csv, "last-digit-of-folder"-stats.csv
#         - wd folder will get "your-wd-name-here"-mean-summary-with-F14C.csv and "your-wd-name-here"-F14C_and_EC-yield-and-charring-summary.pdf
             
####################################################################################
# prep.
####################################################################################
#load libraries & install packages if necessary-------------------------------------
load.lib = c("ggpubr", "dplyr", "data.table", "purrr", "stringr", "pastecs","readxl","MASS", "ggplot2", "pracma", "signal")   
install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)
#clean up environment---------------------------------------------------------------
rm(list=setdiff(ls(), c("result_filename", "csv_raw", "csv_stat", "csv_mean","fitting_type","manual.coef", "r_scripts", "home_wd","F14C_raw_data")))
if(!is.null(dev.list())) dev.off()
#save home wd-----------------------------------------------------------------------
home_wd = getwd()

####################################################################################
# yield calc.
####################################################################################

#run yield calc for each subfolder--------------------------------------------------
parent_folder = getwd()
sub_folders = list.dirs(parent_folder, recursive=TRUE)[-1]
r_scripts <- file.path(parent_folder, "src/yields_calc_io.R")
# Run scripts in sub-folders 
for(j in sub_folders) {
  setwd(j)
  source(r_scripts)
}

####################################################################################
# generate summary boxplots for EC-yield and charring from all measurements
####################################################################################

#cd to parent directory-------------------------------------------------------------
setwd(home_wd)
getwd()

#load csv files from all subdirectories---------------------------------------------
#raw results from each filter
df = list.files(".",pattern = "*raw-results.*csv",  recursive = TRUE) %>% map_df(~fread(.))
colnames(df) = c("filtercounter","EC_yield", "charring_S1", "charring_S2", "charring_S3", "filter_name")
df$filter_name_short=str_sub(df$filter_name,(nchar(df$filter_name)),nchar(df$filter_name))

df_charring = df[,3:7]
df_charring

#stats file
df_stats = list.files(".",pattern = "*-stats.*csv",   recursive = TRUE) %>% map_df(~fread(.))
df_stats$filter_name_short = c(rep(unique(df$filter_name_short),each = 14))
df_stats
# number of filters for each sample
df_stats[is.element(df_stats$V1, "nbr.val"),]
# mean for each  sample and export to csv with partent folder prefix
df_stats_mean=df_stats[is.element(df_stats$V1, "mean"),]
#save csv with EC-yield and charring
#write.csv(df_stats_mean, file =  paste(basename(getwd()),"-mean-summary.csv",sep=""), row.names = F)

#number labs
n_labs <- df_stats[is.element(df_stats$V1, "nbr.val"),] 
n_labs$label = paste0("n=" ,n_labs$EC_yield)

n_labs

#generate plots--------------------------------------------------------------------

#plot: EC-yield
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_EC_yield = ggplot(df, aes(x=filter_name_short, y=EC_yield)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
  xlab("filter")+
  ylab("EC-yield")
plot_EC_yield = plot_EC_yield + theme( plot.margin = margin(1, 0.2, 0.2, 0.2, "cm"))
plot_EC_yield

#plot: charring S1
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S1 = ggplot(df_charring, aes(x=filter_name_short, y=charring_S1)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
  xlab("")+
  ylab("")+
  geom_hline(yintercept = 0, color="red", lty= 5)+
  geom_text(aes(filter_name_short, y = Inf, label = label), data = n_labs, vjust = 1)  +
  coord_cartesian(ylim = c(-0.05, 0.14))
plot_charring_S1 = plot_charring_S1 + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))
plot_charring_S1 

#plot: charring S2
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S2 = ggplot(df_charring, aes(x=filter_name_short, y=charring_S2)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
  xlab("")+
  ylab("")+
  geom_hline(yintercept = 0, color="red", lty= 5)+
  coord_cartesian(ylim = c(-0.05, 0.14))
plot_charring_S2 = plot_charring_S2 + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))
plot_charring_S2

#plot: charring S3
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S3 = ggplot(df_charring, aes(x=filter_name_short, y=charring_S3)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
  xlab("filter")+
  ylab("")+
  geom_hline(yintercept = 0, color="red", lty= 5)+
  coord_cartesian(ylim = c(-0.05, 0.14))
plot_charring_S3 = plot_charring_S3 + theme( plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))

#create a summary figure for charring
fig_summary_charring = ggarrange(
  plot_charring_S1,
  plot_charring_S2, 
  plot_charring_S3, 
  labels = c("S1", "S2", "S3"), 
  align = "hv" ,ncol = 1, 
  nrow = 3)
fig_summary_charring = annotate_figure(fig_summary_charring,left = text_grob("charring", color = "black", size = 14, rot = 90))
fig_summary_charring

#create a summary figure with EC-yield and charring. Export plots as pdf 

#  pdf(file =  paste(basename(getwd()),"-EC-yield-and-charring-summary-boxplot.pdf",sep=""), width = 8.3 , height = 11.7)
# fig_summary = ggarrange(
#   plot_EC_yield,
#   fig_summary_charring,
#   ncol = 1,
#   nrow = 2)
# annotate_figure(fig_summary, top = text_grob(paste("\n",basename(getwd())," EC-yield and charring boxplots",sep=""), color = "black" , face = "bold", size = 16), bottom = text_grob(paste(Sys.info()[["user"]],Sys.time(), "  ",sep=" "), color = "black", hjust = 1, x = 1, face = "italic", size = 10),)
#  dev.off()
 
####################################################################################
#correction to 100% EC-yield
####################################################################################
 
#get EC-yield data
EC_yield_mean_summary = as.data.frame(df_stats_mean[,2])
EC_yield_mean_summary = EC_yield_mean_summary[,]
#import EC F14C raw data: top of file

#Correction of F14C to 100% EC-yield
source("src/corr_100_EC.R")
#export result as csv
F14C_EC100 = F0.all
df_stats_mean_F14C = cbind(df_stats_mean, F14C_EC100)
write.csv(df_stats_mean_F14C, file =  paste(basename(getwd()),"-mean-summary-with-F14C.csv",sep=""), row.names = F)

#plot raw vs corrected F14C values for each filter
F14C_raw = as.data.frame(F14C_raw_data)
F14C_raw$filter_name_short = df_stats_mean_F14C$filter_name_short
F14C_raw$class = "uncorrected"
F14C_raw = F14C_raw[, c("filter_name_short","F14C_raw_data", "class")]
colnames(F14C_raw) = c("filter_name_short","F14C", "class") 
F14C_raw
F14C_corr = df_stats_mean_F14C[,6:7]
F14C_corr$class = "corr. to 100% EC-yield"
colnames(F14C_corr) = c("filter_name_short","F14C", "class") 
F14C_corr
all_F14C_data = rbind(F14C_raw,F14C_corr)
all_F14C_data

#plot the corrected and uncorrected F14C value for each filter. If loop to adjust y axis for supermodern values
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_all_F14C = ggplot(all_F14C_data, aes(x=filter_name_short, y=F14C, color = class )) + 
  geom_point()+
  xlab("filter")+
  ylab(bquote(~ F^14~C))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,
                                                  if ( max(all_F14C_data$F14C) > 1) {
                                                    max(all_F14C_data$F14C)+0.3
                                                  } else {
                                                    1
                                                  }
                                                  
                                                  ))
plot_all_F14C = plot_all_F14C + theme( plot.margin = margin(1, 0.2, 0.2, 0.2, "cm"),legend.title = element_blank(), legend.position="top")
plot_all_F14C

#plot the corrected and uncorrected F14C values against each other
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_F14C_raw_vs_corr = ggplot() + 
  geom_point(aes(x = F14C_raw$F14C, y = F14C_corr$F14C))+
  geom_abline(intercept = 0, slope = 1)+
  xlab(bquote(~ F^14~C~" uncorrected"))+
  ylab(bquote(~ F^14~C~" corr. to 100% EC-yield"))+
  coord_cartesian(xlim = c(0,1),ylim = c(0,1),expand = FALSE)
plot_F14C_raw_vs_corr = plot_F14C_raw_vs_corr + theme( plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot_F14C_raw_vs_corr

#summary figure with F14C and EC-yiled
fig_summary_F14C_top = ggarrange(
  plot_all_F14C,
  plot_EC_yield,
  ncol = 2,
  nrow = 1)
fig_summary_F14C_top

#Version info in pdf
compycalc_version = "pdf generated by COMPYCALC version 1.0.7"

#create a summary figure with F14C, EC-yield, and charring. Export plots as pdf 
pdf(file =  paste(basename(getwd()),"-F14C_and_EC-yield-and-charring-summary.pdf",sep=""), width = 8.3 , height = 11.7)
fig_summary = ggarrange(
  fig_summary_F14C_top,
  fig_summary_charring,
  ncol = 1,
  nrow = 2)
annotate_figure(fig_summary, top = text_grob(paste("\n",basename(getwd())," summary",sep=""), color = "black" , face = "bold", size = 16),
               bottom = text_grob(paste(compycalc_version, "   ", Sys.info()[["user"]],Sys.time(), "  ",sep=" "), color = "black",  face = "italic", size = 10),)
dev.off()

####################################################################################
#end COMPYCALC
####################################################################################


