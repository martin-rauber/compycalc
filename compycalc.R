####################################################################################
####################################################################################
# COMPYCALC: COMprehensive Yield CALCulation 
# A tool for EC yield extrapolation and charring correction 
####################################################################################
#written by Martin Rauber and Gary Salazar
# Github: https://github.com/martin-rauber/compycalc
####################################################################################
#USER:  follow the instructions and run the script
#       details are written in the readme
####################################################################################
####################################################################################

#set wd-----------------------------------------------------------------------------

# 1)  set the working directory (wd) for COMPYCALC: this folder must contain the compycalc.R
#     file and the folder /zsrc containing the scripts. The wd name will be used to name the 
#     result files. Note that no other files should be in the folder, including hidden files.

#setwd("")

#add data---------------------------------------------------------------------------

# 2)  Sunset data
#     Add your folder(s) with the individual Sunset measurement(s) to the wd folder. 

# 3)  enter F14C EC and OC raw data
          ##EC: import from csv
          F14C_EC_raw_data <- read.csv(list.files(".",pattern = "*EC-F14C-raw-data.csv" ,  recursive = FALSE), header = TRUE)
          F14C_EC_u <- F14C_EC_raw_data[,2]
          F14C_EC_raw_data <- F14C_EC_raw_data[,1]
          
          ##OC: import from csv
          F14C_OC_raw_data <- read.csv(list.files(".",pattern = "*OC-F14C-raw-data.csv" ,  recursive = FALSE), header = TRUE)
          F14C_OC_u <- F14C_OC_raw_data[,2]
          F14C_OC_raw_data <- F14C_OC_raw_data[,1]
          
# 4)  run script
          
#OUTPUT          
#         - each folder with measurements will get three files: "last-digit-of-folder"-mean-results.csv, "last-digit-of-folder"-raw-results.csv, "last-digit-of-folder"-stats.csv
#         - wd folder will get "your-wd-name-here"-mean-summary-with-F14C.csv and "your-wd-name-here"-F14C_and_EC-yield-and-charring-summary.pdf

####################################################################################
#USER: DO NOT change the script below this line                      
####################################################################################
####################################################################################
# preparation
####################################################################################
#load libraries & install packages if necessary-------------------------------------
load.lib = c("ggpubr", "dplyr", "data.table", "purrr", "stringr", "pastecs","readxl","MASS", "ggplot2", "pracma", "signal","gridGraphics", "pdftools", "outliers")   
install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)
#clean up environment---------------------------------------------------------------
if(!is.null(dev.list())) dev.off()
#save home wd-----------------------------------------------------------------------
home_wd = getwd()

####################################################################################
# yield calc.
####################################################################################

#run yield calc for each subfolder--------------------------------------------------
parent_folder = getwd()
sub_folders = list.dirs(parent_folder, recursive=TRUE)[-1]
r_scripts <- file.path(parent_folder, "zsrc/yields_calc_io.R")
#omit the hidden .git and .Rproj files
sub_folders <- sub_folders[!str_detect(sub_folders,".git")]
sub_folders <- sub_folders[!str_detect(sub_folders,".Rproj")]
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
#raw results from each filter with outiers for (box)plots
df_all <- list.files(".",pattern = "*raw-results.*csv",  recursive = TRUE) %>% map_df(~fread(.))
colnames(df_all) <- c("filtercounter","EC_yield", "charring_S1", "charring_S2", "charring_S3", "charr_u", "filter_name")
df_all$filter_name_short <- str_sub(df_all$filter_name,(nchar(df_all$filter_name)),nchar(df_all$filter_name))

df_all_charring <- df_all[,c(3:5,7,8)]

#clean results from each filter (w/o outliers)
df <- list.files(".",pattern = "*clean-results.*csv",  recursive = TRUE) %>% map_df(~fread(.))
colnames(df) <- c("filtercounter","EC_yield", "charring_S1", "charring_S2", "charring_S3", "filter_name")
df$filter_name_short <- str_sub(df$filter_name,(nchar(df$filter_name)),nchar(df$filter_name))

df_charring <- df[,c(3:5,7,8)]

#stats file
df_stats = list.files(".",pattern = "*-stats.*csv",   recursive = TRUE) %>% map_df(~fread(.))
df_stats$filter_name_short <- c(rep(unique(df$filter_name_short),each = 14))
df_stats
# number of filters for each sample
df_stats[is.element(df_stats$V1, "nbr.val"),]
# mean for each  sample and export to csv with partent folder prefix
df_stats_mean=df_stats[is.element(df_stats$V1, "mean"),]

#number of filters used for calc
n_labs_calc <- df_stats[is.element(df_stats$V1, "nbr.val"),] 
n_labs_calc$label <- paste0("n=" ,n_labs_calc$`EC-yield`)
n_labs_calc
#total number of filters
n_labs <- df_all %>% group_by(filter_name_short) %>% summarize(count=n())
n_labs$label <- paste0("n=" ,n_labs$count)
n_labs
#number of outliers per filter
n_filter_data <- as.data.frame(cbind(n_labs_calc$filter_name_short,n_labs$count,n_labs_calc$`EC-yield`,n_labs$count-n_labs_calc$`EC-yield`))
colnames(n_filter_data) <- c("Filter name","Total filters","Filters used for calculation","Outliers")
ggtexttable(format(n_filter_data),  theme = ttheme("blank"))

#generate plots--------------------------------------------------------------------

#plot: EC-yield
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_EC_yield = ggplot(df_all, aes(x=filter_name_short, y=EC_yield)) + 
  geom_boxplot(colour = "#424242", outlier.colour="red", outlier.shape=8, outlier.size=2)+
  xlab("Filter")+
  ylab("EC yield")
plot_EC_yield = plot_EC_yield + theme( plot.margin = margin(1, 0.2, 0.2, 0.2, "cm"))
plot_EC_yield

#plot: charring S1
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S1 = ggplot(df_all_charring, aes(x=filter_name_short, y=charring_S1)) + 
  geom_boxplot(colour = "#424242", outlier.colour="red", outlier.shape=8, outlier.size=2)+
  xlab("")+
  ylab("")+
  geom_hline(yintercept = 0, color="red", lty= 5)+
  geom_text(aes(filter_name_short, y = Inf, label = label), data = n_labs , vjust = 1)  +
  coord_cartesian(ylim = c(-0.05, 0.14))
plot_charring_S1 = plot_charring_S1 + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))
plot_charring_S1 

#plot: charring S2
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S2 = ggplot(df_all_charring, aes(x=filter_name_short, y=charring_S2)) + 
  geom_boxplot(colour = "#424242", outlier.colour="red", outlier.shape=8, outlier.size=2)+
  xlab("")+
  ylab("")+
  geom_hline(yintercept = 0, color="red", lty= 5)+
  coord_cartesian(ylim = c(-0.05, 0.14))
plot_charring_S2 = plot_charring_S2 + theme(axis.text.x = element_blank(), plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm"))
plot_charring_S2

#plot: charring S3
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_charring_S3 = ggplot(df_all_charring, aes(x=filter_name_short, y=charring_S3)) + 
  geom_boxplot(colour = "#424242", outlier.colour="red", outlier.shape=8, outlier.size=2)+
  xlab("Filter")+
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
fig_summary_charring = annotate_figure(fig_summary_charring,left = text_grob("Fraction of charring", color = "black", size = 14, rot = 90))
fig_summary_charring

####################################################################################
#correction to 100% EC-yield
####################################################################################
 
#get EC-yield data
EC_yield_mean_summary = as.data.frame(df_stats_mean[,2])
EC_yield_mean_summary = EC_yield_mean_summary[,]

#import EC F14C raw data --> beginning of script

#Extrapolation of F14C to 100% EC-yield
source("zsrc/corr_100_EC.R")

#calculate F14C_EC100_0charr
F14C_EC100 <- F0.all
df_charring$charring_total <- df_charring$charring_S1+df_charring$charring_S2+df_charring$charring_S3
total_charr <- aggregate(x = df_charring$charring_total, by = list(df_charring$filter_name_short), FUN = mean) 
colnames(total_charr) <- c("filter_name_short","total_charr_mean")
#F14C_EC100_0_charr_a: F14C (100% EC yield) with charring (S1+S2+S3) subtraction
F14C_EC100_0_charr_a <- (F14C_EC100-F14C_OC_raw_data*0.5*total_charr$total_charr_mean)/(1-0.5*total_charr$total_charr_mean)
#F14C_EC100_0_charr_b: F14C charring (S1+S2+S3) corrected
F14C_EC100_0_charr_b <- (F14C_EC_raw_data-F14C_OC_raw_data*0.5*total_charr$total_charr_mean)/(1-0.5*total_charr$total_charr_mean)
#F14C_EC100_0_charr_c: F14C with charring (both S1+S2+S3) + Slope correction (charr_corr_slope calculated in yields_calc_ext.R)
F14C_EC100_0_charr_c <- charr_corr_slope*(1-df_stats_mean$`EC-yield`)+F14C_EC100_0_charr_b 
#mean final correction to 0% charring
F14C_EC100_0_charr <- cbind(F14C_EC100_0_charr_a, F14C_EC100_0_charr_c)
F14C_EC100_0_charr <- rowMeans(F14C_EC100_0_charr)
F14C_EC100_0_charr

#estimate uncertainties
# F14C_EC_u: input, measurement uncertainty;  F14C_EC100_u: uncertainty after correction to 100% EC yield; 
# F14C_EC100_0_charr_u uncertainty after correction to 100% EC yield and correction to 0% charring

# F14C_OC_u defined in OC data import
# charr_slope_u defined in yields_calc_ext as fitting_poly.coef_u[2]
#import charr_u
charr_u <- as.data.frame(df_stats_mean[,6]) 
#yield_u approximated with charr_u
yield_u <- charr_u

# uncertainty of the F14C(EC) corrected to 100% EC yield from the corr_100_EC script
F14C_EC100_u<- F0.all_sig
dFcharr_a_FOC<- -(0.5*total_charr$total_charr_mean)/(1-0.5*total_charr$total_charr_mean)
dFcharr_a_FEC100 <- 1/(1-0.5*total_charr$total_charr_mean)    
#uncertainty of the charring a relative to the uncertaintiy of the F(EC) 
dFcharr_a_FEC_raw <- dFcharr_a_FEC100 
dFcharr_a_tcharr <- ((0.5*(F14C_EC100-0.5*F14C_OC_raw_data*total_charr$total_charr_mean))/(1-0.5*total_charr$total_charr_mean)^2)-((0.5*F14C_OC_raw_data)/(1-0.5*total_charr$total_charr_mean))
dFcharr_b_FEC_raw<-dFcharr_a_FEC100
dFcharr_b_tcharr<-((0.5*(F14C_EC_raw_data-0.5*F14C_OC_raw_data*total_charr$total_charr_mean))/(1-0.5*total_charr$total_charr_mean)^2)-((0.5*F14C_OC_raw_data)/(1-0.5*total_charr$total_charr_mean))
dFcharr_c_slope<-(1-df_stats_mean$`EC-yield`) ; dFcharr_c_yield <- -charr_corr_slope
F14C_EC100_0_charr_a_u<-sqrt((dFcharr_a_FOC*F14C_OC_u)^2 + (dFcharr_a_FEC100*F14C_EC100_u)^2 + (dFcharr_a_tcharr*charr_u)^2)
F14C_EC100_0_charr_b_u<-sqrt((dFcharr_a_FOC*F14C_OC_u)^2 + (dFcharr_a_FEC_raw*F14C_EC_u)^2   + (dFcharr_b_tcharr*charr_u)^2)
F14C_EC100_0_charr_c_u<-sqrt((dFcharr_c_slope*charr_slope_u)^2 + (dFcharr_c_yield*yield_u)^2 + (F14C_EC100_0_charr_b_u^2) )
F14C_EC100_0_charr_u<-sqrt(F14C_EC100_0_charr_a_u^2 + F14C_EC100_0_charr_b_u^2 + F14C_EC100_0_charr_c_u^2)  

#linear slope calculation for extrapolation to 100% EC
linear_slope <- (F14C_EC100-F14C_EC_raw_data)/(1-df_stats_mean[,2])
colnames(linear_slope) <- c("linear_slope")

#export result as csv
df_stats_mean_F14C <- cbind(df_stats_mean[,7],df_stats_mean[,2:5],total_charr$total_charr_mean,as.data.frame(F14C_EC_raw_data),F14C_EC_u, F14C_EC100,F14C_EC100_u, linear_slope,F14C_EC100_0_charr,F14C_EC100_0_charr_u)
colnames(df_stats_mean_F14C) <- c("filter_name_short","EC_yield","charring_S1", "charring_S2","charring_S3", "charring_total","F14C_EC","F14C_EC_u","F14C_EC100","F14C_EC100_u", "linear_slope","F14C_EC100_0_charr","F14C_EC100_0_charr_u")
write.csv(df_stats_mean_F14C, file <- paste(basename(getwd()),"-mean-summary-with-F14C.csv",sep=""), row.names = F)

#plot raw vs corrected F14C values for each filter
F14C_raw <- as.data.frame(F14C_EC_raw_data)
F14C_raw$F14C_u <- F14C_EC_u
F14C_raw$filter_name_short <- df_stats_mean_F14C$filter_name_short
F14C_raw$class <- "uncorrected"
F14C_raw <- F14C_raw[, c("filter_name_short","F14C_EC_raw_data","F14C_u", "class")]
colnames(F14C_raw) <- c("filter_name_short","F14C","F14C_u", "class") 
#add values for EC-yield correction to 100%
F14C_EC100_corr <- df_stats_mean_F14C[,c(1,9)]
F14C_EC100_corr$F14C_u <- F14C_EC100_u
F14C_EC100_corr$class <- "corr. to 100% EC-yield"
colnames(F14C_EC100_corr) <- c("filter_name_short","F14C","F14C_u", "class") 

#add values for EC-yield correction to 100% and 0% charring
F14C_EC100_0_charr_corr <-  df_stats_mean_F14C[,c(1,12)]
F14C_EC100_0_charr_corr$F14C_u <- F14C_EC100_0_charr_u
F14C_EC100_0_charr_corr$class <- "corr. to 100% EC-yield at 0% charring"
colnames(F14C_EC100_0_charr_corr) <- c("filter_name_short","F14C","F14C_u", "class") 
F14C_EC100_0_charr_corr
all_F14C_data <- rbind(F14C_raw,F14C_EC100_corr,F14C_EC100_0_charr_corr)

#plot the corrected and uncorrected F14C value for each filter
theme_set(theme_classic(base_size = 13,base_family = "Helvetica"))
plot_all_F14C <- ggplot(all_F14C_data, aes(x=filter_name_short, y=F14C)) +
  geom_point(aes(shape= factor(class, levels = c("uncorrected","corr. to 100% EC-yield","corr. to 100% EC-yield at 0% charring")), color = factor(class, levels = c("uncorrected","corr. to 100% EC-yield","corr. to 100% EC-yield at 0% charring"))))+
  geom_errorbar(mapping=aes(x=filter_name_short, ymin=F14C-F14C_u, ymax=F14C+F14C_u), width=0.15, size=0.2, color="black") +
  xlab("Filter")+
  ylab(bquote(~ F^14~C))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,
                                                  if ( max(all_F14C_data$F14C) > 1) {
                                                    max(all_F14C_data$F14C)+0.3
                                                  } else {
                                                    1
                                                  }

                                                  ))
plot_all_F14C <- plot_all_F14C + theme( plot.margin = margin(1, 0.2, 0.2, 0.2, "cm"),legend.title = element_blank(), legend.position="top")
plot_all_F14C <- plot_all_F14C + guides(shape=guide_legend(nrow=3))

#summary figure with F14C and EC-yiled
fig_summary_F14C_top = ggarrange(
  plot_all_F14C,
  plot_EC_yield,
  ncol = 2,
  nrow = 1)
fig_summary_F14C_top

#create a summary figure with F14C, EC-yield, and charring. Export plots as pdf 
pdf(file =  paste(basename(getwd()),"-F14C_and_EC-yield-and-charring-summary.pdf",sep=""), width = 8.3 , height = 11.7)
fig_summary = ggarrange(
  fig_summary_F14C_top,
  fig_summary_charring,
  ncol = 1,
  nrow = 2)
annotate_figure(fig_summary, top = text_grob(paste("\n",basename(getwd())," summary",sep=""), color = "black" , face = "bold", size = 16),
               bottom = text_grob(paste("pdf generated by COMPYCALC", "   ", Sys.info()[["user"]],Sys.time(), "  ",sep=" "), color = "black",  face = "italic", size = 10),)
dev.off()

####################################################################################
#end COMPYCALC
####################################################################################

