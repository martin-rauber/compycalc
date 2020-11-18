###################################################################################
#Version information: v1.0.3
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

#Set working directory 
#setwd()
#setwd("~/OneDrive/Uni Bern/EC-yield_and_charring_calculation_R/R skript development/Version 1.0.3")

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
load.lib = c("ggpubr", "pastecs", "dplyr", "ggplot2", "pracma", "signal")   
install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

#clean up environment---------------------------------------------------------------
rm(list=setdiff(ls(), c("result_filename", "csv_raw", "csv_stat", "csv_mean","fitting_type","manual.coef", "r_scripts", "home_wd")))
if(!is.null(dev.list())) dev.off()

#load function----------------------------------------------------------------------

data_load_func = function(filename) {
  cooldown = read.csv("../cooldown_data.csv", sep = ",", header = T)
  dat = as.data.frame(read.csv(file = filename, sep = ",", skip = 28, header = T ))[,c(1:18)]
  tabla_complete <<- rbind(dat, cooldown)
  yield_calc = function(tabla_complete, fitting_type, manual.coef) {source("../yields_calc_ext.R")} 
  yield_calc(tabla_complete)
}

#load data, run calculation --------------------------------------------------------
filename = dir(".",pattern="^(.*)txt$")
df = NULL
for (i in filename){
  data_load_func(i)
  df = rbind(df, data.frame(tabla_resultados2$EC_yield,tabla_resultados2$charringS1,tabla_resultados2$charringS2,tabla_resultados2$charringS3))
}
if(!is.null(dev.list())) dev.off()
rm(list=setdiff(ls(), c("df","result_filename", "csv_raw", "csv_stat", "csv_mean", "r_scripts", "home_wd")))
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
theme_set(theme_classic(base_size = 13))
plot_yield = ggboxplot(data=df_yield, x = "dummy", y = "EC.yield",  xlab="", ylab = "EC-yield")
plot_yield

#plot: charring for each Sunset step
plot_charr = ggboxplot(data=df_charr, x="Sunset.step", y="charring.value", xlab="Sunset step", ylab = "charring")+
  geom_hline(yintercept = 0, color="red", lty= 5)
plot_charr 

fig = ggarrange(plot_yield, plot_charr, labels = c("A", "B"),ncol = 2, nrow = 1)
fig

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


####################################################################################
#changelog
#V1.0.2   corrections from Gary 02.04.2020
#V1.0.1   automatic script to load libraries and install package if necessary 
#V1.0.0   code based on dev. version mr20200330a


####################################################################################

