###################################################################################
#Version information:
###################################################################################


####################################################################################
# prep.
####################################################################################

#Set working directory 
#setwd("")

#clean up environment---------------------------------------------------------------
rm(list=ls())
if(!is.null(dev.list())) dev.off()

#load libraries & install packages if necessary-------------------------------------
load.lib = c("ggpubr", "dplyr", "data.table", "purrr", "stringr", "pastecs","readxl")   
install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

####################################################################################
#script
####################################################################################

####################################################################################
#Method 1: allRyield calculation summary (import Sunset txt, calc. for each individual filter , summary)
####################################################################################

#load csv files--------------------------------------------------------------------
#raw results from each filter
df = list.files(".",pattern = "*raw-results.*csv") %>% map_df(~fread(.))
colnames(df) = c("filtercounter","EC_yield", "charring_S1", "charring_S2", "charring_S3", "filter_name")
df$filter_name_short=str_sub(df$filter_name,(nchar(df$filter_name)),nchar(df$filter_name))

df_charring = df[,3:7]
df_charring

#stats file
df_stats = list.files(".",pattern = "*-stats.*csv") %>% map_df(~fread(.))
df_stats$filter_name_short = c(rep(unique(df$filter_name_short),each = 14))
df_stats
# number of filters for each sample
df_stats[is.element(df_stats$V1, "nbr.val"),]
# mean for each  sample and export to csv
df_stats_mean=df_stats[is.element(df_stats$V1, "mean"),]
write.csv(df_stats_mean, "mr01-143-mean-summary.csv", row.names = T)

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
  ylab("EC yield")
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
  align = "h" ,ncol = 1, 
  nrow = 3)
fig_summary_charring = annotate_figure(fig_summary_charring,left = text_grob("charring", color = "black", size = 14, rot = 90))
fig_summary_charring
