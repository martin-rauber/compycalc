# COMPYCALC: COMprehensive Yield CALCulation 
#### a tool to perform 14C OC/EC corrections

##  General usage 

To use the compycalc program, follow the four steps written in the comment section of the `compycalc.R` file. This is the file you want to run, the other files in the subfolder are linked to this script. In the first section, you are asked to set the working directory, either with the command `setwd()` or by going to Session --> Set Working Directory --> To Source File Location in R Studio. 

```
#set wd-----------------------------------------------------------------------------

# 1)  set the working directory (wd) for COMPYCALC: this folder must contain the compycalc R
#     file and the folder /zsrc containing the scripts. The wd name will be used to name the 
#     result files. 

#setwd("")
```

In the second step, you are ask to add your Sunset raw files in folders to the working directory folder. Obviously, you can also do it the other way around by adding the compycalc script to the folder where your data is. Now please be aware that the script will take the last digit of each folder for naming, so at best name your folders with xxx-[letter] (xxx-A, xxx-B, …). 

```
#add data---------------------------------------------------------------------------

# 2)  Sunset data
#     Add your folder(s) with the individual Sunset measurement(s) to the wd folder. 

```

Last but not least you have to add the F14C EC and OC raw data. These need to be in sample order and in a csv file in the working directory, i.e. the folder where your `compycalc.R` file is. 

```
# 3)  enter F14C EC and OC raw data
          ##EC: import from csv
          F14C_EC_raw_data = read.csv(list.files(".",pattern = "*EC-F14C-raw-data.csv" ,  recursive = TRUE), header = TRUE)
          F14C_EC_raw_data = F14C_EC_raw_data[,]
          #OR enter here manually
          #F14C_EC_raw_data = c()
          
          ##OC: import from csv
          F14C_OC_raw_data = read.csv(list.files(".",pattern = "*OC-F14C-raw-data.csv" ,  recursive = TRUE), header = TRUE)
          F14C_OC_raw_data = F14C_OC_raw_data[,]
          #OR enter here manually
          #F14C_OC_raw_data = c()
```

Finally, you are ready to run the code. 

```
# 4)  run script
          
#OUTPUT          
#         - each folder with measurements will get three files: "last-digit-of-folder"-mean-results.csv, "last-digit-of-folder"-raw-results.csv, "last-digit-of-folder"-stats.csv
#         - wd folder will get "your-wd-name-here"-mean-summary-with-F14C.csv and "your-wd-name-here"-F14C_and_EC-yield-and-charring-summary.pdf
```

These For best results, you use the uncorrected EC and OC values, run compycalc and and use the corrected EC and OC F14C values for a second run. After this iteration, the differences should be minuscule and you can use the fine EC and OC data. 

## Other settings

There are a few things you can also change in the first section of the `yields_calc_ext.R` file, however, preferably do not change anything if you do not know what you are doing. Also, if you make any changes to the program locally, make sure that you note that. Because you are working with a script and not with a package or a shiny app, you need to be very careful about that. 

## Repository

The source code is available on [Github](https://github.com/martin-rauber/compycalc).

## About

This tool was written by [Martin Rauber](https://martin-rauber.com) and [Gary Salazar](mailto:gary.salazar@dcb.unibe.ch) for LARA, the Laboratory for the Analysis of Radiocarbon with AMS at the University of Bern. Please get in touch for any bug fixes and suggestions!
