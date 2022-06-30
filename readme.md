[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5958275.svg)](https://doi.org/10.5281/zenodo.5958275) ![COMPYCALC-logo](compycalc-logo.png "COMPYCALC logo")

# COMprehensive Yield CALCulation

#### A tool for EC yield extrapolation and charring correction

## Describtion

COMPYCALC is a R script for EC yield extrapolation and charring correction. The script uses the the raw data output from a thermo-optical OC/EC analyzer (Model 5L, [Sunset Laboratory Inc.](https://www.sunlab.com), OR, United States) running the Swiss_3S protocol for OC/EC separation developed by [Zhang et al. (2012)](https://doi.org/10.5194/acp-12-10841-2012) for EC yield and charring calculation. Using F<sup>14</sup>C(EC) values measured by accelerator mass spectrometry (AMS) and calculated F<sup>14</sup>OC values, the script performs the EC yield extrapolation to 100% EC yield and a charring correction to 0% charring for F<sup>14</sup>C(EC) values.

## Usage

To use the run COMPYCALC program, follow the steps written in the comment section of the `compycalc.R` file. This is the file you want to run, the other files in the subfolder (zsrc) are linked to this script.

### Step 1: set up environment

In the first section, you are asked to set the working directory, either with the command `setwd()` or by going to Session → Set Working Directory → To Source File Location if you are using [R Studio](https://rstudio.com).
    
### Step 2: add OC/EC analyzer files

In the second step, you are ask to add your Sunset raw files in folders to the working directory folder. Obviously, you can also do it the other way around by adding the `compycalc.R` script to the folder where your data is. Please be aware that the script will take the last digit of each folder for naming, so make sure that you name your folders accordingly.

**OPTION**: add a EC, TC or Swiss_4S run raw data file recorded with the same Sunset OC/EC analyzer oven conditions than your samples. Add the file into the *zsrc* folder and rename it to `custom_cooldown.csv`. This is to correct for the opacity of the Sunset OC/EC analyzer to accurately calcualte the EC yield. This option is **highly recommended**. Otherwise, the default `custom_cooldown.csv` file will be used, which is a generic cooldown file reflecting a new analyzer oven. 

**Note:** delete all unnecessary files (including hidden files) in the folder you want to run COMPYCALC. Keep only the Sunset raw file folders as described above, the `compycalc.R` script, the *zsrc* folder containing additional scripts and the cooldown data. 

### Step 3: add radiocarbon data

Last but not least you have to add the F<sup>14</sup>C(EC) and F<sup>14</sup>C(OC) raw data with uncertainties as separate csv files. F<sup>14</sup>C(EC) contains the measured F<sup>14</sup>C(EC) values in the first column and measurement uncertainties in the second column. For OC you do the same: F<sup>14</sup>C(OC) contains the calculated F<sup>14</sup>C(OC) values in the first column and uncertainties in the second column. Note that the files need to be in sample order. The csv files must be in the working directory, i.e. the folder where your `compycalc.R` file is.

### Step 4: run code

Finally, you are ready to run the COMPYCALC script. 
    
As an output, you will get:       
   
*   Each folder with Sunset measurement files will get five output files:
	*   calc-summary-plots.pdf
	*   clean-results.csv
	*   mean-results.csv
	*   raw-results.csv
	*   stats.csv
*   The working directory folder will get two files:
	*   mean-summary-with-F14C.csv
	*   F14C-and-EC-yield-and-charring-summary.pdf
   

## How does COMPYCALC work?

COMPYCALC (COMprehensive Yield CALCulation) consists of three subscripts for data input and output, EC yield and charring, as well as an extrapolation of the F<sup>14</sup>C(EC) values to 100% EC yield. For each sample, the OC/EC analyzer raw data files containing the laser transmission signal for each OC removal run need to be in a designated subfolder. Additionally, the script requires the uncorrected F<sup>14</sup>C(EC) and F<sup>14</sup>C(OC) data in separate files (csv format) in the main folder. The data input and output script loads the OC/EC analyzer raw data files for each sample folder and initiates the calculation with the EC yield and charring script. The results written in each sample folder is then read by the main script and forwarded to the second calculation script for the extrapolation to 100% EC yield. Finally, the F<sup>14</sup>C(EC) value extrapolated to 100% EC yield is corrected for charring in the main script, as this should be regarded as an OC contamination of the measured EC. After all calculations, a summary data file (csv) with overall EC yield, the charring contribution for each OC removal step (S1, S2, S3), the total charring contribution as well as the raw F<sup>14</sup>C(EC), F<sup>14</sup>C(EC) extrapolated to 100% EC yield, and F<sup>14</sup>C(EC) extrapolated to 100% EC yield and corrected for charring is generated as an output. Additionally, a summary pdf is generated with plots for all F<sup>14</sup>C results, EC yields, and charring for each step (S1, S2, S3).

![COMPYCALC scheme](How-does-COMPYCALC-work.png "COMPYCALC scheme")

## Authors

This tool was written by [Martin Rauber](https://www.martin-rauber.com) and [Gary Salazar](mailto:gary.salazar@dcb.unibe.ch) for LARA, the [Laboratory for the Analysis of Radiocarbon with AMS](https://www.14c.unibe.ch) at the University of Bern. Please get in touch for any bug fixes and suggestions!

## Licence

COMPYCALC is released under the [MIT License](LICENCE.txt).
