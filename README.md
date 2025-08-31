# LSM Module Testing Analysis
> Analysis scripts for LSM module testing 

database: https://ccdqc.pha.jhu.edu/QC_production/edit_module_underground.php

## Info tables

**Module layout:** copy from surface table

**Preliminary Grade Assessment:** leave to Alvaro or other reviewers

**Testing:** 

- Test Time: `ls -lt ~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules`
(use the timestamp of the first image)

- Temp: from SC at the first High/Low image 

- sequencer and bcf files: 

`scp -r "daq2:~/Soft/ldaq/config/configfiles/TC_LSM_REF_DM01.bcf" ~/Code/LSM_Module_Testing_Analysis/Data/configfiles`

`scp -r "daq2:~/Soft/ldaq/config/seqfiles/refrence_v4_DM01.seq" ~/Code/LSM_Module_Testing_Analysis/Data/seqfiles`


# Prepare: 


in ~/.zshrc (.bashrc), add:

`alias adc_plotter="python3 /Users/iumicuszhu/Code/LSM_Module_Testing_Analysis/src/adc_plotter.py"`
`alias pdf2plot="python3 /Users/iumicuszhu/Code/LSM_Module_Testing_Analysis/src/pdf_to_plot_dir.py"`

Start by creating a folder named Production_Modules_Data

This will host all of the data for high and low temperature

for fz files:

`scp -r "daq2:~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules/2025-08-05-DM01_DM02_DM03_High_Temp/avg*fz" ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

for trace/psd:

`scp -r "daq2:~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules/2025-08-05-DM01_DM02_DM03_High_Temp/*.csv" ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

<!-- in downloaded data folder (`DM*_DM*_DM*_High_Temp`): -->
<!-- `mkdir trace psd image1 image2 image31 image32 image4` -->

(run all the analysis scripts in this folder)

# analysis of images at high temperature 
----------------------------------------------
go into the folder: `cd ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

## 1. trace
- generate all the trace pdf images from the csv files:
`adc_plotter --trace --delta_t 10000 --fname trace*.csv`

- combine the 4 channels image together to one, and (optional) delete the single pdf files:
`pdf2plot --files 'trace_ch*_110_*.pdf' DM04_trace --cleanup`     

- get reference image names
`ls -l trace*110_*.csv`

## 2. psd
- generate png files from .csv:
`adc_plotter --psd --fname psd*.csv`

## 3. Image1_High_Temp
`python3 ../../src/Image_1.py --files avg_Image_1_High_Temp_106_*_.fz --module DM03`

## 4. Image2_High_Temp
`python3 ../../src/Image_2.py --files avg_Image_2_High_Temp_106_*_.fz --module DM03`

## 5. Image3_High_Temp

for now, use Wader in daq2 computer, then download generated pdfs to local

`setwadersACM`
`cd <datafolder>`
`mkdir temp`
`cd temp`

default:
- `panaSKImg "../avg_Image_3_500*.fz" -j ~/Soft/ccd-cdaq/LSM_Module_Testing_Scripts/Analysis/moduletest_hightemp_image_3_NotHighT.json -o . --acm --save-plots`

if bad fit (usually 1000 skips image):
- `panaSKImg "../avg_Image_3_1000*.fz" -j ~/Soft/ccd-cdaq/LSM_Module_Testing_Scripts/Analysis/moduletest_hightemp_image_3_IsHighT.json -o . --acm --save-plots`


download to local:

`scp -r "daq2:~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules/2025-08-25-DM04_DM05_DM06_High_Temp/temp/*FitDarkCurrentProcess_DCfit*" ~/Code/LSM_Module_Testing_Analysis/Data/DM04_DM05_DM06_High_Temp`

select the fit results to your liking (if similar, use the last set) per each module, each extension, delete the rest, and run:
`pdf2plot --files 'avg_Image_3_1000_*_SR_110_*.pdf' DM04_image32 --cleanup`


------------------------------
cd <Module_data_set>
mkdir Analysis
mkdir Image_1;mkdir Image_2;mkdir Image_3;mkdir Image_4;
> python scripts, these are used for images 1, 2, and 4 mainly for column defects. 

cd Image_1
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_1_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

@@ -47,12 +50,13 @@ python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Imag
cd Image_4
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

# analysis of images at low temperature
--------------------------------
cd <Module_data_set>
mkdir Analysis
mkdir Image_4;

> python scripts, these are used for image 4 mainly for column defects, CTI, noise and sharpness of tracks

cd Image_4
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>







