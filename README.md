# LSM Module Testing Analysis
> Analysis scripts meant for LSM module testing 

# copy the data folder to local

Start by creating a folder named Production_Modules_Data

This will host all of the data for high and low temperature

`scp -r "daq2:~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules/2025-08-05-DM01_DM02_DM03_High_Temp/avg*fz" ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

`scp -r "daq2:~/Soft/ccd-cdaq/LSM_Module_Testing_Data/Production_Modules/2025-08-05-DM01_DM02_DM03_High_Temp/*.csv" ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

# analysis of images at high temperature 
----------------------------------------------
go into the folder: `cd ~/Code/LSM_Module_Testing_Analysis/Data/DM01_DM02_DM03_High_Temp`

## 1. trace

- generate all the trace png images from the csv files:
`python3 ../../src/adc_plotter.py --trace --delta_t 10000 --fname trace*.csv`

- combine the 4 channels image together to one:
`python3 ../../src/pngs_to_plot.py "trace_ch*_109*.png" --outname DM02`     

- reference image
`ls -l trace*110_*.csv`

## 2. psd
- generate png files from .csv:
`python3 ../../src/adc_plotter.py --psd --fname psd*.csv`

## 3. Image1_High_Temp
`python3 ../../src/Image_1.py --files avg_Image_1_High_Temp_106_*_.fz --module DM03`

## 4. Image2_High_Temp
`python3 ../../src/Image_2.py --files avg_Image_2_High_Temp_106_*_.fz --module DM03`

cd <Module_data_set>
mkdir Analysis
mkdir Image_1;mkdir Image_2;mkdir Image_3;mkdir Image_4;
> python scripts, these are used for images 1, 2, and 4 mainly for column defects. 

cd Image_1
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_1_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

cd Image_2
python src/Image_1.py --files Production_Modules_Data/<Module_data_set>/avg_Image_2_High_Temp_<ACM#>_*_*_*.fz --module <ModuleID>

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







