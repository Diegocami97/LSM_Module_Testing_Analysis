#!/bin/bash

# Diego Venegas-Vargas
# LSM Module Testing 
# DAMIC-M Collaboration, JHU 
# This script processes images 3, 4, 5, 6 and 7 for low temperature data
# Usage: bash high_temp_image_analysis.sh

# Define ACM IDs 
export ACM_106="106"
export ACM_109="109"
export ACM_110="110" 

# # Define path for WADERS json files
# export path_json="/home/damicm/Soft/ccd-cdaq/LSM_Module_Testing_Scripts/Analysis"

# Define paths and parameters
# These parameters change for each dataset
export data_dir="Production_Modules/2025-08-21-DM04_DM05_DM06_Low_Temp_LBC_Parameters/" 
export Module_ACM_106="DM06"     
export Module_ACM_109="DM05"     
export Module_ACM_110="DM04" 


# Create directories to store output for each image and module
mkdir $data_dir/Analysis
mkdir $data_dir/Analysis/Image_3
mkdir $data_dir/Analysis/Image_4
mkdir $data_dir/Analysis/Image_5
mkdir $data_dir/Analysis/Image_6
mkdir $data_dir/Analysis/Image_7
mkdir $data_dir/Analysis/Image_3/$Module_ACM_106
mkdir $data_dir/Analysis/Image_3/$Module_ACM_109
mkdir $data_dir/Analysis/Image_3/$Module_ACM_110
mkdir $data_dir/Analysis/Image_4/$Module_ACM_106
mkdir $data_dir/Analysis/Image_4/$Module_ACM_109
mkdir $data_dir/Analysis/Image_4/$Module_ACM_110
mkdir $data_dir/Analysis/Image_5/$Module_ACM_106
mkdir $data_dir/Analysis/Image_5/$Module_ACM_109
mkdir $data_dir/Analysis/Image_5/$Module_ACM_110
mkdir $data_dir/Analysis/Image_6/$Module_ACM_106
mkdir $data_dir/Analysis/Image_6/$Module_ACM_109
mkdir $data_dir/Analysis/Image_6/$Module_ACM_110
mkdir $data_dir/Analysis/Image_7/$Module_ACM_106
mkdir $data_dir/Analysis/Image_7/$Module_ACM_109
mkdir $data_dir/Analysis/Image_7/$Module_ACM_110
mkdir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_106
mkdir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_109
mkdir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_110
mkdir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_106
mkdir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_109
mkdir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_110



# Process Image 4 for each module
echo "Processing Image 4 for Module $Module_ACM_106"
python src/Image_4.py --files $data_dir/avg_Image_4_Low_Temp_"$ACM_106"_*_*_*.fz --output_dir $data_dir/Analysis/Image_4/$Module_ACM_106 --module $Module_ACM_106

echo "Processing Image 4 for Module $Module_ACM_109"
python src/Image_4.py --files $data_dir/avg_Image_4_Low_Temp_"$ACM_109"_*_*_*.fz --output_dir $data_dir/Analysis/Image_4/$Module_ACM_109 --module $Module_ACM_109

echo "Processing Image 4 for Module $Module_ACM_110"
python src/Image_4.py --files $data_dir/avg_Image_4_Low_Temp_"$ACM_110"_*_*_*.fz --output_dir $data_dir/Analysis/Image_4/$Module_ACM_110 --module $Module_ACM_110

root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_4_'$Module_ACM_106'_list.txt",4,"'$Module_ACM_106'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_4_'$Module_ACM_109'_list.txt",4,"'$Module_ACM_109'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_4_'$Module_ACM_110'_list.txt",4,"'$Module_ACM_110'","'$data_dir'/Analysis")'

root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_5_'$Module_ACM_106'_list.txt",5,"'$Module_ACM_106'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_5_'$Module_ACM_109'_list.txt",5,"'$Module_ACM_109'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_5_'$Module_ACM_110'_list.txt",5,"'$Module_ACM_110'","'$data_dir'/Analysis")'

root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_6_'$Module_ACM_106'_list.txt",6,"'$Module_ACM_106'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_6_'$Module_ACM_109'_list.txt",6,"'$Module_ACM_109'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_6_'$Module_ACM_110'_list.txt",6,"'$Module_ACM_110'","'$data_dir'/Analysis")'

root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_7_'$Module_ACM_106'_list.txt",7,"'$Module_ACM_106'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_7_'$Module_ACM_109'_list.txt",7,"'$Module_ACM_109'","'$data_dir'/Analysis")'
root -l -q -b 'src/MakeComposite_Fe55_Clusters_Extensions.cc+("input_lists/Image_7_'$Module_ACM_110'_list.txt",7,"'$Module_ACM_110'","'$data_dir'/Analysis")'








#------------------------------------------------------------------------------------------------#
# # Uncomment the following lines if you have panaSKImg installed and want to process Image 3
# # Make sure to adjust the path to your json files accordingly
## Process Image 3 -- multi-skip Serial resgister image for single electron resolution
# # 500 skip
# panaSKImg "$data_dir/avg_Image_3_500_Low_Temp_SR_"$ACM_106"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_106 --acm --save-plots 
# panaSKImg "$data_dir/avg_Image_3_500_Low_Temp_SR_"$ACM_109"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_109 --acm --save-plots 
# panaSKImg "$data_dir/avg_Image_3_500_Low_Temp_SR_"$ACM_110"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_110 --acm --save-plots 

# 1000 skip
# panaSKImg "$data_dir/avg_Image_3_1000_Low_Temp_SR_"$ACM_106"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_106 --acm --save-plots 
# panaSKImg "$data_dir/avg_Image_3_1000_Low_Temp_SR_"$ACM_109"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_109 --acm --save-plots 
# panaSKImg "$data_dir/avg_Image_3_1000_Low_Temp_SR_"$ACM_110"_*_*_*.fz" -j $path_json/moduletest_lowtemp_image_3.json -o $data_dir/Analysis/Image_3/$Module_ACM_110 --acm --save-plots 
#------------------------------------------------------------------------------------------------#

# Process Image 3 for each module: at this point, the 500 and 1000 skips image files should be in the correct folders
# echo "Processing Image 3 for with 500 skips for Module $Module_ACM_106"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/500_skips/$Module_ACM_106 --output_dir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_106 --module $Module_ACM_106

# echo "Processing Image 3 for with 500 skips for Module $Module_ACM_109"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/500_skips/$Module_ACM_109 --output_dir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_109 --module $Module_ACM_109

# echo "Processing Image 3 for with 500 skips for Module $Module_ACM_110"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/500_skips/$Module_ACM_110 --output_dir $data_dir/Analysis/Image_3/500_skips/$Module_ACM_110 --module $Module_ACM_110

# echo "Processing Image 3 for with 1000 skips for Module $Module_ACM_106"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/1000_skips/$Module_ACM_106 --output_dir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_106 --module $Module_ACM_106

# echo "Processing Image 3 for with 1000 skips for Module $Module_ACM_109"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/1000_skips/$Module_ACM_109 --output_dir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_109 --module $Module_ACM_109

# echo "Processing Image 3 for with 1000 skips for Module $Module_ACM_110"
# python src/pdf_to_plot_dir.py --input_dir $data_dir/Image_3/1000_skips/$Module_ACM_110 --output_dir $data_dir/Analysis/Image_3/1000_skips/$Module_ACM_110 --module $Module_ACM_110

# End of script
