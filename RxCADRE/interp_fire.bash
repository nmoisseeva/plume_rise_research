#!/bin/bash


#interp_set.bash
#=============================
#this script calls ncl to interpolate a set of data

#input-----------------------------
# input_dir="/Users/nadya2/code/plume/RxCADRE"
data_dir="/Users/nadya2/data/plume/RxCADRE"
firename="LG2"

#----------------------------------

echo "==========================================================="
echo "Running RxCADRE interpolation for $firename"
echo "==========================================================="

#local dir
local=${pwd}

#run interpolation code
echo "Running NCL for data interpolation: $firename"

mkdir -p interp_temp				#make temporary directory

rm ./interp_temp/wrfout 			#cleanup
ln -s ${data_dir}/wrfout_${firename} ./interp_temp/wrfout  	#link the wrfout file temporarily
echo ${pwd}
$NCARG_ROOT/ncarg6.3/bin/ncl interp.ncl 					#run ncl code on it


mkdir -p ${data_dir}/interp 		#dir for saving interp data
mv wrfinterp.nc ${data_dir}/interp/wrfinterp_${firename}    #move interpolated data to data dir
rm -R ./interp_temp 				#cleanup

cd ${local}
