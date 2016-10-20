#!/bin/bash


#interp_set.bash
#=============================
#this script calls ncl to interpolate a set of data

#input-----------------------------
input_dir="/Users/nadya2/code/plume/RxCADRE/"
data_dir="/Users/nadya2/data/plume/RxCADRE/"
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
mkdir -p ${data_dir}/interp 		#dir for saving interp data
rm ./interp_temp/wrfout 			#cleanup
ln -s ${data_dir}/wrfout_${firename} ./interp_temp/wrfout  	#link the worout file temporarily
$NCARG_ROOT/ncarg6.3/bin/ncl interp.ncl 					#run ncl code on it
mv wrfinterp.nc ${data_dir}/interp/wrfinterp_${firename}    #move interpolated data to data dir
rm -R ./interp_temp 				#cleanup




# prep for VAPOR
# #move to data directory
# cd ${data_dir}${part}

# for nRun in "${section[@]}"
# do
#     #create vdf files for VAPOR
#     ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC wrfout\_$nRun $testvar\_$nRun.vdf
#     ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC wrfout\_$nRun $testvar\_$nRun.vdf
# done


cd ${local}
