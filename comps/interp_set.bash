#!/bin/bash


#interp_set.bash
#=============================
#this script calls ncl to interpolate a set of data

#input-----------------------------
part="D"
# section=(a b c)
section=(a)
input_dir="/Users/nadya2/code/plume/comps/"
data_dir="/Users/nadya2/data/plume/comps/"

#----------------------------------

echo "==========================================================="
echo "Running comps section bundle: $part"
echo "==========================================================="

#local dir
local=${pwd}

#loop through all subsections
for nRun in "${section[@]}"
do
 
   #run interpolation code
   echo "Running NCL for data interpolation: $nRun"

   mkdir -p interp_temp
   mkdir -p ${data_dir}${part}/interp
   rm ./interp_temp/wrfout
   ln -s ${data_dir}${part}/wrfout_${nRun} ./interp_temp/wrfout
   # $NCARG_ROOT/ncarg6.3/bin/ncl interp.ncl
   $NCARG_ROOT/ncarg6.3/bin/ncl interp_vdf.ncl
   mv wrfinterp.nc ${data_dir}${part}/interp/wrfinterp_${nRun}   
   rm -R ./interp_temp
done



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
