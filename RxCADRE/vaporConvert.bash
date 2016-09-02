#!/bin/bash


#vapor_set.bash
#=============================
#this script converts data to VDF format for vapor

#input-----------------------------
input_dir="/Users/nadya2/code/plume/RxCADRE/"
data_dir="/Users/nadya2/data/plume/RxCADRE/"
wrfout="wrfout_d01_2012-11-10_12:25:00"
firename="LG2"

#----------------------------------

echo "==========================================================="
echo "Running vdf converstion for $wrfout"
echo "==========================================================="

#local dir
local=${pwd}

#run interpolation code
cd $data_dir
# ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
# ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX:HGT $wrfout firename.vdf
ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX:HGT $wrfout firename.vdf

cd ${local}
