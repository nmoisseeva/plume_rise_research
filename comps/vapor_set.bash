#!/bin/bash


#vapor_set.bash
#=============================
#this script converts data to VDF format for vapor

#input-----------------------------
part="B"
section=(a)
#section=(a)
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
   cd $data_dir/$part/interp/
   # ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
   # ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
   # ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfinterp\_$nRun vapor_$nRun.vdf
   # ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfinterp\_$nRun vapor_$nRun.vdf

   ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:W:HGT wrfinterp\_$nRun vapor_$nRun.vdf
   ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:W:HGT wrfinterp\_$nRun vapor_$nRun.vdf

done

cd ${local}
