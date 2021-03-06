#!/bin/bash

#====================================
#preprocess_case.bash
#	created: Oct 2016
#	by: nmoisseeva@eos.ubc.ca
#	-script to copy data from WRF directory and convert to vdf
#====================================

#input-----------------------------
# wrfdir="/Users/nadya2/Applications/WRFV3/test/em_fire/"
wrfdir="/Users/nadya2/Applications/WRF-SFIRE/wrf-fire/WRFV3/test/em_fire/rxcadre/"
datadir="/Users/nadya2/data/plume/RxCADRE/regrid/"
wrfout="wrfout_d01_2012-11-10_12:27:00"
firename="LG2"
#----------------------------------

echo "==========================================================="
echo "Running vdf converstion for $wrfout: $firename"
echo "==========================================================="

#local dir
local=${pwd}
newname=wrfout_${firename}_regrid

# #run interpolation code
cd $datadir
# echo ".....copying data from wrf directory"
# cp ${wrfdir}${wrfout} ${datadir}${newname}

echo ".....converting data to vdf"
# ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
# ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC:HGT wrfout\_$nRun vapor_$nRun.vdf
ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX: $newname $firename.vdf
ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars QVAPOR:GRNHFX: $newname $firename.vdf

cd ${local}
echo "COMPLETE"
