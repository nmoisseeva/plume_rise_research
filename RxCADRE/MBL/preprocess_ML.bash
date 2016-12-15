#!/bin/bash

#====================================
#preprocess_ML.bash
#	created: Dec 2016
#	by: nmoisseeva@eos.ubc.ca
#	-script convert to vdf
#====================================

#input-----------------------------
datadir="/Users/nadya2/data/plume/RxCADRE/MBL/"
wrfout="wrfout_MLbubble"
firename='ML'
#----------------------------------

echo "==========================================================="
echo "Running vdf converstion for $wrfout"
echo "==========================================================="

#local dir
local=${pwd}
cd $datadir

echo ".....converting data to vdf"
ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars T:TKE: $wrfout $firename.vdf
ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars T:TKE: $wrfout $firename.vdf

cd ${local}
echo "COMPLETE"
