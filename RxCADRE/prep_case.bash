#!/bin/bash

local=$pwd
filename="wrfout_L2G_cat3_new"

cd ~/wrf/wrf-fire/wrfv2_fire/test/em_fire/rxcadre

echo ".....copying wrf output to data folder"
cp wrfout_d01_* ~/data/plume/RxCADRE/$filename

echo ".....running regridding"
python2 ~/code/plume/RxCADRE/regrid_interp.py

cd $local
