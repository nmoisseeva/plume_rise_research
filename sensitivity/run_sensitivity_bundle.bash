#!/bin/sh

#run_sensitivity_bundle.bash
#=============================
#this script is a wrapper for performing a subset of sensitivity tests in wrf-fire

#input-----------------------------
testvar="wind"
range=("2v","4v","6v","8v","10v","12v")
input_dir="~/code/plume/sensitivity/"

#----------------------------------

#local dir
local=${pwd}

#switch into wrf-fire directory
cd ~/Applications/WRFV3/test/em_fire/

#clean up original settings files
rm input_sounding
rm namelist.input
rm namelist.fire

#loop through the entire range in sensitivity analysis
for nRun in "${range[@]}"
do
#symlink the proper setting files
ln -s ${input_dir}/${testvar}/namelist.fire ./
ln -s ${input_dir}/${testvar}/namelist.input ./
ln -s ${input_dir}/${testvar}/input\_sounding\_${nRun} ./input\_sounding
done


cd ${local}
