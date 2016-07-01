#!/bin/sh

#run_sensitivity_bundle.bash
#=============================
#this script is a wrapper for performing a subset of sensitivity tests in wrf-fire

#input-----------------------------
testvar="wind"
range=(2v 4v 6v 8v 10v 12v)
input_dir="/Users/nadya2/code/plume/sensitivity/"
data_dir="/Users/nadya2/data/plume/"

#----------------------------------


echo "==========================================================="
echo "Running a bundle sensitivity simulation for case: $testvar"
echo "Testing a range of values: $range"
echo "==========================================================="

#local dir
local=${pwd}

#switch into wrf-fire directory
cd ~/Applications/WRFV3/test/em_fire/

#loop through the entire range in sensitivity analysis
for nRun in "${range[@]}"
do
    #clean up original settings files
    rm input_sounding
    rm namelist.input
    rm namelist.fire

    echo $nRun
    #symlink the proper setting files
    ln -s ${input_dir}${testvar}/namelist.fire ./
    ln -s ${input_dir}${testvar}/namelist.input ./
    ln -s ${input_dir}${testvar}/input\_sounding\_${nRun} ./input\_sounding
    echo "Running WRF for case: $nRun"
    ./ideal.exe
    echo "IDEAL.EXE complete"
    ./wrf.exe
    echo "WRF.EXE complete: $date"

    #rename and move the results to a data directory
    mkdir -p ${data_dir}${testvar}/ 
    mv wrfout\_d01\_0000\-01\-01\_00:00:00 ${data_dir}${testvar}/wrfout\_$nRun

done


#move to data directory
cd ${data_dir}${testvar}

for nRun in "${range[@]}"
do
    #create vdf files for VAPOR
    ncdfvdfcreate -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC wrfout\_$nRun $testvar\_$nRun.vdf
    ncdf2vdf -timedims Time -stagdims bottom_top_stag:south_north_stag:west_east_stag -vars U:V:W:T:QVAPOR:GRNHFX:AVG_FUEL_FRAC wrfout\_$nRun $testvar\_$nRun.vdf
done


cd ${local}

