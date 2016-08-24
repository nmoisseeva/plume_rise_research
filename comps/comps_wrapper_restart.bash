#!/bin/bash


#comps_wrapper.bash
#=============================
#this script is a wrapper for performing performing Rolands comps subsections
#assumes restart run

#input-----------------------------
part="D"
section=(a b c)
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
    #switch into wrf-fire directory

    cd /Users/nadya2/Applications/WRFV3/test/em\_fire

    #clean up original settings files
    rm input_sounding
    rm namelist.input
    rm namelist.fire
   rm wrfrst* 

    echo "..... subsection $nRun"
    #symlink the proper setting files
    ln -s ${input_dir}${part}/namelist.fire ./
    ln -s ${input_dir}${part}/namelist.input ./
    ln -s ${input_dir}${part}/input\_sounding\_${nRun} ./input\_sounding
    ln -s ${data_dir}B/restart/wrfrst_${nRun} ./wrfrst_d01_0000-01-01_01:00:00
    echo "Running WRF for restart case: $nRun"
    ./wrf.exe 
    echo ".....WRF.EXE complete: $date"

    #rename and move the results to a data directory
    mkdir -p ${data_dir}${part}/ 
    mv wrfout\_d01\_0000\-01\-01\_01:01:00 ${data_dir}${part}/wrfout\_$nRun

    # # #run interpolation code
    # echo "Running NCL for data interpolation: $nRun"
    # cd ${input_dir}${part}/
    # mkdir -p interp
    # ln -s  ${data_dir}${part}/wrfout\_$nRun ./interp/wrfout\_temp
    # ncl63 data_prep.ncl 
    # mv intperp\_wrf.nc ./intperp/intperp\_wrf\_${nRun}
    
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
