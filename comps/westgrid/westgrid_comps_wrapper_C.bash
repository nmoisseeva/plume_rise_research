#!/bin/bash


#comps_wrapper.bash
#=============================
#this script is a wrapper for performing performing Rolands comps subsections

#input-----------------------------
part="C"
section=(a b c)
input_dir="/home/moisseev/plume/comps/"
data_dir="/home/moisseev/data/plume/comps/"
wrf_dir="/home/moisseev/WRFV3C/"

#----------------------------------

echo "==========================================================="
echo "Running comps section bundle: $part"
echo "==========================================================="

#local dir
local=${pwd}
source ~/.profile


#loop through all subsections
for nRun in "${section[@]}"
do
    #switch into wrf-fire directory
    cd ${wrf_dir}

    #clean up original settings files
    rm ${wrf_dir}input_sounding
    rm ${wrf_dir}namelist.input
    rm ${wrf_dir}namelist.fire
    rm ${wrf_dir}wrfrst\_d01\_0000\-01\-01\_01:00:00

    echo "..... subsection $nRun"
    #symlink the proper setting files
    ln -s ${input_dir}${part}/namelist.fire ${wrf_dir}
    ln -s ${input_dir}${part}/namelist.input ${wrf_dir}
    ln -s ${input_dir}${part}/input\_sounding\_${nRun} ${wrf_dir}input\_sounding
    ln -s ${data_dir}restart/wrfrst_${nRun} wrfrst\_d01\_0000\-01\-01\_01:00:00
    echo "Running WRF for case: $nRun"
    #${wrf_dir}ideal.exe
    #echo ".....IDEAL.EXE complete"
    ${wrf_dir}wrf.exe 
    echo ".....WRF.EXE complete: $date"

    #rename and move the results to a data directory
    mkdir -p ${data_dir}${part}/ 
    mv ${wrf_dir}wrfout\_d01\_0000\-01\-01\_01:01:00 ${data_dir}${part}/wrfout\_$nRun

    ##run interpolation code
    #echo "Running NCL for data interpolation: $nRun"
    #cd ${input_dir}${part}/
    #mkdir -p interp
    #ln -s  ${data_dir}${part}/wrfout\_$nRun ./interp/wrfout\_temp
    #ncl63 data_prep.ncl 
    #mv `intperp\_wrf.nc ./intperp/intperp\_wrf\_${nRun}
    
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