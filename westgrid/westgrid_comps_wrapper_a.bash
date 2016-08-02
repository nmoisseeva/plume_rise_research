#!/bin/bash


#comps_wrapper.bash
#=============================
#this script is a wrapper for performing performing Rolands comps subsections

#input-----------------------------
part="B"
section=(a)
input_dir="/home/moisseev/plume/comps/"
data_dir="/home/moisseev/plume/comps/westgrid/"

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
    cd $WRFRUNDIR

    #clean up original settings files
    rm ${WRFRUNDIR}input_sounding
    rm ${WRFRUNDIR}namelist.input
    rm ${WRFRUNDIR}namelist.fire

    echo "..... subsection $nRun"
    #symlink the proper setting files
    ln -s ${input_dir}${part}/namelist.fire $WRFRUNDIR
    ln -s ${input_dir}${part}/namelist.input $WRFRUNDIR
    ln -s ${input_dir}${part}/input\_sounding\_${nRun} ${WRFRUNDIR}input\_sounding
    echo "Running WRF for case: $nRun"
    ${WRFRUNDIR}ideal.exe
    echo ".....IDEAL.EXE complete"
    ${WRFRUNDIR}wrf.exe 
    echo ".....WRF.EXE complete: $date"

    #rename and move the results to a data directory
    mkdir -p ${data_dir}${part}/ 
    mv ${WRFRUNDIR}wrfout\_d01\_0000\-01\-01\_00:00:00 ${data_dir}${part}/wrfout\_$nRun

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
