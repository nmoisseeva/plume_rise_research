load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

filename = "/Users/nadya2/data/plume/comps/part1/wrfout_2"
fin      = addfile(filename+".nc","r")
fld = wrf_user_getvar(fin,"tk",-1)     ; The variable to interpolate

vert_coord       = "ght_msl"
interp_levels    = ispan(100,3000,100)
 
  opts             = True
  opts@extrapolate = True 
  opts@field_type  = "T"
  opts@logP        = False 

fld_intrp = wrf_user_vert_interp(fin,fld,vert_coord,interp_levels,opts)
printVarSummary(fld_intrp)
