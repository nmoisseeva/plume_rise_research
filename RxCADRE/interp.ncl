load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;interpolate data in the vertical
filename = "/Users/nadya2/code/plume/RxCADRE/interp_temp/wrfout"
savedir = "/Users/nadya2/data/plume/RxCADRE/interp_temp/"

fout = "wrfinterp.nc"
system("rm -f wrfinterp.nc")
fin      = addfile(filename,"r")
fld_XTIME = wrf_user_getvar(fin,"XTIME",-1)     ; The variable to interpolate
fld_QVAPOR = wrf_user_getvar(fin,"QVAPOR",-1)
;fld_GRNHFX = wrf_user_getvar(fin,"GRNHFX",-1)
;fld_AVG_FUEL_FRAC = wrf_user_getvar(fin,"AVG_FUEL_FRAC",-1)

vert_coord       = "ght_msl"
interp_levels    = (/0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 , 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8 , 1.85, 1.9, 1.95/)
 
  opts             = True
  ;opts@extrapolate = False
  ;opts@field_type  = "t"

interpQVAPOR = wrf_user_vert_interp(fin,fld_QVAPOR,vert_coord,interp_levels,opts)
printVarSummary(interpQVAPOR)

;write to netcdf
setfileoption("nc","Format","LargeFile")
ncdfout = addfile(fout, "c")
ncdfout->QVAPOR = interpQVAPOR

;scale height

;copy as-is variables
;ncdfout->GRNHFX = fld_GRNHFX
;ncdfout->AVG_FUEL_FRAC = fld_AVG_FUEL_FRAC
ncdfout->XTIME = fld_XTIME
