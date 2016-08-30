load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/ncarg6.3/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

;interpolate data in the vertical
filename = "/Users/nadya2/code/plume/comps/interp_temp/wrfout"
savedir = "/Users/nadya2/data/plume/comps/interp_temp/"

fout = "wrfinterp.nc"
system("rm -f wrfinterp.nc")
fin      = addfile(filename,"r")
fld_T = wrf_user_getvar(fin,"T",-1)     ; The variable to interpolate
fld_TKE = wrf_user_getvar(fin,"TKE",-1)
fld_QVAPOR = wrf_user_getvar(fin,"QVAPOR",-1)
fld_HGT = wrf_user_getvar(fin,"HGT",-1)
fld_GRNHFX = wrf_user_getvar(fin,"GRNHFX",-1)
fld_AVG_FUEL_FRAC = wrf_user_getvar(fin,"AVG_FUEL_FRAC",-1)

vert_coord       = "ght_msl"
interp_levels    = (/0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3, 0.33, 0.36, 0.39, 0.42, 0.45, 0.48, 0.51, 0.54, 0.57, 0.6 , 0.63, 0.66, 0.69, 0.72, 0.75, 0.78, 0.81, 0.84, 0.87, 0.9 , 0.93, 0.96, 0.99, 1.02, 1.05, 1.08, 1.11, 1.14, 1.17, 1.2 , 1.23, 1.26, 1.29, 1.32, 1.35, 1.38, 1.41, 1.44, 1.47, 1.5 , 1.53, 1.56, 1.59, 1.62, 1.65, 1.68, 1.71, 1.74, 1.77, 1.8 , 1.83, 1.86, 1.89, 1.92, 1.95, 1.98/)
 
  opts             = True
  ;opts@extrapolate = False
  ;opts@field_type  = "t"

intrpT = wrf_user_vert_interp(fin,fld_T,vert_coord,interp_levels,opts)
printVarSummary(intrpT)
interpTKE = wrf_user_vert_interp(fin,fld_TKE,vert_coord,interp_levels,opts)
printVarSummary(interpTKE)
interpQVAPOR = wrf_user_vert_interp(fin,fld_QVAPOR,vert_coord,interp_levels,opts)
printVarSummary(interpQVAPOR)

;write to netcdf
setfileoption("nc","Format","LargeFile")
ncdfout = addfile(fout, "c")
ncdfout->T = intrpT
ncdfout->TKE = interpTKE
ncdfout->QVAPOR = interpQVAPOR

;scale height
fld_HGT(:,:,:) = fld_HGT(:,:,:)/2000.

;copy as-is variables
ncdfout->HGT = fld_HGT
ncdfout->GRNHFX = fld_GRNHFX
ncdfout->AVG_FUEL_FRAC = fld_AVG_FUEL_FRAC
