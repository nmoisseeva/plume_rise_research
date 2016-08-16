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

vert_coord       = "ght_msl"
interp_levels    = (/0.  ,  0.02,  0.04,  0.06,  0.08,  0.1 ,  0.12,  0.14,  0.16, 0.18,  0.2 ,  0.22,  0.24,  0.26,  0.28,  0.3 ,  0.32,  0.34,0.36,  0.38,  0.4 ,  0.42,  0.44,  0.46,  0.48,  0.5 ,  0.52,0.54,  0.56,  0.58,  0.6 ,  0.62,  0.64,  0.66,  0.68,  0.7 ,0.72,  0.74,  0.76,  0.78,  0.8 ,  0.82,  0.84,  0.86,  0.88,0.9 ,  0.92,  0.94,  0.96,  0.98,  1.  ,  1.02,  1.04,  1.06,1.08,  1.1 ,  1.12,  1.14,  1.16,  1.18,  1.2 ,  1.22,  1.24,1.26,  1.28,  1.3 ,  1.32,  1.34,  1.36,  1.38,  1.4 ,  1.42,1.44,  1.46,  1.48/)
 
  opts             = True
  ;opts@extrapolate = False
  ;opts@field_type  = "t"

intrpT = wrf_user_vert_interp(fin,fld_T,vert_coord,interp_levels,opts)
printVarSummary(intrpT)
interpTKE = wrf_user_vert_interp(fin,fld_TKE,vert_coord,interp_levels,opts)
printVarSummary(interpTKE)

;write to netcdf
ncdfout = addfile(fout, "c")
ncdfout->T = intrpT
ncdfout->TKE = interpTKE
