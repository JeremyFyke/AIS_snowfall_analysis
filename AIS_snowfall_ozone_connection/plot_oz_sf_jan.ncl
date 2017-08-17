begin
jan = False
;if (jan.eq.True)
;load "../common_ncl_utilies/composite_utilities.ncl"
;load "../common_ncl_utilies/general_utilities.ncl"
;load "../common_ncl_utilies/plotting_utilities.ncl"
;else
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"
;end if

LoadFreshData=True


if (jan.eq.True)
NO3SD = "../../../data/fixed-ozone/"
O3SD = "../../../data/lens/"
else
NO3SD="/glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/"
O3SD="/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/"
end if


;Load CAM grid data
if (jan.eq.True)
in:=addfile("../input_data/CAM_base_data.nc","r")
else
in:=addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/CAM_base_data.nc","r")
end if
lev=in->lev
lat=in->lat
lon=in->lon 

nYearClimatology=40
nMonthClimatology=nYearClimatology*12
nYrO3=86
iO3ClimMons=(nYrO3-nYearClimatology)*12
iO3ClimMone=(nYrO3)*12-1
nYrNO3=51
iNO3ClimMons=(nYrNO3-nYearClimatology)*12
iNO3ClimMone=(nYrNO3)*12-1
nLENS=42
NoOzoneNum=(/1,10,12,15,18,20,22,28/)
nNoOzoneEns=dimsizes(NoOzoneNum)

year2sec=31557600.

fname="../input_data/Ozone_analysis_jan.nc"

if (LoadFreshData) then ;Only do this ensemble-wide data load if data processing changes.  Otherwise, for plotting, just load data from file.
   ;get_LandArea/get_LandFrac/get_AIS_mask functions: in general_utilities.ncl
   LandArea = get_LandArea() 
   LandFrac = get_LandFrac()
   AISMask = get_AIS_mask()
   
   AISMaskBinary=where(AISMask.gt.0,1,0)
   
   print("Loading LENS runs.")
   PRECIP_O3=new((/nYrO3*12,192,288,nLENS/),"double")
   MonthlyO3TS=new((/nYrO3*12,nLENS/),"double")

   PRECIP_O3(:,:,:,:)=0.0
   MonthlyO3TS(:,:)=0.0 
   do e=1,42 ; skip first ensemble member for now because it has different file start/end date conventions.	  
      ie=e-1
      if (e.gt.35)
      elong=sprinti("%0.3i", e+65)
      else
      elong=sprinti("%0.3i", e)
      end if
      print(sprinti("%0.3i", e))
      ;LoadMonthlyData function: in composite_utilities.ncl
	  if (e.eq.1)
      PRECIP_O3(:,:,:,ie)=LoadMonthlyData_LENS1("PRECC",str_concat((/"ls ",O3SD,"PRECC/b.e11.B20TRC5CNBDRD.f09_g16.",elong,".cam.h0.PRECC.185001-200512.nc"/)))\
                	 +LoadMonthlyData_LENS1("PRECL",str_concat((/"ls ",O3SD,"PRECL/b.e11.B20TRC5CNBDRD.f09_g16.",elong,".cam.h0.PRECL.185001-200512.nc"/)))		     
	  else
	  PRECIP_O3(:,:,:,ie)=LoadMonthlyData("PRECC",str_concat((/"ls ",O3SD,"PRECC/b.e11.B20TRC5CNBDRD.f09_g16.",elong,".cam.h0.PRECC.192001-200512.nc"/)))\
	                	 +LoadMonthlyData("PRECL",str_concat((/"ls ",O3SD,"PRECL/b.e11.B20TRC5CNBDRD.f09_g16.",elong,".cam.h0.PRECL.192001-200512.nc"/)))		     
	  end if
					 
   end do
   PRECIP_O3:=PRECIP_O3*year2sec; convert to m/yr
   do e=1,42
      ie=e-1
      ;calculate_integrated_timeseries function: in general_utilities.ncl
	  MonthlyO3TS(:,ie)=calculate_integrated_timeseries(PRECIP_O3(:,:,:,ie),LandFrac,LandArea,AISMaskBinary)
   end do
   MonthlyO3TS=MonthlyO3TS/ 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr
   
   ;Make annual-average ensemble-average precip fields
   PRECIP_O3_annual_ensemble_mean=dim_avg_n(PRECIP_O3,0);time
   PRECIP_O3_annual_ensemble_mean:=dim_avg_n(PRECIP_O3_annual_ensemble_mean(:,:,:),2);ensemble
   ;add_lev_lat_lon_arrays function: in general_utilities.ncl 
   PRECIP_O3_annual_ensemble_mean:=add_lev_lat_lon_arrays(PRECIP_O3_annual_ensemble_mean,lev,lat,lon)

   print("Loading no-ozone-change runs.")
   PRECIP_NO3=new((/nYrNO3*12,192,288,nNoOzoneEns/),"double")
   MonthlyNO3TS=new((/nYrNO3*12,nNoOzoneEns/),"double")
   PRECIP_NO3(:,:,:,:)=0.0
   MonthlyNO3TS(:,:)=0.0
   ne=0
   do e=1,8
      ie=e-1
      elong=sprinti("%0.3i", NoOzoneNum(e-1))
      print(sprinti("%0.3i", NoOzoneNum(e-1)))
      PRECIP_NO3(:,:,:,ie)=LoadMonthlyData("PRECC",str_concat((/"ls ",NO3SD,"PRECC/b.e11.B20LE_fixedO3_",elong,".cam.h0.PRECC.195501-200512.nc"/)))\
                	  +LoadMonthlyData("PRECL",str_concat((/"ls ",NO3SD,"PRECL/b.e11.B20LE_fixedO3_",elong,".cam.h0.PRECL.195501-200512.nc"/)))		       
   end do   
   PRECIP_NO3:=PRECIP_NO3 * year2sec; convert to m/yr
   do e=1,8
      ie=e-1
      MonthlyNO3TS(:,ie)=calculate_integrated_timeseries(PRECIP_NO3(:,:,:,ie),LandFrac,LandArea,AISMaskBinary)
   end do
   MonthlyNO3TS=MonthlyNO3TS/ 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr
   
   ;Make annual-average ensemble-average precip fields
   PRECIP_NO3_annual_ensemble_mean=dim_avg_n(PRECIP_NO3,0);time
   PRECIP_NO3_annual_ensemble_mean:=dim_avg_n(PRECIP_NO3_annual_ensemble_mean(:,:,:),2);ensemble
   PRECIP_NO3_annual_ensemble_mean:=add_lev_lat_lon_arrays(PRECIP_NO3_annual_ensemble_mean,lev,lat,lon)

   system("rm -f "+fname)
   ncdf=addfile(fname ,"c")
   fAtt=True
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt ) ; copy file attributes
   filedimdef(ncdf,"time",-1,True)
   ncdf->MonthlyO3TS =MonthlyO3TS
   ncdf->MonthlyNO3TS=MonthlyNO3TS
   ;ncdf->PRECIP_O3_annual_ensemble_mean=PRECIP_O3_annual_ensemble_mean
   ;ncdf->PRECIP_NO3_annual_ensemble_mean=PRECIP_NO3_annual_ensemble_mean
else
   ncdf=addfile(fname,"r")
   MonthlyO3TS =ncdf->MonthlyO3TS
   MonthlyNO3TS=ncdf->MonthlyNO3TS
  ; PRECIP_O3_annual_ensemble_mean=ncdf->PRECIP_O3_annual_ensemble_mean
  ; PRECIP_NO3_annual_ensemble_mean=ncdf->PRECIP_NO3_annual_ensemble_mean   
end if
exit
;;;Analyze/plot stuff

;Calculate ensemble-mean stuff.  Store in arrays that are amenable for multi-line plotting.
Years=fspan(1,nYearClimatology,nYearClimatology)+2005.-nYearClimatology
monthsinDecimalYears=fspan(1,nMonthClimatology,nMonthClimatology)/12.+2005.-nYearClimatology

MonthlyEnsembleMeanStats= new((/6,nMonthClimatology/),double)
MonthlyEnsembleMeanStats(0,:)=dim_min_n(MonthlyO3TS,1)
MonthlyEnsembleMeanStats(1,:)=dim_max_n(MonthlyO3TS,1)
MonthlyEnsembleMeanStats(2,:)=dim_avg_n(MonthlyO3TS,1)

MonthlyEnsembleMeanStats(3,:)=dim_min_n(MonthlyNO3TS,1)
MonthlyEnsembleMeanStats(4,:)=dim_max_n(MonthlyNO3TS,1)
MonthlyEnsembleMeanStats(5,:)=dim_avg_n(MonthlyNO3TS,1)

AnnualO3TS=month_to_annual(MonthlyO3TS,1)
AnnualNO3TS=month_to_annual(MonthlyNO3TS,1)

AnnualEnsembleMeanStats=new((/9,nYearClimatology/),double)
AnnualEnsembleMeanStats(0,:)=dim_min_n(AnnualO3TS,1)
AnnualEnsembleMeanStats(1,:)=dim_max_n(AnnualO3TS,1)
AnnualEnsembleMeanStats(2,:)=dim_avg_n(AnnualO3TS,1)

AnnualEnsembleMeanStats(3,:)=dim_min_n(AnnualNO3TS,1)
AnnualEnsembleMeanStats(4,:)=dim_max_n(AnnualNO3TS,1)
AnnualEnsembleMeanStats(5,:)=dim_avg_n(AnnualNO3TS,1)

AnnualEnsembleMeanStats(6,:)=dim_avg_n(AnnualO3TS,1) - dim_avg_n(AnnualNO3TS,1) ;snowfall change attributable to ozone
AnnualEnsembleMeanStats(7,:)=runave(AnnualEnsembleMeanStats(6,:),10,0)

;Calculate statistical signficance of inter-ensemble differences, at each year.  To be improved on...!
do t=1,nYearClimatology
   tm1=t-1
   ave1=avg(AnnualO3TS(tm1,:))
   ave2=avg(AnnualNO3TS(tm1,:))
   var1=variance(AnnualO3TS(tm1,:))
   var2=variance(AnnualNO3TS(tm1,:))
   AnnualEnsembleMeanStats(8,tm1)=ttest(ave1,var1,ave2,var2,nLENS,nNoOzoneEns,False,False); statistical significance of inter-ensemble difference
end do

;Hard-code simple estimate of recent AIS mass loss (Gt/yr, Forsberg et al., 2017, 10.1007/s10712-016-9398-7) for plot comparisons
ForsbergObsStats=new((/3,4/),double)
ForsbergObsStats(0,:)=fspan(95.,95.,4) ;mean
ForsbergObsStats(1,:)=fspan(95.-50.,95.-50.,4) ;min
ForsbergObsStats(2,:)=fspan(95.+50.,95.+50.,4)
ForsbergTime=fspan(2002.,2005,4)

wks_type=get_default_wks(); function in general_utilities.ncl

;Plot monthy/annual raw ensemble-mean integrated AIS accumulation time series
   wks = gsn_open_wks (wks_type,"figs/PRECIP_O_TS");Plot timeseries and ensemble spread of integrated precip
   res = True                      ; plot mods desired
   res@trXMinF=min(Years)
   res@trXMaxF=max(Years)
   res@trYMinF=500.
   res@trYMaxF=4500.
   res@gsnDraw      = False  
   res@gsnFrame     = False  
   res@xyLineThicknesses = (/  1.0,  1.0, 10.0, 1.0,  1.0, 10.0/) 
   res@xyLineColors      = (/"blue","blue","blue","red","red","red"/) 
   res@xyDashPattern     := (/2,2,0,2,2,0/)
   monthly_plot = gsn_csm_xy(wks,monthsinDecimalYears,MonthlyEnsembleMeanStats,res)
   res@xyLineThicknesses :=(/20.0, 20.0/)
   res@xyLineColors :=(/"blue","red"/)
   res@xyDashPattern     := (/0,0/)
   data=new((/2,nYearClimatology/),double)
   data(0,:)=AnnualEnsembleMeanStats(2,:); grab the mean timeseries, condense into a new temporary array
   data(1,:)=AnnualEnsembleMeanStats(5,:)
   annual_plot = gsn_csm_xy(wks,Years-0.5,data,res)
   lineres=True
   lineres@lgLineThicknesses=(/20.0, 20.0/)
   lineres@lgLineColors=(/"blue","red"/)
   gres = True
   gres@YPosPosition = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   gres@XPosPosition = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
   textres = True
   textres@lgLabels = (/"Dynamic ozone ensemble","Fixed ozone ensemble"/) 
   annual_plot = simple_legend(wks,annual_plot,gres,lineres,textres)
   draw(annual_plot)
   overlay(monthly_plot,annual_plot)
   draw(monthly_plot)  
   frame(wks)

;plot difference in integrated AIS accumulation due to ozone, overlay against estimated loss trend
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_O_TS");Plot timeseries and ensemble spread of integrated precip
   res := True
   res@trXMinF=min(Years)
   res@trXMaxF=max(Years)
   res@trYMinF=-80.
   res@trYMaxF=180.
   res@gsnDraw      = False                        ; don't draw yet
   res@gsnFrame     = False                        ; don't advance frame yet
   res@xyLineThicknesses := (/20.0/)          ; make second line thicker
   res@xyLineColors      := (/"black"/)          ; change line color
   difference_plot = gsn_csm_xy(wks,Years-0.5,AnnualEnsembleMeanStats(6,:),res)

   res@xyLineColors      := (/"grey"/)          ; change line color
   running_average_plot = gsn_csm_xy(wks,Years-0.5,AnnualEnsembleMeanStats(7,:),res)

   res@xyLineThicknesses := (/20.0,  10.0, 10.0/)          ; make second line thicker
   res@xyDashPattern     := (/0,2,2/) 
   res@xyLineColors      := (/"orange","orange","orange"/)          ; change line color
   obs_plot = gsn_csm_xy(wks,ForsbergTime-0.5,ForsbergObsStats,res)
   lineres=True
   lineres@lgLineThicknesses=(/20.0, 20.0/)
   lineres@lgLineColors=(/"black","orange"/)
   gres = True
   gres@YPosPosition = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   gres@XPosPosition = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
   textres = True
   textres@lgLabels = (/"Ozone-attributable AIS accumulation change","2002-2015 AIS mass loss estimate"/) 
   difference_plot = simple_legend(wks,difference_plot,gres,lineres,textres)

   draw(annual_plot)
   overlay(difference_plot,running_average_plot)
   overlay(difference_plot,obs_plot)
   draw(difference_plot)
   frame(wks)


;Plot PRECIP_O3
   wks = gsn_open_wks (wks_type,"figs/PRECIP_O3")
   ;plot_field function: in plotting_utilities.ncl
;   plot=plot_field(PRECIP_O3_annual_ensemble_mean,0.,1.,False,True,-60.,wks)
   draw(plot)
   frame(wks)
   ;Plot PRECIP_NO3
   wks = gsn_open_wks (wks_type,"figs/PRECIP_NO3")
  ; plot=plot_field(PRECIP_NO3_annual_ensemble_mean,0.,1.,False,True,-60.,wks)
   draw(plot)
   frame(wks)

;Plot absolute difference
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_ABS")
   d=add_lev_lat_lon_arrays(PRECIP_O3_annual_ensemble_mean-PRECIP_NO3_annual_ensemble_mean,lev,lat,lon)
 ;  plot=plot_field(d,-0.1,0.1,False,True,-60.,wks)
   draw(plot)
   frame(wks)

;Plot relative difference
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_REL")
   d:=add_lev_lat_lon_arrays((PRECIP_O3_annual_ensemble_mean-PRECIP_NO3_annual_ensemble_mean)/PRECIP_O3_annual_ensemble_mean*100.,lev,lat,lon)
  ; plot=plot_field(d,-15.,15.,False,True,-60.,wks)
   draw(plot)
   frame(wks)

end
