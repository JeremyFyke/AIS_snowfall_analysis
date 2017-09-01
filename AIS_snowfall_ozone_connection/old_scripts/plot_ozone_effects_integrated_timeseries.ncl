load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=True

plot_raw_timeseries=False

EnsNums=(/1,10,12,15,18,20,22,28/)
nEns=dimsizes(EnsNums)

;Load CAM grid data
in:=addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/CAM_base_data.nc","r")
lev=in->lev
lat=in->lat
lon=in->lon

LandArea = get_LandArea() 
LandFrac = get_LandFrac()
AISMask = get_AIS_mask()
AISMaskBinary=where(AISMask.gt.0,1,0)
nAISMask = toint(max(AISMask))

O3SD ="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1955_2005/lens/"
NO3SD="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1955_2005/no-ozone/"
nYears=51
OutputFileName="/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/Ozone_analysis_timeseries.nc"

;Hardcode averaging period indices
nSeason=0
iAnn=nSeason
nSeason=nSeason+1
iDJF=nSeason
nSeason=nSeason+1

Years=fspan(1,nYears,nYears)+2005.-nYears

iO3=0
iNO3=1

iRaw=0
iSmoothed=1

if (LoadFreshData) then
   O3Data=new((/nYears,192,288,nEns,nSeason/),"double")
   NO3Data=new((/nYears,192,288,nEns,nSeason/),"double")      
   O3TS=new((/nYears,nEns,nSeason/),"double")
   O3TS(:,:,:)=0.0
   NO3TS=new((/nYears,nEns,nSeason/),"double")
   NO3TS(:,:,:)=0.0 
   
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))

      ;Load LENS data      
      print("Loading LENS: "+elong)
      in=addfile(str_concat((/O3SD,"year_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn)=in->PRECIP
      in=addfile(str_concat((/O3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF)=in->PRECIP
      
      ;Load no-ozone data      
      print("Loading no-ozone: "+elong)
      in=addfile(str_concat((/NO3SD,"year_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn)=in->PRECIP
      in=addfile(str_concat((/NO3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF)=in->PRECIP

      print("Calculating integrated timeseries for: "+elong)
      ;Calculate AIS-integrated precipitation Annual and DJFtime series     
      O3TS(:,e,iAnn)=calculate_integrated_timeseries(O3Data(:,:,:,e,iAnn),LandFrac,LandArea,AISMaskBinary)
      O3TS(:,e,iDJF)=calculate_integrated_timeseries(O3Data(:,:,:,e,iDJF),LandFrac,LandArea,AISMaskBinary)   
      NO3TS(:,e,iAnn)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iAnn),LandFrac,LandArea,AISMaskBinary)
      NO3TS(:,e,iDJF)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iDJF),LandFrac,LandArea,AISMaskBinary)
   end do

   O3TS=O3TS/1.e12
   NO3TS=NO3TS/1.e12

   nStats=0
   iMean=nStats
   nStats=nStats+1
   iMax=nStats
   nStats=nStats+1
   iMin=nStats
   nStats=nStats+1
   iSmooth=nStats
   nStats=nStats+1 
   
   SmoothPer=11
   SmoothPerMid=SmoothPer/2 ;Note integer division, which truncates decimal part.

   ;Calculate seasonal/annual timeseries of ensemble statistics
   EnsembleStats=new((/nYears,nStats,2,2/),double); time, stats, ensemble, averaging period 
   EnsembleStats(:,iMin,iO3,iAnn)=dim_min_n(O3TS(:,:,iAnn),1)
   EnsembleStats(:,iMax,iO3,iAnn)=dim_max_n(O3TS(:,:,iAnn),1)
   EnsembleStats(:,iMean,iO3,iAnn)=dim_avg_n(O3TS(:,:,iAnn),1)
   EnsembleStats(:,iSmooth,iO3,iAnn)=runave(EnsembleStats(:,iMean,iO3,iAnn),SmoothPer,0)   
   
   EnsembleStats(:,iMin,iNO3,iAnn)=dim_min_n(NO3TS(:,:,iAnn),1)
   EnsembleStats(:,iMax,iNO3,iAnn)=dim_max_n(NO3TS(:,:,iAnn),1)
   EnsembleStats(:,iMean,iNO3,iAnn)=dim_avg_n(NO3TS(:,:,iAnn),1)
   EnsembleStats(:,iSmooth,iNO3,iAnn)=runave(EnsembleStats(:,iMean,iNO3,iAnn),SmoothPer,0)
   
   EnsembleStats(:,iMin,iO3,iDJF)=dim_min_n(O3TS(:,:,iDJF),1)
   EnsembleStats(:,iMax,iO3,iDJF)=dim_max_n(O3TS(:,:,iDJF),1)
   EnsembleStats(:,iMean,iO3,iDJF)=dim_avg_n(O3TS(:,:,iDJF),1)
   EnsembleStats(:,iSmooth,iO3,iDJF)=runave(EnsembleStats(:,iMean,iO3,iDJF),SmoothPer,0)   
   
   EnsembleStats(:,iMin,iNO3,iDJF)=dim_min_n(NO3TS(:,:,iDJF),1)
   EnsembleStats(:,iMax,iNO3,iDJF)=dim_max_n(NO3TS(:,:,iDJF),1)
   EnsembleStats(:,iMean,iNO3,iDJF)=dim_avg_n(NO3TS(:,:,iDJF),1)
   EnsembleStats(:,iSmooth,iNO3,iDJF)=runave(EnsembleStats(:,iMean,iNO3,iDJF),SmoothPer,0)   
   
   ;Calculate ozone minus no-ozone difference time series for all averaging periods
   DEnsembleStats=EnsembleStats(:,:,iO3,:)-EnsembleStats(:,:,iNO3,:)

   DEnsembleStats_pval=new((/nYears,nSeason/),double)
   tEmerged=new((/nSeason/),integer)
   do s=0,nSeason-1
     tEmergedFound=False
     do t=SmoothPerMid,nYears-SmoothPerMid-1
	ts=t-SmoothPerMid
	te=t+SmoothPerMid
	aveX=avg(O3TS(ts:te,:,s))
	varX=variance(O3TS(ts:te,:,s))
	aveY=avg(NO3TS(ts:te,:,s))
	varY=variance(NO3TS(ts:te,:,s))
	sX=SmoothPer*nEns
	sY=SmoothPer*nEns
	DEnsembleStats_pval(t,s)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)
	;Identify first time of significant emergence in smoothed time series
	if (DEnsembleStats_pval(t,s).lt.0.01 .and. .not. tEmergedFound) then
	   tEmerged(s)=t
	   tEmergedFound=True
	end if
     end do
     if (tEmergedFound)
        print(tEmerged(s))
     else
        print("No emergence of significant difference found, stopping so you can figure things out.")
	exit
     end if
   end do
   
end if

;Hard-code simple estimate of recent AIS mass loss (Gt/yr, Forsberg et al., 2017, 10.1007/s10712-016-9398-7) for plot comparisons
ForsbergObsStats=new((/3,4/),double)
ForsbergObsStats(0,:)=fspan(95.,95.,4) ;mean
ForsbergObsStats(1,:)=fspan(95.-50.,95.-50.,4) ;min
ForsbergObsStats(2,:)=fspan(95.+50.,95.+50.,4) ;max
ForsbergTime=fspan(2002.,2005,4)


;;;;;;;;;;;;;;;;;;;;;;;;:PLOTTING BEGINS HERE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



wks_type=get_default_wks(); function in general_utilities.ncl

if (plot_raw_timeseries) then
   ;Plot monthy/annual raw ensemble-mean integrated AIS accumulation time series
   wks = gsn_open_wks (wks_type,"figs/PRECIP_O_TS");Plot timeseries and ensemble spread of integrated precip
   res = True                      ; plot mods desired
   res@trXMinF=min(Years)
   res@trXMaxF=max(Years)
   res@trYMinF=500.
   res@trYMaxF=4500.
   res@gsnDraw      = False  
   res@gsnFrame     = False  
   res@xyLineThicknesses = (/1.0, 1.0, 10.0, 1.0,  1.0, 10.0/) 
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
end if


