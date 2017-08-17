load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=True

plot_raw_timeseries=False
plot_annual_timeseries_difference=False
plot_DJF_timeseries_difference=False
plot_spatial_annual_snowfall=False
plot_climatological_annual_cycle=False
plot_precipitation_PDFs=False
plot_basin_aggreggated_snowfall_change=False
plot_spatial_snowfall_change=True
plot_DJF_climatology_difference=False
plot_spatial_climatological_dPSL=False

;Load CAM grid data
in:=addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/CAM_base_data.nc","r")
lev=in->lev
lat=in->lat
lon=in->lon 

year2sec=31557600.
LandArea = get_LandArea() 
LandFrac = get_LandFrac()
AISMask = get_AIS_mask()
AISMaskBinary=where(AISMask.gt.0,1,0)
nAISMask = toint(max(AISMask))

ZwallyBasinPolygonData=readAsciiTable("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/AIS_basin_masks/Ant_Full_DrainageSystem_Polygons.txt",3,"float",7)

O3SD ="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1991_2005/lens/"
NO3SD="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1991_2005/no-ozone/"
nYears=15
OutputFileName="/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/Ozone_analysis.nc"

nMonths=nYears*12

EnsNums=(/1,10,12,15,18,20,22,28/)
nEns=dimsizes(EnsNums)

;Hardcode dimension numbers for initial full datasets
dTime=0
dLat=1
dLon=2
dEns=3
dAvgPer=4
dVar=5

;Hardcode averaging period indices
nSeason=0
iAnn=nSeason
nSeason=nSeason+1
iDJF=nSeason
nSeason=nSeason+1

;Hardcode variable indices
nVar=0
iPRECIP=nVar
nVar=nVar+1
iPSL=nVar
nVar=nVar+1
iTMQ=nVar
nVar=nVar+1

iO3=0
iNO3=1

if (LoadFreshData) then

   O3Data=new((/nYears,192,288,nEns,nSeason,nVar/),"double")
   NO3Data=new((/nYears,192,288,nEns,nSeason,nVar/),"double")
   
   ;Load LENS data
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))
      print("Loading LENS: "+elong)
      in=addfile(str_concat((/O3SD,"year_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iPRECIP)=in->PRECIP
      in=addfile(str_concat((/O3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iPRECIP)=in->PRECIP     

      in=addfile(str_concat((/O3SD,"year_PSL_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iPSL)=in->PSL
      in=addfile(str_concat((/O3SD,"DJF_PSL_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iPSL)=in->PSL
      
      in=addfile(str_concat((/O3SD,"year_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iTMQ)=in->TMQ       
      in=addfile(str_concat((/O3SD,"DJF_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iTMQ)=in->TMQ          
   end do

   ;Load no-ozone data
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))
      print("Loading no-ozone: "+elong)
      in=addfile(str_concat((/NO3SD,"year_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iPRECIP)=in->PRECIP
      in=addfile(str_concat((/NO3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iPRECIP)=in->PRECIP

      in=addfile(str_concat((/NO3SD,"year_PSL_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iPSL)=in->PSL
      in=addfile(str_concat((/NO3SD,"DJF_PSL_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iPSL)=in->PSL
      
      in=addfile(str_concat((/NO3SD,"year_TMQ_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iTMQ)=in->TMQ       
      in=addfile(str_concat((/NO3SD,"DJF_TMQ_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iTMQ)=in->TMQ        
   end do

   ;Calculate ensemble mean and STD map-view fields for each year.  Dimensions: [time, lat, lon, average period, variable]
   O3DataEnsMean =dim_avg_n   (O3Data,dEns)
   O3DataEnsStd  =dim_stddev_n(O3Data,dEns)  
   NO3DataEnsMean=dim_avg_n   (NO3Data,dEns)
   NO3DataEnsStd =dim_stddev_n(NO3Data,dEns)

   ;Calculate ensemble mean climatologies.  Dimensions: [lat, lon, average period, variable]
   O3DataEnsMeanClimo =dim_avg_n     (O3DataEnsMean,dTime)
   NO3DataEnsMeanClimo=dim_avg_n     (NO3DataEnsMean,dTime)
   
   ;Calculate absolute/relative spatial differences and grid-point difference significance
   dEnsMeanClimoAbs=O3DataEnsMeanClimo-NO3DataEnsMeanClimo  ;Dimensions: [lat, lon, average period, variable]
   dEnsMeanClimoRel=dEnsMeanClimoAbs/O3DataEnsMeanClimo
   
   ;Can't figure out how to do this as a single call, so looping it...
   dEnsMeanClimo_pval=new((/192,288,nSeason,nVar/),double)
   do lati=0,80 ; don't do whole globe, to save processing time.
      print("Processing significance for latitude"+tostring(lati))
      do loni=0,287
	 do peri=0,nSeason-1
	    do vari=0,nVar-1      
	       aveX=avg(O3Data      (:,lati,loni,:,peri,vari))
	       varX=variance(O3Data (:,lati,loni,:,peri,vari))
	       sX=nYears*nEns
	       aveY=avg(NO3Data     (:,lati,loni,:,peri,vari))
	       varY=variance(NO3Data(:,lati,loni,:,peri,vari))
	       sY=nYears*nEns
	       dEnsMeanClimo_pval(lati,loni,peri,vari)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)
	    end do
	 end do      
      end do
   end do
   
   ;Calculate AIS-integrated precipitation Annual and DJFtime series
   O3TS=new((/nYears,nEns,2/),"double")
   O3TS(:,:,:)=0.0
   do e=0,nEns-1
      print("Calculating integrated timeseries for LENS: "+(e+1))
      O3TS(:,e,iAnn)=calculate_integrated_timeseries(O3Data(:,:,:,e,iAnn,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      O3TS(:,e,iDJF)=calculate_integrated_timeseries(O3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,AISMaskBinary)
   end do
   
   NO3TS=new((/nYears,nEns,2/),"double")
   NO3TS(:,:,:)=0.0
   do e=0,nEns-1
      print("Calculating integrated timeseries for no-ozone: "+(e+1))
      NO3TS(:,e,iAnn)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iAnn,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      NO3TS(:,e,iDJF)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,AISMaskBinary)
   end do   

;   print("Saving data...")
;   system("rm -f "+OutputFileName)
;   ncdf=addfile(OutputFileName ,"c")
;
;   ;ncdf->O3Data=O3Data
;   ;ncdf->NO3Data=NO3Data
;
;   ncdf->O3DataEnsMean=O3DataEnsMean
;   ncdf->O3DataEnsStd=O3DataEnsStd
;   ncdf->NO3DataEnsMean=NO3DataEnsMean
;   ncdf->NO3DataEnsStd=NO3DataEnsStd
;   
;   ncdf->O3DataEnsMeanClimo=O3DataEnsMeanClimo
;   ncdf->O3DataEnsVarClimo=O3DataEnsVarClimo
;   ncdf->NO3DataEnsMeanClimo=NO3DataEnsMeanClimo   
;   ncdf->NO3DataEnsVarClimo=NO3DataEnsVarClimo
;
;   ncdf->O3TS=O3TS
;   ncdf->NO3TS=NO3TS
;   print("...done.")  
  
else
;   ncdf=addfile(OutputFileName,"r")
;   
;   ;O3Data=ncdf->O3Data
;   ;NO3Data=ncdf->NO3Data   
;   
;   O3DataEnsMean=ncdf->O3DataEnsMean
;   NO3DataEnsMean=ncdf->NO3DataEnsMean
;   O3DataEnsStd=ncdf->O3DataEnsStd
;   NO3DataEnsStd=ncdf->NO3DataEnsStd
;   O3DataEnsMeanClimo=ncdf->O3DataEnsMeanClimo
;   NO3DataEnsMeanClimo=ncdf->NO3DataEnsMeanClimo
;   O3DataEnsVarClimo=ncdf->O3DataEnsVarClimo
;   NO3DataEnsVarClimo=ncdf->NO3DataEnsVarClimo
;   O3TS=ncdf->O3TS
;   NO3TS=ncdf->NO3TS
;   ds=dimsizes(O3TS)
;   nYears=ds(0)
end if

;;;Analyze/plot stuff

;Calculate ensemble-mean stuff.  Store in arrays that are amenable for multi-line plotting.
Years=fspan(1,nYears,nYears)+2005.-nYears

;monthsinDecimalYears=fspan(1,nMonthClimatology,nMonthClimatology)/12.+2005.-nYearClimatology
;MonthlyEnsembleMeanStats= new((/6,nMonthClimatology/),double)
;MonthlyEnsembleMeanStats(0,:)=dim_min_n(MonthlyO3TS,1)
;MonthlyEnsembleMeanStats(1,:)=dim_max_n(MonthlyO3TS,1)
;MonthlyEnsembleMeanStats(2,:)=dim_avg_n(MonthlyO3TS,1)

;MonthlyEnsembleMeanStats(3,:)=dim_min_n(MonthlyNO3TS,1)
;MonthlyEnsembleMeanStats(4,:)=dim_max_n(MonthlyNO3TS,1)
;MonthlyEnsembleMeanStats(5,:)=dim_avg_n(MonthlyNO3TS,1)

print( "Average dynamic ozone annual mass flux="+avg(O3TS(:,:,iAnn))/1.e12 )
print( "Average fixed ozone annual mass flux="  +avg(NO3TS(:,:,iAnn))/1.e12 )

nStats=0
iMean=nStats
nStats=nStats+1
iMax=1
nStats=nStats+1
iMin=2
nStats=nStats+1
iOzoneDiff=nStats
nStats=nStats+1

;Calculate timeseries of ensemble-average statistics
EnsembleStats=new((/nYears,nStats,2,2/),double); time, stats, ensemble, averaging period 
EnsembleStats(:,iMin,iO3,iAnn)=dim_min_n(O3TS(:,:,iAnn),1)
EnsembleStats(:,iMax,iO3,iAnn)=dim_max_n(O3TS(:,:,iAnn),1)
EnsembleStats(:,iMean,iO3,iAnn)=dim_avg_n(O3TS(:,:,iAnn),1)
EnsembleStats(:,iMin,iNO3,iAnn)=dim_min_n(NO3TS(:,:,iAnn),1)
EnsembleStats(:,iMax,iNO3,iAnn)=dim_max_n(NO3TS(:,:,iAnn),1)
EnsembleStats(:,iMean,iNO3,iAnn)=dim_avg_n(NO3TS(:,:,iAnn),1)
EnsembleStats(:,iMin,iO3,iDJF)=dim_min_n(O3TS(:,:,iDJF),1)
EnsembleStats(:,iMax,iO3,iDJF)=dim_max_n(O3TS(:,:,iDJF),1)
EnsembleStats(:,iMean,iO3,iDJF)=dim_avg_n(O3TS(:,:,iDJF),1)
EnsembleStats(:,iMin,iNO3,iDJF)=dim_min_n(NO3TS(:,:,iDJF),1)
EnsembleStats(:,iMax,iNO3,iDJF)=dim_max_n(NO3TS(:,:,iDJF),1)
EnsembleStats(:,iMean,iNO3,iDJF)=dim_avg_n(NO3TS(:,:,iDJF),1)

;Calculate ozone minus no-ozone difference time series
DAnnualEnsembleStats=new((/nYears,nStats,2/),double); time, stats, averaging period
DAnnualEnsembleStats(:,iMin,iAnn) =EnsembleStats(:,iMin,iO3,iAnn) -EnsembleStats(:,iMin,iNO3,iAnn)
DAnnualEnsembleStats(:,iMax,iAnn) =EnsembleStats(:,iMax,iO3,iAnn) -EnsembleStats(:,iMax,iNO3,iAnn)
DAnnualEnsembleStats(:,iMean,iAnn)=EnsembleStats(:,iMean,iO3,iAnn)-EnsembleStats(:,iMean,iNO3,iAnn)
DAnnualEnsembleStats(:,iMin,iDJF) =EnsembleStats(:,iMin,iO3,iDJF) -EnsembleStats(:,iMin,iNO3,iDJF)
DAnnualEnsembleStats(:,iMax,iDJF) =EnsembleStats(:,iMax,iO3,iDJF) -EnsembleStats(:,iMax,iNO3,iDJF)
DAnnualEnsembleStats(:,iMean,iDJF)=EnsembleStats(:,iMean,iO3,iDJF)-EnsembleStats(:,iMean,iNO3,iDJF)

;;Here, just calculate post-depletion climatologies.  I.e. 15==1991->2005
;printVarSummary(MonthlyO3TS)
;print(nMonthClimatology-nMonthFinalComparisonPeriod)
;MonthlyClimoEnsembleMeanStats=new((/6,12/),double)
;MonthlyClimo=rm_single_dims( clmMonTLL(reshape(MonthlyO3TS(nMonthClimatology-nMonthFinalComparisonPeriod:nMonthClimatology-1,:) , (/nMonthFinalComparisonPeriod,nLENS,1/))) )
;MonthlyClimoEnsembleMeanStats(0,:)=dim_avg_n( MonthlyClimo , 1 )
;MonthlyClimoEnsembleMeanStats(1,:)=MonthlyClimoEnsembleMeanStats(0,:)-dim_stddev_n( MonthlyClimo , 1 )
;MonthlyClimoEnsembleMeanStats(2,:)=MonthlyClimoEnsembleMeanStats(0,:)+dim_stddev_n( MonthlyClimo , 1 )
;MonthlyClimo:=rm_single_dims( clmMonTLL(reshape(MonthlyNO3TS(nMonthClimatology-nMonthFinalComparisonPeriod:nMonthClimatology-1,:) , (/nMonthFinalComparisonPeriod,nEns,1/))) )
;MonthlyClimoEnsembleMeanStats(3,:)=dim_avg_n( MonthlyClimo , 1 )
;MonthlyClimoEnsembleMeanStats(4,:)=MonthlyClimoEnsembleMeanStats(3,:)-dim_stddev_n( MonthlyClimo , 1 )
;MonthlyClimoEnsembleMeanStats(5,:)=MonthlyClimoEnsembleMeanStats(3,:)+dim_stddev_n( MonthlyClimo , 1 )

;Hard-code simple estimate of recent AIS mass loss (Gt/yr, Forsberg et al., 2017, 10.1007/s10712-016-9398-7) for plot comparisons
ForsbergObsStats=new((/3,4/),double)
ForsbergObsStats(0,:)=fspan(95.,95.,4) ;mean
ForsbergObsStats(1,:)=fspan(95.-50.,95.-50.,4) ;min
ForsbergObsStats(2,:)=fspan(95.+50.,95.+50.,4) ;max
ForsbergTime=fspan(2002.,2005,4)

;Load Brooke's reconstructed differences
f=addfile("../input_data/Medley_reconstructions/final_recon_v1_cesm.nc","r")
dPRECIP_Recon=new((/192,288/),double)
dPRECIP_Recon_pval=new((/192,288/),double)

tmp=lonPivot(f->rel_change_1986_2000(lat|:,lon|:),0.)
tmp2=lonPivot(f->pval_1986_2000(lat|:,lon|:),0.)
do lati=0,42
   do loni=0,287
      if (tmp(lati,loni).ge.-100.) then
	 if (tmp(lati,loni).le.100.) then
	    dPRECIP_Recon(lati,loni)=tmp(lati,loni)
	    dPRECIP_Recon_pval(lati,loni)=tmp2(lati,loni)
	 end if
      end if
   end do
end do

if (plot_basin_aggreggated_snowfall_change) then
   ;Get CESM and reconstruction climatological accumulation changes on basin-scale
   Basin_integrated_dSNOWFALL_O3  = new((/nYears,nEns,nAISMask/),double)
   Basin_integrated_dSNOWFALL_NO3 = new((/nYears,nEns,nAISMask/),double)
   dPRECIP_CESM_basin_pval = new((/nAISMask/),double)
   dPRECIP_CESM_basin = new((/nAISMask/),double)
   dPRECIP_CESM_mask=new(dimsizes(AISMask),"double")
   dPRECIP_CESM_pval_mask=new(dimsizes(AISMask),"double")
   do b = 0,nAISMask-1
      print("***")
      BasinMask=where(AISMask.eq.b+1,1,0)
      do e = 0,nEns-1
	 Basin_integrated_dSNOWFALL_O3(:,e,b)=calculate_integrated_timeseries(O3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,BasinMask)
	 Basin_integrated_dSNOWFALL_NO3(:,e,b)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,BasinMask)      
      end do
      aveX=avg(Basin_integrated_dSNOWFALL_O3(:,:,b))
      varX=variance(Basin_integrated_dSNOWFALL_O3(:,:,b))
      aveY=avg(Basin_integrated_dSNOWFALL_NO3(:,:,b))
      varY=variance(Basin_integrated_dSNOWFALL_NO3(:,:,b))
      dPRECIP_CESM_basin(b)=(aveX-aveY)/aveX ;Relative change
      dPRECIP_CESM_basin_pval(b)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)  

      dPRECIP_CESM_mask=where(AISMask.eq.b,dPRECIP_CESM_basin(b),dPRECIP_CESM_mask)
      dPRECIP_CESM_pval_mask=where(AISMask.eq.b,dPRECIP_CESM_basin_pval(b),dPRECIP_CESM_pval_mask)
   end do
end if


;;;;;;;;;;;;;;;;;;;;;;;;:PLOTTING BEGINS HERE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
wks_type=get_default_wks(); function in general_utilities.ncl

if (plot_basin_aggreggated_snowfall_change) then


     
end if

if (plot_spatial_snowfall_change) then
   
   wks = gsn_open_wks (wks_type,"figs/dTMQ_REL")
   PRECIP_plot = plot_field(alll(dEnsMeanClimoRel(:,:,iDJF,iTMQ),lev,lat,lon),-.15,.15,False,"cmp_b2r",True,-50.,wks)
   significance_plot=plot_significance(alll(dEnsMeanClimo_pval(:,:,iDJF,iTMQ),lev,lat,lon),wks)   
   PRECIP_plot = plot_basin_polygons(ZwallyBasinPolygonData,PRECIP_plot,wks,nAISMask,-9999)
   PSL_plot = plot_field(alll(dEnsMeanClimoAbs(:,:,iDJF,iPSL),lev,lat,lon),-1000.,10.,True,(/"blue","black","red"/),False,-50.,wks)
   overlay(PRECIP_plot,significance_plot)
   overlay(PRECIP_plot,PSL_plot)
   draw(PRECIP_plot)
   frame(wks)    
   
end if

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

if (plot_annual_timeseries_difference) then
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
end if

if (plot_DJF_timeseries_difference) then
   ;plot difference in integrated AIS accumulation due to ozone, overlay against estimated loss trend
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_DJF_O_TS");Plot timeseries and ensemble spread of integrated precip
   res := True
   res@trXMinF=min(Years)
   res@trXMaxF=max(Years)
   res@trYMinF=-150.
   res@trYMaxF=300.
   res@gsnDraw      = False                        ; don't draw yet
   res@gsnFrame     = False                        ; don't advance frame yet
   res@xyLineThicknesses := (/20.0/)          ; make second line thicker
   res@xyLineColors      := (/"black"/)          ; change line color   
   res@tiYAxisString = "DJF precipitation difference (Gt/yr)"
   res@tiXAxisString = "Year"
   difference_plot = gsn_csm_xy(wks,Years-0.5,DJFEnsembleMeanStats(6,:),res)

   res@xyLineColors      := (/"red"/)          ; change line color
   running_average_plot = gsn_csm_xy(wks,Years-0.5,DJFEnsembleMeanStats(7,:),res)

   ;res@xyLineThicknesses := (/20.0,  10.0, 10.0/)          ; make second line thicker
   ;res@xyDashPattern     := (/0,2,2/) 
   ;res@xyLineColors      := (/"orange","orange","orange"/)          ; change line color
   ;obs_plot = gsn_csm_xy(wks,ForsbergTime-0.5,ForsbergObsStats,res)
   ;lineres=True
   ;lineres@lgLineThicknesses=(/20.0, 20.0/)
   ;lineres@lgLineColors=(/"black","orange"/)
   ;gres = True
   ;gres@YPosPosition = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   ;gres@XPosPosition = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
   ;textres = True
   ;textres@lgLabels = (/"Ozone-attributable AIS accumulation change","2002-2015 AIS mass loss estimate"/) 
   ;difference_plot = simple_legend(wks,difference_plot,gres,lineres,textres)

   overlay(difference_plot,running_average_plot)
   ;overlay(difference_plot,obs_plot)
   draw(difference_plot)
   frame(wks)   
end if

if (plot_spatial_climatological_dPSL) then
   wks = gsn_open_wks (wks_type,"figs/dPSL_DJF")
   field=alll(O3DataEnsMeanClimo(:,:,iDJF,iPSL)-NO3DataEnsMeanClimo(:,:,iDJF,iPSL),lev,lat,lon)
   plot=plot_field(field,-250.,250.,True,(/"blue","black","red"/),True,-50.,wks)
   draw(plot)
   frame(wks)   
end if

if (plot_spatial_annual_snowfall) then
   ;Plot PRECIP_O3
   wks = gsn_open_wks (wks_type,"figs/PRECIP_O3")
   ;plot_field function: in plotting_utilities.ncl
   plot=plot_field(PRECIP_O3_annual_ensemble_mean,0.,1.,False,"rainbow",True,-60.,wks)
   draw(plot)
   frame(wks)
   ;Plot PRECIP_NO3
   wks = gsn_open_wks (wks_type,"figs/PRECIP_NO3")
   plot=plot_field(PRECIP_NO3_annual_ensemble_mean,0.,1.,False,"rainbow",True,-60.,wks)
   draw(plot)
   frame(wks)

   ;Plot absolute difference
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_ABS")
   d=alll(PRECIP_O3_annual_ensemble_mean-PRECIP_NO3_annual_ensemble_mean,lev,lat,lon)
   plot=plot_field(d,-0.1,0.1,False,"rainbow",True,-60.,wks)
   draw(plot)
   frame(wks)

   ;Plot relative difference
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_REL")
   d:=alll((PRECIP_O3_annual_ensemble_mean-PRECIP_NO3_annual_ensemble_mean)/PRECIP_O3_annual_ensemble_mean*100.,lev,lat,lon)
   plot=plot_field(d,-15.,15.,False,"rainbow",True,-60.,wks)
   draw(plot)
   frame(wks)
end if

if (plot_climatological_annual_cycle) then
   months=ispan(0,11,1)
   wks = gsn_open_wks (wks_type,"figs/Monthly_climatology");Plot timeseries and ensemble spread of integrated precip
   miny=1300.
   maxy=3100.
   res = True                      ; plot mods desired   
   res@gsnDraw          = False
   res@gsnFrame         = False
   res@xyLineThicknesses = (/  10.0,  5.0, 5.0, 10.0,  5.0, 5.0/) 
   res@xyLineColors      = (/"blue","blue","blue","red","red","red"/) 
   res@xyDashPattern     := (/0,2,2,0,2,2/)
   res@tmXBValues        = months
   res@tmXBLabels        = (/"J","A","S","O","N","D",\
                                "J","F","M","A","M","J"/)
   res@tmXBMode             = "Explicit"      ; label independently
   res@tiYAxisString = "Integrated accumulation (Gt/yr)"
   res@tiXAxisString = "Month"
   res@trXMaxF = 11
   res@trYMinF = miny
   res@trYMaxF = maxy      
   tmp=new(dimsizes(MonthlyClimoEnsembleMeanStats),double)
   tmp(:,0:5)=MonthlyClimoEnsembleMeanStats(:,6:11);this makes cycle start in July
   tmp(:,6:11)=MonthlyClimoEnsembleMeanStats(:,0:5);this makes cycle start in July

   monthly_climatology_plot = gsn_csm_xy(wks,months,tmp,res)
   pgres                  = True
   pgres@gsFillOpacityF   = 0.2 
   ;Draw box around emerged region
   EmergedIndex=where(tmp(3,:).gt.tmp(1,:),0,1)
   xp=new((/4/),float)
   yp=new((/4/),float)
   xp(0)=maxind(EmergedIndex);Find first index of emerged conditions
   xp(3)=xp(0)   
   xp(1)=11-maxind(EmergedIndex(::-1));Cryptically, find last index of emerged conditions
   xp(2)=xp(1)

   yp(0)=maxy
   yp(1)=yp(0)
   yp(2)=miny
   yp(3)=yp(2) 

   pgres@gsFillColor      = "Grey"
   dummy0 = gsn_add_polygon (wks,monthly_climatology_plot,xp,yp,pgres)
   ;Make shaded ranges
   xp:=new((/24/),float)
   yp:=new((/24/),double)
   xp(0:11)=months
   xp(12:23)=months(::-1)
   yp(0:11)=tmp(5,:)
   yp(12:23)=tmp(4,::-1)
   pgres@gsFillColor      = "Red"
   pgres@gsFillOpacityF        = 0.2  
   dummy1 = gsn_add_polygon (wks,monthly_climatology_plot,xp,yp,pgres)
   yp(0:11)=tmp(2,:)
   yp(12:23)=tmp(1,::-1)   
   pgres@gsFillColor      = "Blue"
   dummy2 = gsn_add_polygon (wks,monthly_climatology_plot,xp,yp,pgres)
   lineres=True
   lineres@lgLineThicknesses=(/20.0, 20.0/)
   lineres@lgLineColors=(/"blue","red"/)   
   gres = True
   gres@YPosPosition = 95.    ; expressed as %, 0->100, sets position of top border of legend 
   gres@XPosPosition = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
   textres = True
   textres@lgLabels = (/"Dynamic ozone ensemble","Fixed ozone ensemble"/) 
   monthly_climatology_plot = simple_legend(wks,monthly_climatology_plot,gres,lineres,textres)
   draw(monthly_climatology_plot)  
   frame(wks)
end if

if (plot_precipitation_PDFs) then
   wks = gsn_open_wks (wks_type,"figs/Accumulation_PDF");Plot timeseries and ensemble spread of integrated precip   
   res:=True
   res@gsnDraw          = False
   res@gsnFrame         = False
   res@xyLineThicknessF       = 5
   res@xyLineColor            = "blue" 
   plot0=gsn_csm_xy (wks,O3Annualpdf@bin_center,O3Annualpdf/sum(O3Annualpdf),res)
   res@xyLineColor            = "red" 
   plot1=gsn_csm_xy (wks,NO3Annualpdf@bin_center,NO3Annualpdf/sum(NO3Annualpdf),res)
   overlay(plot0,plot1)
   draw(plot0)
   frame(wks)   
end if
