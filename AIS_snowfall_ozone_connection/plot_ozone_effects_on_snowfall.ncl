load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=False
   process_monthly_data=False
plot_spatial_annual_snowfall=False
plot_precipitation_PDFs=False
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

O3SD ="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1986_2005/lens/"
NO3SD="/glade/p/cesm/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/files_1986_2005/no-ozone/"
nYears=20
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
iZ500=nVar
nVar=nVar+1

iO3=0
iNO3=1

if (LoadFreshData) then

   O3Data=new((/nYears,192,288,nEns,nSeason,nVar/),"double")
   NO3Data=new((/nYears,192,288,nEns,nSeason,nVar/),"double")
   O3DataMonthly =new((/nMonths,192,288,nEns,nVar/),"double")
   NO3DataMonthly=new((/nMonths,192,288,nEns,nVar/),"double")
      
   ;Load LENS data
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))
      print("Loading LENS: "+elong)
      in=addfile(str_concat((/O3SD,"year_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iPRECIP)=in->PRECIP
      in=addfile(str_concat((/O3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iPRECIP)=in->PRECIP
      if (process_monthly_data) then
         in=addfile(str_concat((/O3SD,"mon_PRECIP_",elong,".nc"/)),"r")
         O3DataMonthly(:,:,:,e,iPRECIP)=in->PRECIP
      end if
      in=addfile(str_concat((/O3SD,"year_PSL_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iPSL)=in->PSL
      in=addfile(str_concat((/O3SD,"DJF_PSL_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iPSL)=in->PSL
      
      in=addfile(str_concat((/O3SD,"year_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iTMQ)=in->TMQ       
      in=addfile(str_concat((/O3SD,"DJF_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iTMQ)=in->TMQ
      
      in=addfile(str_concat((/O3SD,"year_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iTMQ)=in->TMQ       
      in=addfile(str_concat((/O3SD,"DJF_TMQ_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iTMQ)=in->TMQ
      
      in=addfile(str_concat((/O3SD,"year_Z500_850_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iAnn,iZ500)=in->Z500       
      in=addfile(str_concat((/O3SD,"DJF_Z500_850_",elong,".nc"/)),"r")
      O3Data(:,:,:,e,iDJF,iZ500)=in->Z500

   end do

   ;Load no-ozone data
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))
      print("Loading no-ozone: "+elong)
      in=addfile(str_concat((/NO3SD,"year_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iPRECIP)=in->PRECIP
      in=addfile(str_concat((/NO3SD,"DJF_PRECIP_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iPRECIP)=in->PRECIP
      if (process_monthly_data) then
         in=addfile(str_concat((/NO3SD,"mon_PRECIP_",elong,".nc"/)),"r")
         NO3DataMonthly(:,:,:,e,iPRECIP)=in->PRECIP
      end if
      
      in=addfile(str_concat((/NO3SD,"year_PSL_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iPSL)=in->PSL
      in=addfile(str_concat((/NO3SD,"DJF_PSL_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iPSL)=in->PSL
      
      in=addfile(str_concat((/NO3SD,"year_TMQ_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iTMQ)=in->TMQ       
      in=addfile(str_concat((/NO3SD,"DJF_TMQ_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iTMQ)=in->TMQ
      
      in=addfile(str_concat((/NO3SD,"year_Z500_850_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iAnn,iZ500)=in->Z500       
      in=addfile(str_concat((/NO3SD,"DJF_Z500_850_",elong,".nc"/)),"r")
      NO3Data(:,:,:,e,iDJF,iZ500)=in->Z500          
   end do

   ;Calculate ensemble mean fields for each year.  Dimensions: [time, lat, lon, average period, variable]
   O3DataEnsMean =dim_avg_n   (O3Data,dEns) 
   NO3DataEnsMean=dim_avg_n   (NO3Data,dEns)

   ;Calculate ensemble mean climatologies.  Dimensions: [lat, lon, average period, variable]
   O3DataEnsMeanClimo =dim_avg_n     (O3DataEnsMean,dTime)
   NO3DataEnsMeanClimo=dim_avg_n     (NO3DataEnsMean,dTime)
   
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
   O3TS=new((/nYears,nEns,nSeason/),"double")
   O3TS(:,:,:)=0.0
   NO3TS=new((/nYears,nEns,nSeason/),"double")
   NO3TS(:,:,:)=0.0   
   O3TSMonthly =new((/nMonths,nEns/),"double")
   O3TSMonthly (:,:)=0.0
   NO3TSMonthly =new((/nMonths,nEns/),"double")
   NO3TSMonthly (:,:)=0.0

   do e=0,nEns-1
      print("Calculating integrated timeseries for LENS: "+(e+1))
      O3TS(:,e,iAnn)=calculate_integrated_timeseries(O3Data(:,:,:,e,iAnn,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      O3TS(:,e,iDJF)=calculate_integrated_timeseries(O3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      if (process_monthly_data) then
        print("...monthly timeseries...")      
        O3TSMonthly(:,e)=calculate_integrated_timeseries(O3DataMonthly(:,:,:,e,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      end if

      print("Calculating integrated timeseries for no-ozone: "+(e+1))
      NO3TS(:,e,iAnn)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iAnn,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      NO3TS(:,e,iDJF)=calculate_integrated_timeseries(NO3Data(:,:,:,e,iDJF,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      if (process_monthly_data) then
         print("...monthly timeseries...")
         NO3TSMonthly(:,e)=calculate_integrated_timeseries(NO3DataMonthly(:,:,:,e,iPRECIP),LandFrac,LandArea,AISMaskBinary)
      end if      
   end do   
   
   print("Saving data...")
   system("rm -f "+OutputFileName)
   ncdf=addfile(OutputFileName ,"c")

   ncdf->O3DataEnsMean=O3DataEnsMean
   ncdf->NO3DataEnsMean=NO3DataEnsMean
   ncdf->O3DataEnsMeanClimo=O3DataEnsMeanClimo  
   ncdf->NO3DataEnsMeanClimo=NO3DataEnsMeanClimo 
   ncdf->dEnsMeanClimo_pval=dEnsMeanClimo_pval 

   O3TS=O3TS/1.e12
   NO3TS=NO3TS/1.e12
   O3TSMonthly=O3TSMonthly/1.e12
   NO3TSMonthly=NO3TSMonthly/1.e12

   ncdf->O3TS=O3TS
   ncdf->NO3TS=NO3TS
   ncdf->O3TSMonthly=O3TSMonthly
   ncdf->NO3TSMonthly=NO3TSMonthly
      
   print("...done.")  
  
else
   
   print("Loading data...")
   ncdf=addfile(OutputFileName,"r")

   O3DataEnsMean=ncdf->O3DataEnsMean
   NO3DataEnsMean=ncdf->NO3DataEnsMean
   O3DataEnsMeanClimo=ncdf->O3DataEnsMeanClimo  
   NO3DataEnsMeanClimo=ncdf->NO3DataEnsMeanClimo 
   dEnsMeanClimo_pval=ncdf->dEnsMeanClimo_pval 

   O3TS=ncdf->O3TS
   NO3TS=ncdf->NO3TS
   O3TSMonthly=ncdf->O3TSMonthly
   NO3TSMonthly=ncdf->NO3TSMonthly
   print("...done.")   

end if

;Calculate absolute/relative spatial differences and grid-point difference significance
dEnsMeanClimoAbs=O3DataEnsMeanClimo-NO3DataEnsMeanClimo  ;Dimensions: [lat, lon, average period, variable]
dEnsMeanClimoRel=dEnsMeanClimoAbs/O3DataEnsMeanClimo

;Calculate monthly timeseries of ensemble statistics.

nStats=0
iMean=nStats
nStats=nStats+1
iMax=nStats
nStats=nStats+1
iMin=nStats
nStats=nStats+1

Years=fspan(1,nYears,nYears)+2005.-nYears
monthsinDecimalYears=fspan(1,nMonths,nMonths)/12.+2005.-nYears

MonthlyEnsembleMeanStats= new((/nMonths,nStats,2/),double)
MonthlyEnsembleMeanStats(:,iMin,iO3)=dim_min_n(O3TSMonthly,1)
MonthlyEnsembleMeanStats(:,iMax,iO3)=dim_max_n(O3TSMonthly,1)
MonthlyEnsembleMeanStats(:,iMean,iO3)=dim_avg_n(O3TSMonthly,1)
MonthlyEnsembleMeanStats(:,iMin,iNO3)=dim_min_n(NO3TSMonthly,1)
MonthlyEnsembleMeanStats(:,iMax,iNO3)=dim_max_n(NO3TSMonthly,1)
MonthlyEnsembleMeanStats(:,iMean,iNO3)=dim_avg_n(NO3TSMonthly,1)

;Calculate seasonal/annual timeseries of ensemble statistics
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

;Calculate monthly climatologies, for all ensemble members.  Dimensions: time, ensemble number
MonthlyClimoEnsembleMeanStats=new((/12,6/),double) ;months-in-year,ensemble
MonthlyClimo=rm_single_dims( clmMonTLL(reshape(O3TSMonthly, (/nMonths,nEns,1/) )) )

;Get ensemble-mean monthly climatology stats for integrated snowfall flux
MonthlyClimoEnsembleMeanStats(:,0)=dim_avg_n( MonthlyClimo , 1 )
MonthlyClimoEnsembleMeanStats(:,1)=dim_avg_n( MonthlyClimo , 1 )-dim_stddev_n( MonthlyClimo , 1 )
MonthlyClimoEnsembleMeanStats(:,2)=dim_avg_n( MonthlyClimo , 1 )+dim_stddev_n( MonthlyClimo , 1 )

MonthlyClimo:=rm_single_dims( clmMonTLL(reshape(NO3TSMonthly, (/nMonths,nEns,1/) )) )
MonthlyClimoEnsembleMeanStats(:,3)=dim_avg_n( MonthlyClimo , 1 )
MonthlyClimoEnsembleMeanStats(:,4)=dim_avg_n( MonthlyClimo , 1 )-dim_stddev_n( MonthlyClimo , 1 )
MonthlyClimoEnsembleMeanStats(:,5)=dim_avg_n( MonthlyClimo , 1 )+dim_stddev_n( MonthlyClimo , 1 )


;Load Brooke's reconstructed differences
    ;f=addfile("../input_data/Medley_reconstructions/final_recon_v1_cesm.nc","r")
    ;dPRECIP_Recon=new((/192,288/),double)
    ;dPRECIP_Recon_pval=new((/192,288/),double)
    ;tmp=lonPivot(f->rel_change_1986_2000(lat|:,lon|:),0.)
    ;tmp2=lonPivot(f->pval_1986_2000(lat|:,lon|:),0.)
    ;do lati=0,42
    ;   do loni=0,287
    ;      if (tmp(lati,loni).ge.-100.) then
    ;	 if (tmp(lati,loni).le.100.) then
    ;	    dPRECIP_Recon(lati,loni)=tmp(lati,loni)
    ;	    dPRECIP_Recon_pval(lati,loni)=tmp2(lati,loni)
    ;	 end if
    ;      end if
    ;   end do
    ;end do


;;;;;;;;;;;;;;;;;;;;;;;;:PLOTTING BEGINS HERE;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


wks_type=get_default_wks(); function in general_utilities.ncl


if (plot_spatial_snowfall_change) then
   
   latbound=-58.
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_REL")
   PRECIP_plot = plot_field(alll(dEnsMeanClimoRel(:,:,iDJF,iPRECIP),lev,lat,lon),-.2,.2,False,"MPL_BrBG",True,latbound,wks)
   significance_plot=plot_significance(alll(dEnsMeanClimo_pval(:,:,iDJF,iPRECIP),lev,lat,lon),wks)   
   PRECIP_plot = plot_basin_polygons(ZwallyBasinPolygonData,PRECIP_plot,wks,nAISMask,-9999)
   Z500_plot = plot_field(alll(dEnsMeanClimoAbs(:,:,iDJF,iZ500),lev,lat,lon),-25.,25.,True,"MPL_bwr",False,latbound,wks)
   getvalues Z500_plot
     "cnLevels" : levels
   end getvalues
   print(levels)
   overlay(PRECIP_plot,significance_plot)
   overlay(PRECIP_plot,Z500_plot)
   draw(PRECIP_plot)
   
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
