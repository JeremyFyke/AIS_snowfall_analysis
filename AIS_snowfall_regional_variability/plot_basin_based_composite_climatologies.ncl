load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"
begin   

   plot_AIS_SMB_spatial_fields=0
   plot_AISMask=0
   plot_dSIC=0
   plot_dSIC_plus_dZ500=0
   plot_dMoisture_plus_SIC=0
   plot_dMoisture_dZ500_transport=0
   plot_non_annular_difference=0
   plot_dMoisture_plus_basin_correlation=0
   plot_dz500_plus_basin_correlation=0
   plot_dMoisture_plus_PE=0
   plot_evap_to_transport_ratio=0
   plot_spatial_correlation_CESM=0
   plot_matrix_of_CESM_interbasin_correlations=0
   plot_spatial_correlation_RACMO=0
   plot_matrix_of_RACMO_interbasin_correlations=0
   plot_spatial_correlation_PAGES=0
   plot_matrix_of_PAGES_interbasin_correlations=0 
   plot_seasonal_cycle_difference=0
   plot_non_annular_difference=0
   plot_ice_core_correlations=0
   plot_ice_core_variability=0  
   plot_basin_overview=0
   
   year2sec=31557600.
   wks_type=get_default_wks()   
   default_res = get_default_res()
   SD="./"

   ;Load CAM grid data
   in=addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/CAM_base_data.nc","r")
   P0=in->P0
   hyam=in->hyam
   hybm=in->hybm
   lev=in->lev
   lat=in->lat
   lon=in->lon 
   AISMask = get_AIS_mask()
   nAISMask = max(AISMask)
   Area=get_LandArea() ;m^2
   LandFrac = get_LandFrac()
   LandArea = get_LandArea()
   Topog = get_Topog()
   
   ;Load ice drainage boundary polygons
   ZwallyBasinPolygonData=readAsciiTable("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/AIS_basin_masks/Ant_Full_DrainageSystem_Polygons.txt",3,"float",7)    
       
   ;;;Load RACMO AIS data
   infile = addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/AIS_basin_masks/AIS_Full_basins_Zwally_ANT27grid.nc","r")
   AISMask_RACMO=infile->Zwallybasins
   infile = addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/smb_monthly_ANT27_ERAIN_r490_1979_2015.nc","r")
   PRECIP_annual_RACMO=month_to_annual(infile->smb,1)
   PRECIP_annual_RACMO=dtrend_n(PRECIP_annual_RACMO,False,0)
   PRECIP_climatology_RACMO=dim_avg_n(PRECIP_annual_RACMO,0)
   PRECIP_variability_RACMO=dim_stddev_n(PRECIP_annual_RACMO,0)
   PRECIP_CV_RACMO=PRECIP_variability_RACMO/where(PRECIP_climatology_RACMO.ne.0.0,PRECIP_climatology_RACMO,PRECIP_climatology_RACMO@_FillValue); avoid divide by zer
   
   ds_RACMO=dimsizes(PRECIP_annual_RACMO)
   nt_RACMO=ds_RACMO(0)
   infile = addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/maskers_RACMO2.nc","r")
   LandFrac_RACMO=infile->mask2d
   infile = addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/gridarea_RACMO2.nc","r")
   LandArea_RACMO=infile->cell_area ;m^2  
   time_RACMO=fspan(1,nt_RACMO,nt_RACMO)
   Basin_integrated_TS_RACMO = new((/nt_RACMO,toint(nAISMask)/),"double")
   do b = 1,nAISMask
      bm1=b-1
      BasinMask_RACMO=where(AISMask_RACMO.eq.b,1,0)
      Basin_integrated_TS_RACMO(:,bm1)=calculate_integrated_timeseries(PRECIP_annual_RACMO,LandArea_RACMO,LandFrac_RACMO*BasinMask_RACMO)
   end do

   ;Load ice core data, stack by basin
   out:=composite_PAGES_data(lat,lon,AISMask,nAISMask)
   PAGES_core_lat=out[0]
   PAGES_core_lon=out[1]
   time_PAGES=out[3]
   Basin_integrated_TS_PAGES=out[5]
      
   if (plot_AISMask.eq.1) then
      wks = gsn_open_wks (wks_type,"AISMask")	 	 
      plot  = plot_AISMask_overview(AISMask,wks,1,-60)
      plot  = plot_core_locations(plot,wks,PAGES_core_lat,PAGES_core_lon) 	 
      draw(plot)
      frame(wks)
   end if
   
   if (plot_AIS_SMB_spatial_fields.eq.1)
      AISMaskBinary=where(AISMask.gt.0,1,0)
      in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
      lat=in->lat
      lon=in->lon
      SNOW_annual = (in->SNOW) * year2sec / 1000.    ;mm/s -> m/yr
      SNOW_annual=replace_mask_values(SNOW_annual,0.)
      SNOW_AIS_integrated=calculate_integrated_timeseries(SNOW_annual,LandArea,LandFrac*AISMaskBinary)
      SNOW_AIS_integrated=SNOW_AIS_integrated / 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr
      RAIN_annual = (in->RAIN) * year2sec / 1000.    ;mm/s -> m/yr
      RAIN_annual=replace_mask_values(RAIN_annual,0.)
      RAIN_AIS_integrated=calculate_integrated_timeseries(RAIN_annual,LandArea,LandFrac*AISMaskBinary)
      RAIN_AIS_integrated=RAIN_AIS_integrated / 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr

      PRECIP_annual=SNOW_annual+RAIN_annual
      PRECIP_annual@_FillValue=default_fillvalue(typeof(PRECIP_annual))
      PRECIP_annual=where(PRECIP_annual.eq.0.,PRECIP_annual@_FillValue,PRECIP_annual)
      PRECIP_climatology=mask(dim_avg_n(PRECIP_annual,0),AISMaskBinary,1)
      
      PRECIP_climatology=alll(PRECIP_climatology,lev,lat,lon)
      PRECIP_variability=mask(dim_stddev_n(PRECIP_annual,0),AISMaskBinary,1)
      PRECIP_variability=alll(PRECIP_variability,lev,lat,lon)
      PRECIP_AIS_integrated=RAIN_AIS_integrated+SNOW_AIS_integrated
      PRECIP_AIS_integrated_avg=avg(PRECIP_AIS_integrated)
      PRECIP_AIS_integrated_std=stddev(PRECIP_AIS_integrated)
      PRECIP_CV=PRECIP_variability/PRECIP_climatology
      PRECIP_CV=alll(PRECIP_CV,lev,lat,lon)
      
      print("CV averaged from point-wise CVs="+wgt_areaave2(PRECIP_CV,Area,0))
      print("Climatological AIS integrated mean snow:")
      print(avg(SNOW_AIS_integrated))
      print("Climatological AIS integrated mean rain:")
      print(avg(RAIN_AIS_integrated))      
      print("Total integrated AIS precip (average/std:")
      print(PRECIP_AIS_integrated_avg)
      print(PRECIP_AIS_integrated_std)   

      wks=gsn_open_wks(wks_type,"figs/AIS_mean_precip_map")
      cnLevels = array_append_record(fspan(0.0,0.3,31),fspan(0.4,1.2,9),0)
      AIS_climatology_plot = plot_precip_climatology(PRECIP_climatology,cnLevels,wks,1)
      AIS_climatology_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_climatology_plot,wks,nAISMask,-9999)
      draw (AIS_climatology_plot)
      frame(wks)      

      wks=gsn_open_wks(wks_type,"figs/AIS_std_precip_map")
      cnLevels := array_append_record(fspan(0.0,0.04,21),fspan(0.05,0.22,18),0)
      AIS_variability_plot = plot_precip_climatology(PRECIP_variability,cnLevels,wks,1)
      AIS_variability_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_variability_plot,wks,nAISMask,-9999)      
      draw (AIS_variability_plot)
      frame(wks)
      
      wks=gsn_open_wks(wks_type,"figs/AIS_CV_map")
      cnLevels := fspan(0.1,0.37,37)      
      AIS_CV_plot = plot_precip_climatology(PRECIP_CV,cnLevels,wks,1)
      AIS_CV_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_CV_plot,wks,nAISMask,-9999)      
      draw (AIS_CV_plot)
      frame(wks)
      
      ;wks=gsn_open_wks(wks_type,"figs/AIS_CV_map_RACMO")
      ;cnLevels := fspan(0.1,0.37,37)      
      ;AIS_CV_plot = plot_precip_climatology(PRECIP_CV_RACMO,cnLevels,wks,0)     
      ;draw (AIS_CV_plot)
      ;frame(wks)            
   end if
   
   ;Load CVDP diagnostics for CESM LE run
   in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/CVDP/Controls/CESM1_LENS_Coupled_Control.cvdp_data.401-2200.nc","r")
   nino34=month_to_annual(in->nino34,1)
   sam=in->sam_pc_ann
   psa1=in->psa1_pc_ann
   psa2=in->psa2_pc_ann
   sam_eof=-(in->sam_ann)*100. ; convert to Pa from hPa, flip sign
   psa1_eof=-(in->psa1_sst_regression_ann) ;
   lat_CVDP=in->lat
   lon_CVDP=in->lon
   sam_eof=alll(sam_eof,lev,lat_CVDP,lon_CVDP)
   psa1_eof=alll(psa1_eof,lev,lat_CVDP,lon_CVDP)

   ;Load basin accumulation time series
   in = addfile("compositing_work/output/AISAccumulationTimeSeries.nc","r")
   AISAccum=in->time_series
   in := addfile("compositing_work/output/BasinAccumulationTimeSeries.nc","r")
   BasinAccum=in->time_series
   in := addfile("compositing_work/output/TimeSeriestime.nc","r")
   time=in->time_series
   BasinAccum!0="time"
   BasinAccum!1="basin"
   ds=dimsizes(BasinAccum)
   cc=escorc(AISAccum,BasinAccum(basin|:,time|:))
   
   ;Figure out which basins tend to contribute most to AIS-wide integrated accumulation anomaly
   BasinImportance=new(ds(1),"integer")
   BasinImportance(:)=0 
   ;get anomaly in total AIS snowfall
   AISAnom=dim_rmvmean(AISAccum)
   BasinAnoms=dim_rmvmean_n(BasinAccum,0)
   BasinAnoms!0="time"
   BasinAnoms!1="basin"   
   cc=escorc(AISAnom,BasinAnoms(basin|:,time|:))
   do t=1,ds(0)
      tm1=t-1
      AbsDif=abs(BasinAnoms(tm1,:)-AISAnom(tm1))
      ib=minind(AbsDif)
      BasinImportance(ib)=BasinImportance(ib)+1
   end do
   do b=1,ds(1)
     print("***Basin "+b+" Contribution to total AIS variability***")
     print("    #of years where basin leads in contribution to AIS variability: "+BasinImportance(b-1))
     print("    Correlation between basin variability and AIS variability: "+cc(b-1))
   end do

;;;Calculate basin-specific annual area-integrated precipitation correlations
   out:=calculate_interbasin_correlations(BasinAccum,nAISMask)
   Xcor_basin_by_basin_CESM=out[4]
   Rtest_basin_by_basin_CESM=out[5]
   print("***CESM: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***CESM: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***CESM: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***CESM: Average number of significant negative correlations per basin="+avg(out[3]))   

   out:=calculate_interbasin_correlations(Basin_integrated_TS_RACMO,nAISMask)
   Xcor_basin_by_basin_RACMO=out[4]
   Rtest_basin_by_basin_RACMO=out[5]
   print("***RACMO: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***RACMO: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***RACMO: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***RACMO: Average number of significant negative correlations per basin="+avg(out[3]))  

   out:=calculate_interbasin_correlations(Basin_integrated_TS_PAGES,nAISMask)
   Xcor_basin_by_basin_PAGES=out[4]
   Rtest_basin_by_basin_PAGES=out[5]
   print("***PAGES: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***PAGES: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***PAGES: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***PAGES: Average number of significant negative correlations per basin="+avg(out[3]))   

   Xcor_basin_by_basin_CESM_copy =Xcor_basin_by_basin_CESM
   Rtest_basin_by_basin_CESM_copy=Rtest_basin_by_basin_CESM       
   Xcor_basin_by_basin_RACMO_copy =Xcor_basin_by_basin_RACMO
   Rtest_basin_by_basin_RACMO_copy=Rtest_basin_by_basin_RACMO
   Xcor_basin_by_basin_PAGES_copy=Xcor_basin_by_basin_PAGES
   Rtest_basin_by_basin_PAGES_copy=Rtest_basin_by_basin_PAGES
   
   do b=1,toint(nAISMask)
      endm1=toint(nAISMask)-1
      bm1=b-1
      
      Xcor_basin_by_basin_CESM_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_CESM@_FillValue
      Rtest_basin_by_basin_CESM_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_CESM@_FillValue	     
      
      Xcor_basin_by_basin_RACMO_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_RACMO@_FillValue
      Rtest_basin_by_basin_RACMO_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_RACMO@_FillValue
      
      Xcor_basin_by_basin_PAGES_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_PAGES@_FillValue
      Rtest_basin_by_basin_PAGES_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_PAGES@_FillValue   
         		  
   end do
   wks = gsn_open_wks (wks_type,"figs/RACMO_vs_CESM_correlation_scatter")
   x=ndtooned(Xcor_basin_by_basin_CESM_copy)
   y=ndtooned(Xcor_basin_by_basin_RACMO_copy)
   test=regline(x,y)
   rc =  regline_stats(x,y) ; linear regression coef
   print(rc)
    
   res=True
   res@xyMarkLineModes     = (/"Markers"/)  ; choose which have markers
   res@xyMarkers           = 16                     ; choose type of marker 
   res@xyMarkerColor       = "red"                  ; Marker color
   res@tiYAxisString = "RACMO"
   res@tiXAxisString = "CESM"   
   
   plot=gsn_csm_xy (wks,x,y,res)
 
   print(test)
   exit   
   
   
   if (plot_matrix_of_CESM_interbasin_correlations.eq.1) then
       wks = gsn_open_wks (wks_type,"figs/CESM_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_CESM_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_CESM_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if
   
   if (plot_matrix_of_RACMO_interbasin_correlations.eq.1) then
       wks = gsn_open_wks (wks_type,"figs/RACMO_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_RACMO_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_RACMO_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if
   
   if (plot_matrix_of_PAGES_interbasin_correlations.eq.1) then
       wks = gsn_open_wks (wks_type,"figs/PAGES_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_PAGES_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_PAGES_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if 
   
   if (plot_basin_overview.eq.1) then
       wks = gsn_open_wks (wks_type,"figs/basin_overview")
       ;Plot topography
       Topog=alll(Topog,lev,lat,lon)
       Topog@_FillValue=default_fillvalue(typeof(Topog))
       ;Topog=where(Topog.gt.1,Topog,Topog@_FillValue)
       Topog_plot=plot_field(Topog,0,4000,False,"rainbow",True,-60,wks)
       Topog_plot=plot_basin_polygons(ZwallyBasinPolygonData,Topog_plot,wks,nAISMask,-9999)
       ;Plot basin numbers
       latAvg=new((/toint(nAISMask)/),"double")
       lonAvg=new((/toint(nAISMask)/),"double") 
       nCount=new((/toint(nAISMask)/),"integer")
       latAvg(:)=0.0
       lonAvg(:)=0.0
       nCount(:)=0     
       do i=1,dimsizes(lat)
          do j=1,dimsizes(lon)
	     im1=i-1
	     jm1=j-1
	     MaskVal=toint(AISMask(im1,jm1))-1
	     if (MaskVal.ge.0) then
		latAvg(MaskVal)=latAvg(MaskVal)+lat(im1)
		lonAvg(MaskVal)=lonAvg(MaskVal)+lon(jm1)
		nCount(MaskVal)=nCount(MaskVal)+1
	     end if
	  end do
       end do
       latAvg=latAvg/nCount
       lonAvg=lonAvg/nCount
       lonAvg=where(lonAvg.gt.180,lonAvg-360,lonAvg)
       lonAvg(1)=-10.;Basin 2
       latAvg(1)=-89.
       lonAvg(2)=0.;Basin 3
       latAvg(2)=-80.
       lonAvg(4)=0.;Basin 5
       latAvg(4)=-72.5
       latAvg(17)=-82.5;Basin 18
       latAvg(22)=-71.  ;Basin 23
       lonAvg(23)=-69.   ;Basin 24
       latAvg(23)=-71.   ;Basin 24       
       lonAvg(24)=-71.   ;Basin 25
       lonAvg(26)=-55   ;Basin 27                
       txres = True                     ; text mods desired
       txres@txFontHeightF = 0.015                     ; font smaller. default big
       do b=1,nAISMask
           bm1=b-1
           dum=gsn_add_text(wks,Topog_plot,tostring(b),lonAvg(bm1),latAvg(bm1),txres)
       end do
       
       draw(Topog_plot)
       frame(wks)           
   end if

   do b=1,nAISMask
      bm1=b-1
      
      ;;;Load basin-specific composite data
      
      out:=load_composites(SD,"aice",b,0.01,lev,lat,lon)
      SIC_mean_low=out[0]
      SIC_mean_high=out[1]
      SIC_diff=out[2] 
      SIC_alpha=out[3]
      
      out:=load_composites(SD,"PSL",b,1.,lev,lat,lon)
      SLP_mean_low=out[0]
      SLP_mean_high=out[1]
      SLP_diff=out[2]
      SLP_alpha=out[3]      
      
      out:=load_composites(SD,"QFLX",b,1.,lev,lat,lon)
      EVAP_mean_low=out[0]
      EVAP_mean_high=out[1]
      EVAP_diff=out[2]
      EVAP_alpha=out[3]
      
      out:=load_composites(SD,"SST",b,1.,lev,lat,lon)
      SST_mean_low=out[0]
      SST_mean_high=out[1]
      SST_diff=out[2]
      SST_alpha=out[3]    
                  
      out:=load_composites(SD,"Z3",b,1.,lev,lat,lon)
      Z3_mean_low=out[0]
      Z3_mean_high=out[1]
      Z3_diff=out[2]
      Z3_alpha=out[3]
   
      out:=load_composites(SD,"PRECIP",b, (year2sec/1000.),lev,lat,lon);mm/s -> m/yr 
      PRECIP_mean_low=out[0]
      PRECIP_mean_high=out[1]
      PRECIP_diff=out[2]
      PRECIP_alpha=out[3]      

      out:=load_composites(SD,"PS",b, 1,lev,lat,lon);mm/s -> m/yr 
      PS_mean_low=out[0]
      PS_mean_high=out[1]
      PS_diff=out[2]
      PS_alpha=out[3]  
      PS_mean_mean=(PS_mean_low+PS_mean_high)/2.

      out:=load_Hailong_composites(b,"low",lev,lat,lon)
      UQ_VERT_INTEGRAL_mean_low=out[0]
      VQ_VERT_INTEGRAL_mean_low=out[1]
      CONVERGENCE_mean_low=out[2]
      PE_mean_low=out[3]
      out:=load_Hailong_composites(b,"high",lev,lat,lon)
      UQ_VERT_INTEGRAL_mean_high=out[0]
      VQ_VERT_INTEGRAL_mean_high=out[1]
      CONVERGENCE_mean_high=out[2]
      PE_mean_high=out[3]      
      
      UQ_VERT_INTEGRAL_diff=(UQ_VERT_INTEGRAL_mean_high-UQ_VERT_INTEGRAL_mean_low)
      VQ_VERT_INTEGRAL_diff=(VQ_VERT_INTEGRAL_mean_high-VQ_VERT_INTEGRAL_mean_low)
      UQ_VERT_INTEGRAL_diff=alll(UQ_VERT_INTEGRAL_diff,lev,lat,lon)
      VQ_VERT_INTEGRAL_diff=alll(VQ_VERT_INTEGRAL_diff,lev,lat,lon)
      
      CONVERGENCE_diff=CONVERGENCE_mean_high-CONVERGENCE_mean_low
      PE_diff=PE_mean_high-PE_mean_low
      CONVERGENCE_diff=alll(CONVERGENCE_diff,lev,lat,lon)
      PE_diff=alll(PE_diff,lev,lat,lon)      
             
      ;Calculate z500
      Z500_diff =rm_single_dims (vinth2p (Z3_diff,hyam,hybm,(/500.0/),PS_mean_mean,1,P0/100.,1,True) )
      Z500_alpha=rm_single_dims (vinth2p (Z3_alpha,hyam,hybm,(/500.0/),PS_mean_mean,1,P0/100.,1,True))     
      
      rc=regline_stats(ndtooned(Z500_diff),ndtooned(SLP_diff))
      print("***Relationship between Z500 and SLP variability for basin "+b+"***")
      print("  regression coefficient: "+rc@r)

      BasinMask=new(dimsizes(AISMask),"integer")
      BasinMask@_FillValue=default_fillvalue(typeof(BasinMask))
      BasinMask=where(AISMask.eq.b,1,0);set basinmask to one where basin equals 1
      BasinMask=alll(BasinMask,lev,lat,lon)
      out:=make_correlation_mask(AISMask,nAISMask,Xcor_basin_by_basin_CESM,time,bm1)
      CorrMask_CESM=alll(out[0],lev,lat,lon)
      RtestMask_CESM=alll(out[1],lev,lat,lon)
      out:=make_correlation_mask(AISMask_RACMO,nAISMask,Xcor_basin_by_basin_RACMO,time_RACMO,bm1)
      CorrMask_RACMO=out[0]
      RtestMask_RACMO=out[1]
      out:=make_correlation_mask(AISMask,nAISMask,Xcor_basin_by_basin_PAGES,time_PAGES,bm1)
      CorrMask_PAGES=alll(out[0],lev,lat,lon)
      RtestMask_PAGES=alll(out[1],lev,lat,lon)    
      
       ;;;;Various basin-specific plots!;;;;

      if(plot_dMoisture_dZ500_transport.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_dZ500",sprinti("%0.2i", b)/)))
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,1,-30)
	 moisture_plot=plot_vector_scalar_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,Z500_diff,wks,0,-30)
	 significance_plot=plot_significance(Z500_alpha,wks)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay(correlation_plot,moisture_plot)
	 overlay (correlation_plot,significance_plot)
	 overlay (correlation_plot,basin_plot)	  
	 draw(correlation_plot)
	 frame(wks)	    
      end if

      if(plot_dz500_plus_basin_correlation.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dz500_transport_plus_basin_correlation",sprinti("%0.2i", b)/)))
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,1,-30)
	 correlation_plot=plot_basin_polygons(ZwallyBasinPolygonData,correlation_plot,wks,nAISMask,b)
	 ;range=max((/abs(min(Z500_diff)),abs(max(Z500_diff))/))
	 range=30.
	 Z500_plot=plot_field(Z500_diff,-range,range,True,(/"blue","black","red"/),False,0,wks)
	 significance_plot=plot_significance(Z500_alpha,wks)
	 overlay(correlation_plot,Z500_plot)
	 overlay (correlation_plot,significance_plot) 
	 draw(correlation_plot)
	 frame(wks)	    
      end if

      if(plot_dMoisture_plus_basin_correlation.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_plus_basin_correlation",sprinti("%0.2i", b)/)))
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,1,-60)
	 correlation_plot=plot_basin_polygons(ZwallyBasinPolygonData,correlation_plot,wks,nAISMask,b)
	 moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 significance_plot = plot_significance(RtestMask_CESM,wks)	 	 
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay(correlation_plot,moisture_plot)
	 overlay (correlation_plot,significance_plot)
	 draw(correlation_plot)
	 frame(wks)	    
      end if
      
      if(plot_dMoisture_plus_PE.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_dSIC_dPE",sprinti("%0.2i", b)/)))
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,-30)
	 ;correlation_plot=plot_basin_polygons(ZwallyBasinPolygonData,correlation_plot,wks,nAISMask,b)
	 dSIC_plot=plot_SIC(SIC_diff,wks,-55)
         PE_plot=plot_field(PE_diff,-10.,10.,True,(/"saddlebrown","transparent","darkgreen"/),False,0,wks)
	 ;moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 basin_plot=plot_basin(BasinMask,wks)
	 ;overlay (dSIC_plot,moisture_plot)
	 overlay (dSIC_plot,correlation_plot)	 
	 overlay (dSIC_plot,PE_plot)
	 overlay (dSIC_plot,basin_plot)	 	 
	 draw(dSIC_plot)
	 frame(wks)	 
      end if
      
      if(plot_dSIC_plus_dZ500.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC_plus_dZ500",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks,-45)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
	 dZ500_plot=plot_SLP(Z500_diff,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot) 
         overlay(dSIC_plot,dZ500_plot)	 
	 draw(dSIC_plot)
	 frame(wks)	      
      end if  
      
      if(plot_dSIC.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks,-45)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot)   
	 draw(dSIC_plot)
	 frame(wks)	      
      end if        
      
      if(plot_dMoisture_plus_SIC.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_plus_SIC",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks,-45)
         moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot)
         overlay(dSIC_plot,moisture_plot)	   
	 draw(dSIC_plot)
	 frame(wks)	 
      end if
      
      if (plot_evap_to_transport_ratio.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/dEVAP",sprinti("%0.2i", b)/)))
	 res = default_res   
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = -50 
	 res@tiMainString         = "EVAP Difference (high - low)"
	 res@cnFillOn             = False       
	 dEVAP_plot = gsn_csm_contour_map_polar(wks,EVAP_diff,res)  
	 draw (dEVAP_plot)
	 frame(wks)      
         ;TO BE CONTINUED!  See Burrows email, and Tsukernik 2013 for approach
      end if

      if (plot_spatial_correlation_CESM.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation",sprinti("%0.2i", b)/)))
	 correlation_plot  = plot_basin_correlation(CorrMask_CESM,wks,1,-60) 
	 significance_plot = plot_significance(RtestMask_CESM,wks)
	 basin_plot        = plot_basin(BasinMask,wks)
	 overlay (correlation_plot,significance_plot) 
	 overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if
      
      if (plot_spatial_correlation_RACMO.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation_RACMO",sprinti("%0.2i", b)/)))	 	 
	 correlation_plot  = plot_basin_correlation(CorrMask_RACMO,wks,0,0)	 
	 significance_plot = plot_significance(RtestMask_RACMO,wks)	 
	 ;basin_plot        = plot_basin(BasinMask_RACMO,wks)		  
	 overlay (correlation_plot,significance_plot) 
	 ;overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if      
      
      if (plot_spatial_correlation_PAGES.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation_PAGES",sprinti("%0.2i", b)/)))	 	 
	 correlation_plot  = plot_basin_correlation(CorrMask_PAGES,wks,1,-60)
	 correlation_plot  = plot_core_locations(correlation_plot,wks,PAGES_core_lat,PAGES_core_lon)
	 significance_plot = plot_significance(RtestMask_PAGES,wks)	 
	 basin_plot        = plot_basin(BasinMask,wks)		  
	 overlay (correlation_plot,significance_plot) 
	 overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if            
           
      if (plot_seasonal_cycle_difference.eq.1) then
      
         Monthly_integrated_cycle_low=calculate_integrated_timeseries(PRECIPm_low,LandArea,LandFrac*BasinMask)
	 Monthly_integrated_cycle_low=Monthly_integrated_cycle_low / 1.e9 ; integrated PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
         Monthly_integrated_cycle_high=calculate_integrated_timeseries(PRECIPm_high,LandArea,LandFrac*BasinMask)
	 Monthly_integrated_cycle_high=Monthly_integrated_cycle_high / 1.e9 ; integrated PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr

         dMonthly_integrated_cycle=Monthly_integrated_cycle_high-Monthly_integrated_cycle_low
	 
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Seasonal_cycle" ,sprinti("%0.2i", b)/)))
	 months=ispan(0,11,1)
	 res=default_res
	 res@xyLineThicknessF = 4      
	 res@trYMinF = 0.
	 res@trXMaxF = 11	 
	 res@xyLineColor      = "red"    
	 res@tmXBMode             = "Explicit"      ; label independently
	 res@tmXBValues           = months
         res@tmXBLabels           = (/"J","F","M","A","M","J",\
                                      "J","A","S","O","N","D"/)
	 res@tiYAxisString = "Integrated accumulation (Gt/yr)"
	 res@tiXAxisString = "Month"
         plot0=gsn_csm_xy (wks,months,Monthly_integrated_cycle_high,res) 
	 res@xyLineColor      = "blue" 	 
	 plot1=gsn_csm_xy (wks,months,Monthly_integrated_cycle_low,res)
	 res@xyLineColor      = "purple" 
	 res@xyLineThicknessF = 12	 	 
	 plot2=gsn_csm_xy (wks,months,dMonthly_integrated_cycle,res)	 
	 
	 overlay(plot0,plot1)
	 overlay(plot0,plot2)	 
	 draw(plot0)
	 frame(wks)

      end if

   end do
 
end
