load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=False
plot_Annual_climatology_difference=False
plot_DJF_climatology_difference=True

year2sec=31557600.
LandArea = get_LandArea() 
LandFrac = get_LandFrac()
AISMask = get_AIS_mask()
AISMaskBinary=where(AISMask.gt.0,1,0)
nAISMask = toint(max(AISMask))
;Load CAM grid data

in:=addfile("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/CAM_base_data.nc","r")
lev=in->lev
lat=in->lat
lon=in->lon

ZwallyBasinPolygonData=readAsciiTable("/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/AIS_basin_masks/Ant_Full_DrainageSystem_Polygons.txt",3,"float",7)

EnsNums=(/1,10,12,15,18,20,22,28/)
nEns=dimsizes(EnsNums)

YearStart=1986
YearEnd=2005
MonStart=YearStart*12
MonEnd=YearEnd*12+12
nYear=YearEnd-YearStart+1
nMon=MonEnd-MonStart

Years=fspan(1,nYear,nYear)+YearEnd-nYear
OutputFileName="/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/Ozone_analysis_spatial"+tostring(YearStart)+"_"+tostring(YearEnd)+".nc"

iO3=0
iNO3=1

iDJF=0
iAnn=1
PeriodName=(/"DJF", "Annual"/)

if (LoadFreshData) then;
   PRECIP=load_no_ozone_lens_data(MonStart,MonEnd,EnsNums,2,"PRECC")+\
          load_no_ozone_lens_data(MonStart,MonEnd,EnsNums,2,"PRECL")       
   Z3=load_no_ozone_lens_data(MonStart,MonEnd,EnsNums,3,"Z3")
   PSFC=load_no_ozone_lens_data(MonStart,MonEnd,EnsNums,2,"PS")

   print("Loading hyam,hybm")
   fin="/glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/PS/b.e11.B20LE_fixedO3_001.cam.h0.PS.195501-200512.nc"
   in=addfile(fin,"r")
   hyam = in->hyam
   hybm = in->hybm

   P0mb=(100000.)/100. ; Apparently this is the standard CESM reference pressure.  convert to hPa/mb
   
   PRECIP_Annually_Resolved=new((/nYear,192,288,nEns,2,2/),double)  
   PRECIP_Annually_Resolved!0="time"
   PRECIP_Annually_Resolved!1="lat" 
   PRECIP_Annually_Resolved!2="lon"
   PRECIP_Annually_Resolved!3="ensemble_number"
   PRECIP_Annually_Resolved!4="experiment" 
   PRECIP_Annually_Resolved!5="time_period"   
   
   print("PRECIP averaging...")
   do en=0,nEns-1
      do e=0,1;O3,NO3
	    ;Do Annual and DJF averaging here.
	    tmp:=(/PRECIP(:,:,:,en,e)/)
	    tmp!0="time" ;just so month_to_season works
	    tmp!1="lat" ;just so month_to_season works.  This is actually eEns.   
	    tmp!2="lon" ;just so month_to_season works.  This is actually exper.
	    PRECIP_Annually_Resolved(:,:,:,en,e,iDJF)=(/month_to_season(tmp,"DJF")/)
	    PRECIP_Annually_Resolved(:,:,:,en,e,iAnn)=(/month_to_annual(tmp,1) /)
      end do
   end do

   Z500_Annually_Resolved=new((/nYear,192,288,nEns,2,1/),double)  ;Just do DJF
   Z500_Annually_Resolved!0="time"
   Z500_Annually_Resolved!1="lat" 
   Z500_Annually_Resolved!2="lon"
   Z500_Annually_Resolved!3="ensemble_number"
   Z500_Annually_Resolved!4="experiment" 
   Z500_Annually_Resolved!5="time_period"

   print("Z500 averaging...")     
   do en=0,nEns-1
      do e=0,1;O3,NO3
         Z500_Annually_Resolved(:,:,:,en,e,iDJF) = (/month_to_season( vinth2p (Z3(:,:,:,:,en,e),hyam,hybm,(/ 500.0/),PSFC(:,:,:,en,e),1,P0mb,1,True)  ,"DJF" )/)
      end do
   end do
   
   print("Saving data...")
   system("rm -f "+OutputFileName)
   ncdf=addfile(OutputFileName ,"c")
   ncdf->PRECIP_Annually_Resolved=PRECIP_Annually_Resolved
   ncdf->Z500_Annually_Resolved  =Z500_Annually_Resolved
   
else

   print("Loading data...")
   ncdf=addfile(OutputFileName,"r")
   PRECIP_Annually_Resolved=ncdf->PRECIP_Annually_Resolved
   Z500_Annually_Resolved  =ncdf->Z500_Annually_Resolved
   
end if

PRECIP_Annually_Resolved_Climo=dim_avg_n   (PRECIP_Annually_Resolved,0) ;Time average
PRECIP_Annually_Resolved_Climo_Ensemble_Mean=dim_avg_n   (PRECIP_Annually_Resolved_Climo,2) ;Ensemble average
dPRECIP_Annually_Resolved_Climo_Ensemble=(PRECIP_Annually_Resolved_Climo_Ensemble_Mean(:,:,iO3,:) -\
                                          PRECIP_Annually_Resolved_Climo_Ensemble_Mean(:,:,iNO3,:)) /\
                                          PRECIP_Annually_Resolved_Climo_Ensemble_Mean(:,:,iO3,:);Relative change

dPRECIP_Annually_Resolved_Climo_Ensemble_Pval=new((/192,288,2/),double)
;Can't figure out how to do this as a single call, so looping it...
do lati=0,80 ; don't do whole globe, to save processing time.
   do loni=0,287
      do p=0,1
	  aveX=avg(PRECIP_Annually_Resolved(:,lati,loni,:,iO3,p))
	  varX=variance(PRECIP_Annually_Resolved(:,lati,loni,:,iO3,p))
	  sX=nYear*nEns
	  aveY=avg(PRECIP_Annually_Resolved(:,lati,loni,:,iNO3,p))
	  varY=variance(PRECIP_Annually_Resolved(:,lati,loni,:,iNO3,p))
	  sY=nYear*nEns
	  dPRECIP_Annually_Resolved_Climo_Ensemble_Pval(lati,loni,p)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)         
      end do      
   end do
end do

Z500_Annually_Resolved_Climo=dim_avg_n (Z500_Annually_Resolved,0) ;Time average
Z500_Annually_Resolved_Climo_Ensemble_Mean=dim_avg_n (Z500_Annually_Resolved_Climo,2) ;Time average
dZ500_Annually_Resolved_Climo_Ensemble=Z500_Annually_Resolved_Climo_Ensemble_Mean(:,:,iO3,:) -\
                                       Z500_Annually_Resolved_Climo_Ensemble_Mean(:,:,iNO3,:)

;TO DO:
;-get basin-scale Gt/yr changes
;-get dot locations
;-plot basin-scale Gt/yr changes
;-integrate Jan figure comments

;;;Plotting below

wks_type=get_default_wks(); function in general_utilities.ncl

if (plot_DJF_climatology_difference) then
   latbound=-58.
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_REL_DJF")
   PRECIP_plot = plot_field(alll(dPRECIP_Annually_Resolved_Climo_Ensemble(:,:,iDJF),lev,lat,lon),-.2,.2,False,"MPL_BrBG",True,latbound,wks)
   significance_plot=plot_significance(alll(dPRECIP_Annually_Resolved_Climo_Ensemble_Pval(:,:,iDJF),lev,lat,lon),wks)   
   PRECIP_plot = plot_basin_polygons(ZwallyBasinPolygonData,PRECIP_plot,wks,nAISMask,-9999)
   Z500_plot = plot_field(alll(dZ500_Annually_Resolved_Climo_Ensemble(:,:,iDJF),lev,lat,lon),-25.,25.,True,"MPL_bwr",False,latbound,wks)
   overlay(PRECIP_plot,significance_plot)
   overlay(PRECIP_plot,Z500_plot)
   draw(PRECIP_plot)
   frame(wks)    
end if

if (plot_Annual_climatology_difference) then
   latbound=-58.
   wks = gsn_open_wks (wks_type,"figs/dPRECIP_REL_Ann")
   PRECIP_plot = plot_field(alll(dPRECIP_Annually_Resolved_Climo_Ensemble(:,:,iAnn),lev,lat,lon),-.08,.08,False,"MPL_BrBG",True,latbound,wks)
   PRECIP_plot = plot_basin_polygons(ZwallyBasinPolygonData,PRECIP_plot,wks,nAISMask,-9999)
   draw(PRECIP_plot)
   frame(wks)    
end if
