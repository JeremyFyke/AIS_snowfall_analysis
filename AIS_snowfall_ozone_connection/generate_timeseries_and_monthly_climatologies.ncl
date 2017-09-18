load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=False

plot_climatological_annual_cycle=False
plot_DJF_Annual_timeseries_difference=True

EnsNums=(/1,10,12,15,18,20,22,28/)
nEns=dimsizes(EnsNums)

YearStart=1955
YearEnd=2005
MonStart=YearStart*12
MonEnd=YearEnd*12+12
nYear=YearEnd-YearStart+1
nMon=MonEnd-MonStart
Years=fspan(1,nYear,nYear)+YearEnd-nYear

OutputFileName="/glade/p/cesm/liwg/AIS_snowfall_analysis/input_data/Ozone_analysis_"+tostring(YearStart)+"_"+tostring(YearEnd)+".nc"

iO3=0
iNO3=1

iDJF=0
iAnn=1
PeriodName=(/"DJF", "Annual"/)

LandArea = get_LandArea() 
LandFrac = get_LandFrac()
AISMask = get_AIS_mask()
AISMaskBinary=where(AISMask.gt.0.,1,0)
nAISMask = toint(max(AISMask))

TotalArea=sum(LandArea*LandFrac*AISMaskBinary)/1.e12
print("TotalArea according to CESM masks (including shelves)="+TotalArea+" million km^2")

AISMaskGrounded=get_AIS_mask_grounded()
AISMaskGroundedBinary=where(AISMaskGrounded.gt.0.,1.,0.)

TotalAreaGrounded=sum(AISMaskGrounded(27,:,:)+AISMaskGrounded(28,:,:))/1.e12
print("TotalAreaGrounded according to Brooke masks="+TotalAreaGrounded+" million km^2")


year2sec=31557600.
m3ps_to_GtpYr=1. * year2sec *1000. / 1.e12;

if (LoadFreshData) then; Regenerate timeseries from raw ensemble data

   PRECIP_Integrated_TS                = new((/nMon,nEns,2/),double)
   PRECIP_Integrated_TS_Grounded       = new((/nMon,nEns,2/),double)   
   PRECIP_BasinIntegrated_TS           = new((/nMon,nEns,nAISMask,2/),double)
   PRECIP_BasinIntegrated_TS_Grounded  = new((/nMon,nEns,nAISMask,2/),double)  
   
   ;Load ensemble monthly precipitation data, for all members of dynamic and fixed ozone ensembles.
   ;Function load_no_ozone_lens_composites is in ../common_ncl_utilies/general_utilities.ncl
   ;PRECIP dimension: time, lat, lon, ensemble #, experiment (dynamic ozone, no ozone)
   PRECIP=load_no_ozone_lens_composites(MonStart,MonEnd,EnsNums,"PRECC")+\
          load_no_ozone_lens_composites(MonStart,MonEnd,EnsNums,"PRECL")

   ;dPRECIP_CESM_basin_pval = new((/nAISMask/),double)
   ;dPRECIP_CESM_basin = new((/nAISMask/),double)
   ;dPRECIP_CESM_mask=new(dimsizes(AISMask),"double")
   ;dPRECIP_CESM_pval_mask=new(dimsizes(AISMask),"double")

;      aveX=avg(Basin_integrated_dSNOWFALL_O3(:,:,b))
;      varX=variance(Basin_integrated_dSNOWFALL_O3(:,:,b))
;      aveY=avg(Basin_integrated_dSNOWFALL_NO3(:,:,b))
;      varY=variance(Basin_integrated_dSNOWFALL_NO3(:,:,b))
;      dPRECIP_CESM_basin(b)=(aveX-aveY)/aveX ;Relative change
;      dPRECIP_CESM_basin_pval(b)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)  
;
;      dPRECIP_CESM_mask=where(AISMask.eq.b,dPRECIP_CESM_basin(b),dPRECIP_CESM_mask)
;      dPRECIP_CESM_pval_mask=where(AISMask.eq.b,dPRECIP_CESM_basin_pval(b),dPRECIP_CESM_pval_mask)

   do e=0,nEns-1
      print("Calculating integrated time series for "+tostring(e)+"...")
      do exper=0,1
	 ;Calculate integrated time series.
	 Mask=new((/192,288/),double)
	 Mask(:,:)=1.0
	 ;Function calculate_integrated_timeseries is in ../common_ncl_utilies/general_utilities.ncl
	 GroundedArea=dim_sum_n(AISMaskGrounded(27:28,:,:),0)
	 PRECIP_Integrated_TS(:,e,exper)         =calculate_integrated_timeseries(PRECIP(:,:,:,e,exper),LandArea,LandFrac*AISMaskBinary)
	 PRECIP_Integrated_TS_Grounded(:,e,exper)=calculate_integrated_timeseries(PRECIP(:,:,:,e,exper),GroundedArea,Mask)

	 do b = 0,nAISMask-1
	    BasinMask=where(AISMask.eq.b+1,1.,0.)
	    PRECIP_BasinIntegrated_TS(:,e,b,exper)=calculate_integrated_timeseries(PRECIP(:,:,:,e,exper),LandArea,BasinMask*LandFrac)
	    PRECIP_BasinIntegrated_TS_Grounded(:,e,b,exper)=calculate_integrated_timeseries(PRECIP(:,:,:,e,exper),AISMaskGrounded(b,:,:),Mask)
	 end do
      end do
   end do
   
   print("Saving data...")
   system("rm -f "+OutputFileName)
   ncdf=addfile(OutputFileName ,"c")
   ncdf->PRECIP_Integrated_TS=PRECIP_Integrated_TS
   ncdf->PRECIP_BasinIntegrated_TS=PRECIP_BasinIntegrated_TS   
   ncdf->PRECIP_Integrated_TS_Grounded=PRECIP_Integrated_TS_Grounded
   ncdf->PRECIP_BasinIntegrated_TS_Grounded=PRECIP_BasinIntegrated_TS_Grounded  
   
else;Load previously-generated time series data

   print("Loading data...")
   ncdf=addfile(OutputFileName,"r")
   PRECIP_Integrated_TS=ncdf->PRECIP_Integrated_TS
   PRECIP_BasinIntegrated_TS=ncdf->PRECIP_BasinIntegrated_TS
   PRECIP_Integrated_TS_Grounded=ncdf->PRECIP_Integrated_TS_Grounded
   PRECIP_BasinIntegrated_TS_Grounded=ncdf->PRECIP_BasinIntegrated_TS_Grounded  
   
end if

;;;CHOOSE EITHER GROUNDED OR WHOLE ICE SHEET time series to analyze
;TS_Array=PRECIP_Integrated_TS
;TS_Array_Basin=PRECIP_BasinIntegrated_TS
nEns=4
TS_Array=PRECIP_Integrated_TS_Grounded(:,0:nEns-1,:)
TS_Array_Basin=PRECIP_BasinIntegrated_TS_Grounded(:,0:nEns-1,:,:)

TS_Array!0="time" ;just so month_to_season works
TS_Array!1="lat" ;just so month_to_season works.  This is actually eEns.   
TS_Array!2="lon" ;just so month_to_season works.  This is actually exper.

TS_Array_Basin!0="time" ;just so month_to_season works
TS_Array_Basin!2="lev" ;just so month_to_season works.  This is actually eEns.
TS_Array_Basin!2="lat" ;just so month_to_season works.  This is actually basin.   
TS_Array_Basin!3="lon" ;just so month_to_season works.  This is actually exper.

PRECIP_Integrated_TS_Annually_Resolved=new((/nYear,nEns,2,2/),double); time, nEns, exper, season
PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iDJF)=month_to_season(TS_Array,"DJF")
PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iDJF)=PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iDJF)*m3ps_to_GtpYr/4.
PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iAnn)=(/month_to_annual(TS_Array,1)/) ;brackets omit metadata
PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iAnn)=PRECIP_Integrated_TS_Annually_Resolved(:,:,:,iAnn)*m3ps_to_GtpYr

PRECIP_Integrated_TS_Annually_Resolved_Basin=new((/nYear,nEns,nAISMask,2,2/),double); time, nEns, basin, exper, season
PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iDJF)=month_to_season(TS_Array_Basin,"DJF")
PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iDJF)=PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iDJF)*m3ps_to_GtpYr/4.
PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iAnn)=(/month_to_annual(TS_Array_Basin,1)/) ;brackets omit metadata
PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iAnn)=PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,:,:,iAnn)*m3ps_to_GtpYr

SmoothPer=11
SmoothPerMid=SmoothPer/2 ;Note integer division, which truncates decimal remainder.

PRECIP_Integrated_TS_Annually_Resolved_Mean         =dim_avg_n(PRECIP_Integrated_TS_Annually_Resolved,1);ensemble average, dimenions: time,exper,season
PRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed=runave_n(PRECIP_Integrated_TS_Annually_Resolved_Mean,SmoothPer,0,0); times-smoothed, same dimensions as above

;Calculated differences in raw/smoothed timeseries, for both DJF and annual
dPRECIP_Integrated_TS_Annually_Resolved_Mean         =PRECIP_Integrated_TS_Annually_Resolved_Mean(:,iO3,:)-\
                                                      PRECIP_Integrated_TS_Annually_Resolved_Mean(:,iNO3,:); time, season
dPRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed=PRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed(:,iO3,:)-\
                                                      PRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed(:,iNO3,:)

dPRECIP_Integrated_TS_Annually_Resolved_pval=new((/nYear,2/),double); time, season
tEmerged=new((/2/),integer)

do s=0,1;iDJF,iANN
  ;Calculate running significance and smoothed time series (used for time series plots)
  tEmergedFound=False
  do t=SmoothPerMid,nYear-SmoothPerMid-1
     ts=t-SmoothPerMid
     te=t+SmoothPerMid
     aveX=     avg(PRECIP_Integrated_TS_Annually_Resolved(ts:te,:,iO3,s))
     varX=variance(PRECIP_Integrated_TS_Annually_Resolved(ts:te,:,iO3,s))
     aveY=     avg(PRECIP_Integrated_TS_Annually_Resolved(ts:te,:,iNO3,s))
     varY=variance(PRECIP_Integrated_TS_Annually_Resolved(ts:te,:,iNO3,s))
     sX=SmoothPer*nEns
     sY=SmoothPer*nEns
     dPRECIP_Integrated_TS_Annually_Resolved_pval(t,s)=ttest(aveX,varX,sX,aveY,varY,sY,False,False)
     ;Identify first time of significant emergence in smoothed time series
     if (dPRECIP_Integrated_TS_Annually_Resolved_pval(t,s).lt.0.01 .and. .not. tEmergedFound) then
	tEmerged(s)=t
	tEmergedFound=True
     end if
  end do
  if (tEmergedFound)
     print("Year emerged="+tostring(tEmerged(s)))
  else
     print("No emergence of significant difference found, stopping so you can figure things out.")
     ;exit
  end if
end do

;Print out scalar diagnostics to terminal
sX=nYear*nEns
sY=sX
do s=0,1;iDJF,iANN
  print("*********************************")
  aveX=avg     (PRECIP_Integrated_TS_Annually_Resolved(:,:,iO3,s))
  varX=variance(PRECIP_Integrated_TS_Annually_Resolved(:,:,iO3,s))
  stderrX=stddev  (PRECIP_Integrated_TS_Annually_Resolved(:,:,iO3,s))/sqrt(sX)
  aveY=avg     (PRECIP_Integrated_TS_Annually_Resolved(:,:,iNO3,s))
  varY=variance(PRECIP_Integrated_TS_Annually_Resolved(:,:,iNO3,s))
  stderrY=stddev  (PRECIP_Integrated_TS_Annually_Resolved(:,:,iNO3,s))/sqrt(sY)
  pval=ttest(aveX,varX,sX,aveY,varY,sY,False,False)
  print( "Average dynamic ozone "+PeriodName(s)+" mass flux="+tostring(aveX))
  print( "Average fixed ozone "+PeriodName(s)+" mass flux="  +tostring(aveY))
  print( "Mean "+PeriodName(s)+" ozone-forced snowfall change=" + tostring((aveX-aveY)))
  print( "Pval of "+PeriodName(s)+" ozone-forced snowfall change=" +tostring(pval))
  stderr=sqrt(stderrX^2 + stderrY^2)
  print( "Standard error of the mean of "+PeriodName(s)+" ozone-forced snowfall change=" + stderr)
  do b = 0,nAISMask-1
     ;print("   ********Basin"+tostring(b+1)+" values*************")  
     aveX=avg     (PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iO3,s))
     varX=variance(PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iO3,s))
     stderrX=stddev  (PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iO3,s))/sqrt(sX)
     aveY=avg     (PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iNO3,s))
     varY=variance(PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iNO3,s))
     stderrY=stddev  (PRECIP_Integrated_TS_Annually_Resolved_Basin(:,:,b,iNO3,s))/sqrt(sY)
     pval=ttest(aveX,varX,sX,aveY,varY,sY,False,False)
     ;print( "   Average dynamic ozone "+PeriodName(s)+" mass flux="+tostring(aveX))
     ;print( "   Average fixed ozone "+PeriodName(s)+" mass flux="  +tostring(aveY))
     ;print( "   Mean "+PeriodName(s)+" ozone-forced snowfall change=" + tostring((aveX-aveY)))
     ;print( "   Pval of "+PeriodName(s)+" ozone-forced snowfall change=" +tostring(pval))
     ;stderr=sqrt(stderrX^2 + stderrY^2)
     ;print( "   Standard error of the mean of "+PeriodName(s)+" ozone-forced snowfall change=" + stderr)   
  end do

end do
print("*********************************")

;Generate integrated seasonal cycles
PRECIP_Integrated_monthlyClimo=clmMonTLL(PRECIP_Integrated_TS); month, eEns, exper
PRECIP_Integrated_monthlyClimo=PRECIP_Integrated_monthlyClimo*m3ps_to_GtpYr/12. ; -> Gt/mon
PRECIP_Integrated_monthlyClimo=array_append_record(PRECIP_Integrated_monthlyClimo(6:11,:,:),PRECIP_Integrated_monthlyClimo(0:5,:,:),0) ; shift annual cycle
PRECIP_Integrated_monthlyClimo_Mean=dim_avg_n(PRECIP_Integrated_monthlyClimo,1); Month, exper
PRECIP_Integrated_monthlyClimo_STD=dim_stddev_n(PRECIP_Integrated_monthlyClimo,1)

;Calculate seasonal ranges
O3Range =max(PRECIP_Integrated_monthlyClimo_Mean(:,iO3))-min(PRECIP_Integrated_monthlyClimo_Mean(:,iO3))
NO3Range=max(PRECIP_Integrated_monthlyClimo_Mean(:,iNO3))-min(PRECIP_Integrated_monthlyClimo_Mean(:,iNO3))
RangeChange=(NO3Range-O3Range)/NO3Range
print("Change in seasonal cycle range due to ozone depletion)="+tostring(RangeChange))

;Load ozone timeseries
;ozoneFile="/glade/scratch/jfyke/ozone_data/ozone.nc"
;f=addfile(ozoneFile ,"r")
;ozoneTimeseries=f->

                               ;;;Plotting below;;;
			       
			       
			       
			       

wks_type=get_default_wks(); function in general_utilities.ncl

if (plot_climatological_annual_cycle) then
   months=ispan(0,11,1)
   wks = gsn_open_wks (wks_type,"figs/Monthly_climatology");Plot timeseries and ensemble spread of integrated precip
   miny=1300./12.
   maxy=3100./12.
   res = True                      ; plot mods desired   
   res@gsnDraw          = False
   res@gsnFrame         = False
   res@xyLineThicknesses = (/  10.0,  10.0 /) 
   res@xyLineColors      = (/"blue","red"/)
   res@xyDashPattern      = (/0,0/)   
   res@tmXBValues        = months
   res@tmXBLabels        = (/"J","A","S","O","N","D",\
                             "J","F","M","A","M","J"/)
   res@tmXBMode             = "Explicit"      ; label independently
   res@tiYAxisString = "Integrated accumulation (Gt/month)"
   res@tiXAxisString = "Month"
   res@trXMaxF = 11
   res@trYMinF = miny
   res@trYMaxF = maxy      
   tmp:=new((/2,12/),double)
   tmp(0,:)=PRECIP_Integrated_monthlyClimo_Mean(:,iO3)
   tmp(1,:)=PRECIP_Integrated_monthlyClimo_Mean(:,iNO3)
   base_plot = gsn_csm_xy(wks,months,tmp,res)
   
   ;;Make shaded ranges
   pgres                  = True
   pgres@gsFillOpacityF        = 0.2  
   xp:=new((/24/),float)
   yp:=new((/24/),double)
   xp(0:11)=months
   xp(12:23)=months(::-1)
   ;LENS
      yp(0:11)=PRECIP_Integrated_monthlyClimo_Mean(:,iO3)-PRECIP_Integrated_monthlyClimo_STD(:,iO3)
      yp(12:23)=PRECIP_Integrated_monthlyClimo_Mean(::-1,iO3)+PRECIP_Integrated_monthlyClimo_STD(::-1,iO3)
      pgres@gsFillColor	  = "Blue"
      dummy1 = gsn_add_polygon (wks,base_plot,xp,yp,pgres)
   ;Fixed-ozone      
      yp(0:11)=PRECIP_Integrated_monthlyClimo_Mean(:,iNO3)-PRECIP_Integrated_monthlyClimo_STD(:,iNO3)
      yp(12:23)=PRECIP_Integrated_monthlyClimo_Mean(::-1,iNO3)+PRECIP_Integrated_monthlyClimo_STD(::-1,iNO3)
      pgres@gsFillColor	  = "Red"
      dummy2 = gsn_add_polygon (wks,base_plot,xp,yp,pgres)
   ;Legend
      lineres=True
      lineres@lgLineThicknesses=(/20.0, 20.0/)
      lineres@lgLineColors=(/"blue","red"/)   
      gres = True
      gres@YPosPosition = 95.    ; expressed as %, 0->100, sets position of top border of legend 
      gres@XPosPosition = 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
      textres = True
      textres@lgLabels = (/"Dynamic ozone ensemble","Fixed ozone ensemble"/) 
      base_plot = simple_legend(wks,base_plot,gres,lineres,textres) 
  
   draw(base_plot)
   frame(wks)
end if

if (plot_DJF_Annual_timeseries_difference) then
   ;plot difference in integrated AIS accumulation due to ozone, overlay against estimated loss trend in the case of annual
   
   fname=(/"dPRECIP_DJF_O_TS","dPRECIP_Ann_O_TS"/)
   do s=0,1
      wks = gsn_open_wks (wks_type,"figs/"+fname(s)+"_"+tostring(nEns));Plot timeseries and ensemble spread of integrated precip
      res := True
      res@trXMinF=min(Years)
      res@trXMaxF=max(Years)
      if (s.eq.iAnn) then
         res@trYMinF=-100.
         res@trYMaxF=180.         
      else
	 res@trYMinF=-50.
	 res@trYMaxF=100.
      end if
      res@gsnDraw      = False                        ; don't draw yet
      res@gsnFrame     = False                        ; don't advance frame yet
      res@xyLineThicknesses := (/20.0/)          ; make second line thicker
      res@xyLineColors      := (/"grey"/)          ; change line color
      res@xyDashPattern = 2
      res@tiXAxisString =  "Year"
      res@tiYAxisString =  "Gt/yr"
      base_plot = gsn_csm_xy(wks,Years-0.5,dPRECIP_Integrated_TS_Annually_Resolved_Mean(:,s),res)

      res@xyDashPattern := 0
      res@xyLineColors      := (/"black"/)          ; change line color
      running_average_plot = gsn_csm_xy(wks,Years-0.5,dPRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed(:,s),res)

      res@xyLineThicknesses := (/60.0/)
      ts:=tEmerged(iDJF)
      te:=nYear-1
      ;running_average_plot_significant_diff = gsn_csm_xy(wks,Years(ts:te)-0.5,dPRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed(ts:te,s),res)
      
      lineres := True
      textres := True
      gres := True
      gres@YPosPosition := 95.    ; expressed as %, 0->100, sets position of top border of legend 
      gres@XPosPosition := 5      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
      te:=nYear-1 
      if (s.eq.iDJF) 	 
         ts:=tEmerged(iDJF)
	 ;Construct legend
	 lineres@lgLineThicknesses:=(/20.0, 20.0, 60.0/)
	 lineres@lgLineColors:=(/"grey","black","black"/)
	 lineres@lgDashPatterns:=(/2,0,0/)
	 textres@lgLabels := (/"DJF ensemble difference","11-year running mean ensemble difference","99% Stat. sig. 11-yr running mean ensemble difference"/) 
      else
         res@xyLineThicknesses := (/60.0/)
	 ts:=tEmerged(iAnn)
	 ;running_average_plot_significant_diff = gsn_csm_xy(wks,Years(ts:te)-0.5,dPRECIP_Integrated_TS_Annually_Resolved_Mean_Smoothed(ts:te,s),res)         
	 ;Add in observational mass loss estimate graphic
         ;Hard-code simple estimate of recent AIS mass loss (Gt/yr, Forsberg et al., 2017, 10.1007/s10712-016-9398-7) for plot comparisons
	 ForsbergObsStats=new((/3,4/),double)
	 ForsbergObsStats(0,:)=fspan(95.,95.,4) ;mean
	 ForsbergObsStats(1,:)=fspan(95.-50.,95.-50.,4) ;min
	 ForsbergObsStats(2,:)=fspan(95.+50.,95.+50.,4) ;max
	 ForsbergTime=fspan(2002.,2005,4)
         res@xyLineThicknesses := (/60.0,  40.0, 40.0/)          ; make second line thicker
	 res@xyDashPattern     := 0
	 res@xyLineColors      := (/"orange","orange","orange"/)          ; change line color
	 obs_plot = gsn_csm_xy(wks,ForsbergTime-0.5,ForsbergObsStats,res)
	 overlay(base_plot,obs_plot)
	 ;Construct legend
	 lineres@lgLineThicknesses=(/20.0, 20.0, 60.0, 60.0/)
	 lineres@lgLineColors=(/"grey","black","black","orange"/)
	 lineres@lgDashPatterns=(/2,0,0,0/)
	 textres@lgLabels = (/"Annual ensemble difference","11-year running mean ensemble difference","99% Stat. sig. 11-yr running mean ensemble difference","2002-2015 AIS mass loss estimate range"/) 
      end if
      te:=nYear-1
      res@xyLineColors:="black"
      base_plot = simple_legend(wks,base_plot,gres,lineres,textres)
      overlay(base_plot,running_average_plot)
      ;overlay(base_plot,running_average_plot_significant_diff)
      draw(base_plot)
      frame(wks)
   
   end do 
    
end if




