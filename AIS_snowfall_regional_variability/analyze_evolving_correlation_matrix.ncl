load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/composite_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/general_utilities.ncl"
load "/glade/p/cesm/liwg/AIS_snowfall_analysis/common_ncl_utilies/plotting_utilities.ncl"

LoadFreshData=False

OutputFileName="EvolvingCorrMats.nc"
SmoothPer=31
SmoothPerMid=SmoothPer/2 ;Note integer division, which truncates decimal remainder.

AISMask = get_AIS_mask()
nAISMask = max(AISMask)

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

if LoadFreshData

   Xcor=new((/ds(0)-SmoothPer+1,ds(1),ds(1)/),double)
   R=new((/ds(0)-SmoothPer+1,ds(1),ds(1)/),double)
   ;Call plot_basin_correlation_matrix for 30-year-windowed segments of 1800-year record.
   do t=SmoothPerMid,ds(0)-SmoothPerMid-1 
      tm=t-SmoothPerMid
      print(tm-1+1)
      out:=calculate_interbasin_correlations(BasinAccum(t-15:t+15,:),nAISMask)
      Xcor(tm,:,:)=out[4]
      R(tm,:,:)=out[5]     
   end do
  
   print("Saving data...")
   system("rm -f "+OutputFileName)
   ncdf=addfile(OutputFileName ,"c")
   ncdf->Xcor=Xcor
   ncdf->R=R
   
else 

   print("Loading data...")
   ncdf=addfile(OutputFileName,"r") 
   Xcor=ncdf->Xcor
   R   =ncdf->R
   
end if

ds:=dimsizes(Xcor)
Years=fspan(1,ds(0),ds(0))
res = True
res@trXMinF = min(Years)
res@trXMaxF = max(Years)
res@tiXAxisString = "Year"
res@tiYAxisString = tostring(SmoothPer)+"-year windowed correlation"
res@gsnYRefLine           = 0.0
polyres                  = True
polyres@gsLineThicknessF = 3.0
res@gsnAboveYRefLineColor = "red"              ; above ref line fill red
res@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue  

;Minimum correlation evolution
inds=rm_single_dims(ind_resolve(minind(ndtooned(Xcor)),ds))
TS=Xcor(:,inds(1),inds(2))
wks  = gsn_open_wks("png","figs/EvolvingMinCorrelation")
res@tiYAxisString = tostring(SmoothPer)+"-year correlation, basins "+tostring(inds(1)+1)+"/"+tostring(inds(2)+1)
plot=gsn_csm_xy(wks,Years,TS,res)
gsn_polyline(wks,plot,Years,(/TS/),polyres)

;Max correlation evolution
ds:=dimsizes(Xcor)
Xcor_copy=where(Xcor.gt..9999,0.0,Xcor)
inds=rm_single_dims(ind_resolve(maxind(ndtooned(Xcor_copy)),ds))
TS=Xcor(:,inds(2),inds(1))
wks  = gsn_open_wks("png","figs/EvolvingMaxCorrelation")
res@tiYAxisString = tostring(SmoothPer)+"-year correlation, basins "+tostring(inds(2)+1)+"/"+tostring(inds(1)+1)
plot=gsn_csm_xy(wks,Years,TS,res)
gsn_polyline(wks,plot,Years,(/TS/),polyres)

;Mixed correlation evolution
;Find correlation time series with both + and - values.

maxrange=todouble(0.0)
do b1 = 0, ds(1)-1
   do b2 = 0, ds(2)-1
      minval=min(Xcor_copy(b1,b2,:))
      maxval=max(Xcor_copy(b1,b2,:))
      if (maxval.gt.0. .and. minval.lt.0) then
          if (maxval-minval .gt. maxrange) then
	    maxrange=maxval-minval
            inds(1)=b1
	    inds(2)=b2
	  end if
      end if
   end do
end do
TS=Xcor(:,inds(1),inds(2))
wks  = gsn_open_wks("png","figs/EvolvingMixedCorrelation")
res@tiYAxisString = tostring(SmoothPer)+"-year correlation, basins "+tostring(inds(1)+1)+"/"+tostring(inds(2)+1)
plot=gsn_csm_xy(wks,Years,TS,res)
gsn_polyline(wks,plot,Years,(/TS/),polyres)   
   
do b1 = 0, ds(1)-1
   do b2 = 0, ds(2)-1
      Fname="figs/EvolvingCorrelationBasins_"+tostring(b1+1)+"-"+tostring(b2+1)
      InverseFname="figs/EvolvingCorrelationBasins_"+tostring(b2+1)+"-"+tostring(b1+1)+".png"
      if (b1.ne.b2) .and. .not. fileexists(InverseFname) then
         print("b1="+b1+", b2="+b2)      
	 TS=Xcor(:,b1,b2)
	 wks  = gsn_open_wks("png",Fname)
	 res@tiYAxisString = tostring(SmoothPer)+"-year correlation, basins "+tostring(b1+1)+"/"+tostring(b2+1)
	 plot=gsn_csm_xy(wks,Years,TS,res)
	 gsn_polyline(wks,plot,Years,(/TS/),polyres)        
      end if
   end do
end do






