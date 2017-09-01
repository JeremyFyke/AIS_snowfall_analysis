#!/bin/bash

vars=(PRECC PRECL)
vars2=(ICEFRAC TREFHT PSL TMQ PS Z3 UQ VQ)
vars2=()

permdir=/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/AIS_snowfall_ozone_connection/
syr=1955
eyr=2005


for nrr in $1; do

   echo '           #########Processing ensemble number: '$nrr '#########'

#FIXED-O3 LE PROCESSING
   
   rm -f tmp*.nc

   for vr in ${vars[@]}; do
      cdo -a setcalendar,standard /glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/$vr/b.e11.B20LE_fixedO3_$nrr.cam.h0.$vr.195501-200512.nc tmp$vr".nc"
   done

   cdo merge tmpPRECC.nc tmpPRECL.nc tmp3.nc

   cdo expr,'PRECIP=PRECL+PRECC' tmp3.nc tmp4.nc

   cdo -b F64 mulc,86400000. tmp4.nc tmp5.nc # m/s to mm/day
   cdo -b F64 muldpm tmp5.nc tmp.nc # to mm/month
   cdo -b F64 selyear,"$syr"/"$eyr" tmp.nc mon_PRECIP_$nrr".nc" 
   cdo -b F64 yearsum -selseas,DJF tmp.nc tmp2.nc
   cdo -b F64 selyear,"$syr"/"$eyr" tmp2.nc DJF_PRECIP_$nrr".nc"
   cdo -b F64 yearsum mon_PRECIP_$nrr".nc" year_PRECIP_$nrr".nc" # to mm/year

   mv *PRECIP_"$nrr".nc $permdir/files_"$syr"_"$eyr"/no-ozone

   for vr in ${vars2[@]}; do
      cdo -a setcalendar,standard /glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/$vr/b.e11.B20LE_fixedO3_$nrr.cam.h0.$vr.195501-200512.nc tmp$vr".nc"
      cdo -b F64 selyear,"$syr"/"$eyr" tmp$vr".nc" mon_$vr"_"$nrr".nc" 
      cdo -b F64 yearmean -selseas,DJF mon_$vr"_"$nrr".nc" DJF_$vr"_"$nrr".nc"
      cdo -b F64 yearmean mon_$vr"_"$nrr".nc" year_$vr"_"$nrr".nc"
      mv -f *$vr"_"$nrr".nc" $permdir/files_"$syr"_"$eyr"/no-ozone
   done
   
#CESM LE PROCESSING

   rm -f tmp*.nc

   for vr in ${vars[@]}; do      
      if [ $nrr == 001 ]; then
         cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.185001-200512.nc tmp$vr".nc"
      else 
         cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.192001-200512.nc tmp$vr".nc"
      fi
   done

   cdo merge tmpPRECC.nc tmpPRECL.nc tmp3.nc

   cdo expr,'PRECIP=PRECL+PRECC' tmp3.nc tmp4.nc

   cdo -b F64 mulc,86400000. tmp4.nc tmp5.nc # m/s to mm/day
   cdo -b F64 muldpm tmp5.nc tmp.nc # to mm/month
   cdo -b F64 selyear,"$syr"/"$eyr" tmp.nc mon_PRECIP_$nrr".nc" 
   cdo -b F64 yearsum -selseas,DJF tmp.nc tmp2.nc
   cdo -b F64 selyear,"$syr"/"$eyr" tmp2.nc DJF_PRECIP_$nrr".nc"
   cdo -b F64 yearsum mon_PRECIP_$nrr".nc" year_PRECIP_$nrr".nc" # to mm/year

   mv *PRECIP_"$nrr".nc $permdir/files_"$syr"_"$eyr"/lens

   for vr in ${vars2[@]}; do
      if [ $nrr == 001 ]; then
         cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.185001-200512.nc tmp$vr".nc"
      else    
         cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.192001-200512.nc tmp$vr".nc"
      fi      
      cdo -b F64 selyear,"$syr"/"$eyr" tmp$vr".nc" mon_$vr"_"$nrr".nc" 
      cdo -b F64 yearmean -selseas,DJF mon_$vr"_"$nrr".nc" DJF_$vr"_"$nrr".nc"
      cdo -b F64 yearmean mon_$vr"_"$nrr".nc" year_$vr"_"$nrr".nc"
      mv -f *$vr"_"$nrr".nc" $permdir/files_"$syr"_"$eyr"/lens
   done
   
done

rm -f tmp*.nc
