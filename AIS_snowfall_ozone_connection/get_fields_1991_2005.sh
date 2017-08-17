#!/bin/csh -fx


#loop over ensemble members
set nr = 1
while ( $nr < 100 )
   if ( $nr < 10) then
      set nrr = "00"$nr
   else
      set nrr = "0"$nr
   endif

   rm -f tmp*.nc

   set vars = ( PRECC PRECL )
   set vars2 = ( ICEFRAC TREFHT PSL TMQ )

   foreach vr ( $vars )
      #cdo shifttime,-1mo ../../../data/fixed-ozone/$vr/b.e11.B20LE_fixedO3_$nrr.cam.h0.$vr.195501-200512.nc tmp$vr".nc"
      cdo -a setcalendar,standard /glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/$vr/b.e11.B20LE_fixedO3_$nrr.cam.h0.$vr.195501-200512.nc tmp$vr".nc"
   end
   cdo merge tmpPRECC.nc tmpPRECL.nc tmp3.nc

   cdo expr,'PRECIP=PRECL+PRECC' tmp3.nc tmp4.nc

   cdo -b F64 mulc,86400000. tmp4.nc tmp5.nc # m/s to mm/day
   cdo -b F64 muldpm tmp5.nc tmp.nc # to mm/month
   cdo -b F64 selyear,1991/2005 tmp.nc mon_PRECIP_$nrr".nc" 
   cdo -b F64 yearsum -selseas,DJF tmp.nc tmp2.nc
   cdo -b F64 selyear,1991/2005 tmp2.nc DJF_PRECIP_$nrr".nc"
   cdo -b F64 yearsum mon_PRECIP_$nrr".nc" year_PRECIP_$nrr".nc" # to mm/year

   foreach vr ( $vars2 )
   cdo -a setcalendar,standard /glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/$vr/b.e11.B20LE_fixedO3_$nrr.cam.h0.$vr.195501-200512.nc tmp$vr".nc"
   cdo -b F64 selyear,1991/2005 tmp$vr".nc" mon_$vr"_"$nrr".nc" 
   cdo -b F64 yearmean -selseas,DJF mon_$vr"_"$nrr".nc" DJF_$vr"_"$nrr".nc"
   cdo -b F64 yearmean mon_$vr"_"$nrr".nc" year_$vr"_"$nrr".nc"
   end

   mv *$nrr*nc files_1991_2005/no-ozone/

   @ nr = $nr + 1

end

#loop over ensemble members
set nr = 1
while ( $nr < 200 )
if ( $nr < 10) then
set nrr = "00"$nr
else
if ( $nr < 100) then
set nrr = "0"$nr
else
set nrr = ""$nr
endif
endif

rm -f tmp*.nc

set vars = ( PRECC PRECL )
set vars2 = ( ICEFRAC TREFHT PSL TMQ )

foreach vr ( $vars )

if ( $nr == 1 ) then
	cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.185001-200512.nc tmp$vr".nc"
else
	cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.192001-200512.nc tmp$vr".nc"
endif

end

cdo merge tmpPRECC.nc tmpPRECL.nc tmp3.nc

cdo expr,'PRECIP=PRECL+PRECC' tmp3.nc tmp4.nc

cdo -b F64 mulc,86400000. tmp4.nc tmp5.nc # m/s to mm/day
cdo -b F64 muldpm tmp5.nc tmp.nc # to mm/month
cdo -b F64 selyear,1991/2005 tmp.nc mon_PRECIP_$nrr".nc" 
cdo -b F64 yearsum -selseas,DJF tmp.nc tmp2.nc
cdo -b F64 selyear,1991/2005 tmp2.nc DJF_PRECIP_$nrr".nc"
cdo -b F64 yearsum mon_PRECIP_$nrr".nc" year_PRECIP_$nrr".nc" # to mm/year

foreach vr ( $vars2 )
if ( $nr == 1 ) then
	cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.185001-200512.nc tmp$vr".nc"
else
	cdo -a setcalendar,standard /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/$vr/b.e11.B20TRC5CNBDRD.f09_g16.$nrr.cam.h0.$vr.192001-200512.nc tmp$vr".nc"
endif
cdo -b F64 selyear,1991/2005 tmp$vr".nc" mon_$vr"_"$nrr".nc" 
cdo -b F64 yearmean -selseas,DJF mon_$vr"_"$nrr".nc" DJF_$vr"_"$nrr".nc"
cdo -b F64 yearmean mon_$vr"_"$nrr".nc" year_$vr"_"$nrr".nc"
end
mv *$nrr*nc files_1991_2005/lens/

@ nr = $nr + 1


end
