Dirs=   (/"lens","no-ozone"/)
EnsNums=(/1,10,12,15,18,20,22,28/)
AvgPer= (/"year","DJF"/)
nEns=dimsizes(EnsNums)
do d=0,dimsizes(Dirs)-1
   do e=0,nEns-1
      elong=sprinti("%0.3i", EnsNums(e))
      do s=0,dimsizes(AvgPer)-1

         print("Loading Z3")
         fin=str_concat((/"files_1986_2005/",Dirs(d),"/",AvgPer(s),"_Z3_",elong,".nc"/))
	 in=addfile(fin,"r")
	 Z=in->Z3
	 
	 print("Loading PS")
	 fin=str_concat((/"files_1986_2005/",Dirs(d),"/",AvgPer(s),"_PS_",elong,".nc"/))
	 in=addfile(fin,"r")
	 PSFC=in->PS ;Pa
	 
	 print("Loading hyam,hybm")
	 fin="/glade/p/umit0006/FESD/fixedO3/atm/proc/tseries/monthly/PS/b.e11.B20LE_fixedO3_001.cam.h0.PS.195501-200512.nc"
         in=addfile(fin,"r")
	 hyam = in->hyam
	 hybm = in->hybm

	 P0mb=(100000.)/100. ; Apparently this is the standard CESM reference pressure.  convert to hPa/mb

	 Znew = vinth2p (Z,hyam,hybm,(/ 850.0,500.0 /),PSFC,1,P0mb,1,True)

	 fout=str_concat((/"files_1986_2005/",Dirs(d),"/",AvgPer(s),"_Z500_850_",elong,".nc"/))
	 print("Writing "+fout)
	 system("rm -f "+fout)
         ncdf=addfile(fout ,"c")
	 ncdf->Z580=Znew(:,0,:,:)
	 ncdf->Z500=Znew(:,1,:,:)

      end do
   end do
end do
