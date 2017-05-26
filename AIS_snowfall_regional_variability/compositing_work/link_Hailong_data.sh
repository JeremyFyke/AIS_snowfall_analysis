#Copy and rename composite files made by Hailong, in a way that allows them to be read by analysis scripts.
for b in $(seq -f "%02g" 1 27); do
   for c in high low; do
      rm output/"$c"/"$b"/Ann_column_UQ.nc
      rm output/"$c"/"$b"/Ann_column_VQ.nc      
   done
done
