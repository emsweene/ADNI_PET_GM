#!/bin/bash

cd /dexter/disk1/smart/AD/ADNI/Greven_ADNI/Time_Mem

subjects='10  
100  
1000
10000'


for subject in $subjects 
do 

dimxs='1
10
100
1000'

for dimx in $dimxs 
do  

dimys='1
10
100
1000'

for dimy in $dimys
do  

dimzs='1
10
100
1000'

for dimz in $dimzs
do  


echo $subject
echo $dimx 
echo $dimy
echo $dimz

qsub -cwd -l mem_free=100G,h_vmem=100G, batch.sh $subject $dimx $dimy $dimz

done
done 
done
done
