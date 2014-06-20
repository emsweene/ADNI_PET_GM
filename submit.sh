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

qsub -cwd -l mem_free=25G,h_vmem=30G,h_fsize=10G,h_stack=256M batch.sh $subject $dimx $dimy $dixz

done
done 
done
done
