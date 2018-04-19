#!/bin/bash
pkill gnuplot

gnuplot -p << EOP
#set terminal wxt enhanced font 'Verdana,10' persist
set terminal wxt 


set pm3d map
#set palette gray
set samples 100; set isosamples 100



#set zrange [2990:3010]
splot "matrix.txt" matrix  
EOP

#eog data_1.jpg 
