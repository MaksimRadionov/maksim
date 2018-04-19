#!/bin/bash
pkill gnuplot_qt
rm ./matrix.txt
make all_max && ./all_max < date.dat $1 

gnuplot -p << EOP

#set terminal jpeg size 640,480
set terminal wxt
#set output "data_1.jpg"

splot "matrix.txt" matrix with dots 
EOP

#eog data_1.jpg 
