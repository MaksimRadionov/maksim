#!/bin/bash
rm ./graph.txt
pkill gnuplot_qt
make script && ./script  < date2.dat $1 $2 $3
title="height${1}wigth$2"
gnuplot -p << EOP

#set terminal jpeg size 640,480
#set output "data_1.jpg"
set terminal wxt
set title "$title"
plot "graph2.txt" with lines,\
        "graph.txt" with lines 
EOP
#eog data_1.jpg 
