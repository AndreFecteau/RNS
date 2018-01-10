set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement1_30_700_80000_0.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement1_30_700_80000_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement1_30_700_80000_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement1_30_700_80000_3.dat"   using 1:2 title " " with lines,\

pause -1
