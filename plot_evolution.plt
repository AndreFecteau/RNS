set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Test2_500_95290_0.dat"           using 1:4 title " " with lines,\
"./Movie/Test2_500_95290_1.dat"           using 1:4 title " " with lines,\
"./Movie/Test2_500_95290_2.dat"           using 1:4 title " " with lines,\
"./Movie/Test2_500_95290_3.dat"           using 1:4 title " " with lines,\

pause -1
