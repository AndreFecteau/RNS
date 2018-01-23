set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Plot_Euler008_500_95654_0.dat"          using 1:4 title " " with lines,\
"./Movie/Plot_Euler008_500_95654_1.dat"          using 1:4 title " " with lines,\

pause -1
