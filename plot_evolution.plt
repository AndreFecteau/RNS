set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Plot_20_1000_95654_134.dat"          using 1:4 title " " with lines,\
"./Movie/Plot_20_1000_95654_135.dat"          using 1:4 title " " with lines,\
"./Movie/Plot_20_1000_95654_136.dat"          using 1:4 title " " with lines,\

pause -1
