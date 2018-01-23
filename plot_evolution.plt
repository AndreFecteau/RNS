set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Plot_20_1000_95654_155.dat"          using 1:4 title " " with lines,\
"./Movie/Plot_20_1000_95654_160.dat"          using 1:4 title " " with lines,\
"./Movie/Plot_20_1000_95654_165.dat"          using 1:4 title " " with lines,\

#"./Movie/Test_500_95290_0.dat"           using 1:4 title " " with lines,\

pause -1
