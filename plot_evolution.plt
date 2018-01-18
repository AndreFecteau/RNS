set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement9_10_6000_95654_0.dat"    using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_3.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_4.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_8.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement9_10_6000_95654_9.dat"   using 1:2  title " " with lines,\

pause -1
