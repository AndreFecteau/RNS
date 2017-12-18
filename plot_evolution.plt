set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement_4_6000_94500_0.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_4_6000_94500_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_4_6000_94500_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_4_6000_94500_3.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_4_6000_94500_4.dat"   using 1:2 title " " with lines,\

pause -1
