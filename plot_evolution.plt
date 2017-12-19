set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement_2_1000_94500_0.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_93250_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_93250_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92875_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92875_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92687_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92687_2.dat"   using 1:2 title " " with lines,\

pause -1
