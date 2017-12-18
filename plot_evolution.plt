set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement_5_6000_92500_0.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_92500_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_92500_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_92500_3.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_93750_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_93750_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_5_6000_93750_3.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_6000_94350_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_6000_94350_2.dat"   using 1:2 title " " with lines,\


pause -1
