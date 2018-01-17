set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Refinement9_100_6000_94500_0.dat"    using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_20.dat"   using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_21.dat"   using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_22.dat"   using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_23.dat"   using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_24.dat"   using 1:5 title " " with lines,\
"./Movie/Refinement9_100_6000_94500_25.dat"   using 1:5  title " " with lines,\

pause -1
