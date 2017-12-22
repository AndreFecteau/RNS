set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
<<<<<<< HEAD
"./Movie/Refinement_2_1000_94500_0.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_93250_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_93250_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92875_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92875_2.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92687_1.dat"   using 1:2 title " " with lines,\
"./Movie/Refinement_2_1000_92687_2.dat"   using 1:2 title " " with lines,\
=======
"./Movie/Hello10_6000_93250_0.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_93250_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_93250_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_91625_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_91625_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90812_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90812_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90406_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90406_4.dat"   using 1:2 title " " with lines,\
>>>>>>> 0872bb68967b239678d75c20c1faa3a2225e94f0

pause -1
