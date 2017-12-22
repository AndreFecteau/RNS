set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Hello10_6000_93250_0.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_93250_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_93250_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_91625_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_91625_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90812_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90812_4.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90406_2.dat"   using 1:2 title " " with lines,\
"./Movie/Hello10_6000_90406_4.dat"   using 1:2 title " " with lines,\

pause -1
