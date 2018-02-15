set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Plot_HLLE_Resolution_16_2000_3_0.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_Resolution_16_2000_3_1.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_Resolution_16_2000_3_2.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_Resolution_16_2000_3_3.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_Resolution_16_2000_3_4.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_Resolution_16_2000_3_5.dat"   using 1:4 title "Initial Solution " with lines,\





#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\










pause -1
