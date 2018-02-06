set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Test_Explicit_Residual_1500_0.dat"           using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_1500_4.dat"           using 1:4 title "1500  Cells" with lines,\
"./Movie/Test_Explicit_Residual_1500_5.dat"           using 1:4 title "1500  Cells" with lines,\
"./Movie/Test_Explicit_Residual_1500_6.dat"           using 1:4 title "1500  Cells" with lines,\
"./Movie/Test_Explicit_Residual_1500_7.dat"           using 1:4 title "1500  Cells" with lines,\




#"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\






pause -1
