set title "Effects of changing the CFL"
set xlabel "x"
set ylabel "p"
set key title "CFL"
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Test_Implicit_Residual_750_0.dat"           using 1:4 title "Exact Solution " with lines,\
"./Movie/Test_Implicit_Residual_750_2.dat"           using 1:4 title "500 Cells 10's " with lines,\
"./Movie/Test_Implicit_Residual_750_4.dat"           using 1:4 title "500 Cells 20's " with lines,\
"./Movie/Test_Implicit_Residual_750_6.dat"           using 1:4 title "500 Cells 20's " with lines,\
"./Movie/Test_Implicit_Residual_750_8.dat"           using 1:4 title "500 Cells 20's " with lines,\
"./Movie/Test_Implicit_Residual_750_10.dat"          using 1:4 title "500 Cells 20's " with lines,\
"./Movie/Test_Implicit_Residual_750_11.dat"          using 1:4 title "500 Cells 20's " with lines,\
"./Movie/Test_Implicit_Residual_750_12.dat"          using 1:4 title "500 Cells 20's " with lines,\




#"./Tests/Problem4_HLLE.dat"               using 1:3  title " " with lines,\
"./Tests/Problem4_Exact.dat"               using 1:3  title " " with lines,\


pause -1
