set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
<<<<<<< HEAD
"./Movie/Test_Explicit_Residual_20000_0.dat"           using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_20000_10.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_20000_20.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_20000_30.dat"          using 1:4 title "Initial Solution " with lines,\
=======
"./Movie/Test_Explicit_Residual_10000_0.dat"           using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_10000_0.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_10000_1.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_10000_2.dat"          using 1:4 title "Initial Solution " with lines,\
>>>>>>> d444da6efe4b69421316ad4a1d5fab52dd459cec




#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\










pause -1
