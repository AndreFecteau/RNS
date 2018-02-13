set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Test_Explicit_Residual_100000_1000000_0.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_100000_1000000_10.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_100000_1000000_20.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_100000_1000000_30.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_100000_1000000_40.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_100000_1000000_50.dat"          using 1:4 title "Initial Solution " with lines,\




#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\










pause -1
