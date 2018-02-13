set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1]
#set xrange [0:1]
plot \
"./Movie/Test_Explicit_Residual_16_0_0.dat"           using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_15.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_16.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_17.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_18.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_19.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_16_0_110.dat"          using 1:4 title "Initial Solution " with lines,\
"./Movie/Test_Explicit_Residual_0_10000_20.dat"          using 1:4 title "Initial Solution " with lines,\





#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\










pause -1
