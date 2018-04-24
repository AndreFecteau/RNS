set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
#set xrange [-10:100]

plot \
"./Movie/Case_3_D500_L250_R16_2000.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Case_1_D500_L250_R256_2000.dat"               using ($1):7    title "Initial Solution " with lines,\
#"./Movie/Case_1_D500_L250_R256_0.dat"                  using ($1):6    title "Initial Solution " with lines,\
"./Movie/Case_1_D500_L250_R256_1.dat"                  using ($1):6    title "Initial Solution " with lines,\

pause -1
