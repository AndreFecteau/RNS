set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
#set xrange [-10:100]

plot \
"./Movie/Case_1_D500_L300_R256_0.dat"                  using ($1):6    title "Initial Solution " with lines,\
"./Movie/Case_1_D500_L300_R256_1.dat"                  using ($1):6    title "Initial Solution " with lines,\
"./Movie/Case_1_D500_L300_R256_2.dat"                  using ($1):6    title "Initial Solution " with lines,\

pause -1
