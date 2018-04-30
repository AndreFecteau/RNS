set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
#set xrange [-10:100]

plot \
"./dat_saves/CJ_Point_studies/Pr075_Le03_g14_Q9_b5/D2000_F1000_R16/Plot_0.dat"               using ($1):7    title "Initial Solution " with lines,\
"./dat_saves/CJ_Point_studies/Pr075_Le03_g14_Q9_b5/D2000_F1000_R16/Plot_49.dat"               using ($1):7    title "Initial Solution " with lines,\
"./dat_saves/CJ_Point_studies/Pr075_Le03_g14_Q9_b5/D2000_F1000_R16/Plot_69.dat"               using ($1):7    title "Initial Solution " with lines,\
#"./Movie/Case_1_D500_L250_R256_0.dat"                using ($1):6    title "Initial Solution " with lines,\
"./Movie/Case_1_D500_L250_R256_124.dat"               using ($1):3    title "Initial Solution " with lines,\

pause -1
