set title ""
set xlabel "Cells per Flame Legnth"
set ylabel "Lambda"
set key title ""
# set yrange [0:1.05]
set xrange [-10:350]

plot \
"./Movie/Delete_0.dat"      using ($1):4    title "_D500" with lines ,\
"./Movie/Delete_1.dat"      using ($1):4    title "_D500" with lines ,\
"./Movie/Delete_2.dat"      using ($1):4    title "_D500" with lines ,\
"./Movie/Delete_3.dat"      using ($1):4    title "_D500" with lines ,\
#"./Movie/Case_1_D500_L250_R256_0.dat"                using ($1):6    title "Initial Solution " with lines,\
"./Lambda_Plot_D500.dat"      using ($1):3    title "_D2500",\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R16/Solution_N.dat"      using ($1):7    title "_D500_F250_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R32/Solution_N.dat"      using ($1):7    title "_D500_F250_R32" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R64/Solution_N.dat"      using ($1):7    title "_D500_F250_R64" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R128/Solution_N.dat"     using ($1):7    title "_D500_F250_R128" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R256/Solution_N.dat"     using ($1):7    title "_D500_F250_R256" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D2000_F250_R16/Solution_N.dat"      using ($1):3    title "_D2000_F250_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D2000_F250_R16/Solution_N.dat"      using ($1):4    title "_D2000_F250_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D2000_F250_R16/Solution_N.dat"      using ($1):5    title "_D2000_F250_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D2000_F1000_R16/Solution_N.dat"     using ($1):4    title "_D2000_F1000_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R16/Solution_N.dat"       using ($1):4    title "_D500_F250_R16" with lines,\
"./dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D1000_F250_R16/Solution_N.dat"      using ($1):4    title "_D1000_F250_R16" with lines,\
"min_mod.dat"                using ($1):5    title "Initial Solution " with lines,\
"van_albada.dat"             using ($1):5    title "Initial Solution " with lines,\
"solution_vector.dat"        using ($1):5    title "Initial Solution " with lines,\
"solution_vector_rho.dat"    using ($1):5    title "Initial Solution " with lines,\
"./idontcare.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_3.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_6.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_9.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_12.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_15.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_18.dat"               using ($1):4    title "Initial Solution " with lines,\
"./Movie/Low_Mach_R16_21.dat"               using ($1):4    title "Initial Solution " with lines,\

pause -1
