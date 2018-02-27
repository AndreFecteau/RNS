set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
set xrange [-10:10]

plot \
"./Movie/Plot_4th_4096_250_8_0.dat"    using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_4th_4096_250_8_7.dat"    using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_4th_4096_250_8_8.dat"    using 1:4 title "Initial Solution " with lines,\


#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./dat_saves/Implicit_R512_D500_VC/Plot_6_0_N.dat"             using 1:5 title "Solution" with lines,\
"./dat_saves/Implicit_HLLE_R128_D250_L92073/Plot_75_N.dat"     using 1:5 title "128-92073 " with lines,\
"./dat_saves/Implicit_HLLE_R512_D250_L93200/Plot_49_N.dat"     using 1:5 title "512-93200 " with lines,\
"./dat_saves/Implicit_HLLE_R4096_D250_L94900/Plot_15_N.dat"    using 1:5 title "4096-94900 " with lines,\
"./dat_saves/Implicit_HLLE_Resolution_D250/Plot_256_49_N.dat"  using 1:5 title "Flame Moving " with lines,\
"./Movie/Plot_HLLE_128_250_5_67.dat"   using 1:5 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_5_68.dat"   using 1:5 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_5_69.dat"   using 1:5 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_5_70.dat"   using 1:5 title "Initial Solution " with lines,\
"./dat_saves/Implicit_R512_D500_VC/Plot_6_0_N.dat"     using 1:6 title "Initial_Conditions " with lines,\
"./dat_saves/Implicit_HLLE_Resolution_D250/Plot_128_16_N.dat"    using 1:4 title "Domaine 250 " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\

g = -2
f(x) = g*x+9
#plot "./dat_saves/Implicit_HLLE_Convergence/Convergence_Plot.dat"   using (log($1)):(log($2)) title "Convergence Plot" with lines,\
f(x)










pause -1
