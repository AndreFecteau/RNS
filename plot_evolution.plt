set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
set xrange [-10:10]

plot \
"./Movie/Plot_HLLE_64_250_3_0.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_3_4.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_3_8.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_3_13.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_128_250_3_16.dat"   using 1:4 title "Initial Solution " with lines,\


#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./dat_saves/Implicit_R512_D500_VC/Plot_6_300_N.dat"   using 1:6 title "Solution" with lines,\
"./dat_saves/Implicit_R512_D500_VC/Plot_6_0_N.dat"     using 1:6 title "Initial_Conditions " with lines,\
"./dat_saves/Implicit_HLLE_Resolution_D250/Plot_64_28_N.dat"     using 1:4 title "Domaine 100 " with lines,\
"./dat_saves/Implicit_HLLE_Resolution_D250/Plot_128_16_N.dat"    using 1:4 title "Domaine 250 " with lines,\
"./dat_saves/Implicit_HLLE_Resolution_D250/Plot_256_49_N.dat"    using 1:4 title "Domaine 500 " with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\

g = -2
f(x) = g*x+9
#plot "./dat_saves/Implicit_HLLE_Convergence/Convergence_Plot.dat"   using (log($1)):(log($2)) title "Convergence Plot" with lines,\
f(x)










pause -1
