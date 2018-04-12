set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [1:100000]
#set xrange [-10:100]

plot \
"./Movie/Plot_256_100_0.dat"                    using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_10.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_20.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_30.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_40.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_50.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_60.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_70.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_80.dat"                   using 1:3    title "Initial Solution " with lines,\
"./Movie/Plot_256_100_90.dat"                   using 1:3    title "Initial Solution " with lines,\
#"./Movie/Plot13_256_500_0.dat"                   using 1:2    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1019.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1195.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1395.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1495.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1595.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1695.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1795.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_2000.dat"                using ($1):7    title "Initial Solution " with lines,\
"./Movie/Plot15_2560_100_0.dat"                    using 1:7    title "Initial Solution " with lines,\
"./Movie/Plot15_2560_100_41.dat"                   using 1:7    title "Initial Solution " with lines,\
"./Movie/Plot15_2560_100_42.dat"                   using 1:7    title "Initial Solution " with lines,\
"./Movie/Plot15_2560_100_43.dat"                   using 1:7    title "Initial Solution " with lines,\
"./Movie/Plot15_2560_100_2000.dat"                 using 1:7    title "Initial Solution " with lines,\
"./dat_saves/Pr075_Le03_Q9_B5_G14_mfcj_L124889_9453/Plot12_256_500_1260.dat" using 1:7    title "Initial Solution " with lines,\
"./Movie/Plot13_1024000_1_5.dat"                   using 1:4    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1260.dat"                using ($1-500):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1261.dat"                using ($1-500):7    title "Initial Solution " with lines,\
"./Movie/Plot12_256_500_1262.dat"                using ($1-500):7    title "Initial Solution " with lines,\



#"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_0.dat"      using 1:6  title "256-0" with lines,\
"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_2.dat"      using 1:6  title "256-2" with lines,\
"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_4.dat"      using 1:6  title "256-4" with lines,\
"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_6.dat"      using 1:6  title "256-6" with lines,\
"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_8.dat"      using 1:6  title "256-8" with lines,\
"./dat_saves/Implicit_CD_R256_D250_CE9/Plot_10.dat"     using 1:6  title "256-10" with lines,\
"./dat_saves/Implicit_CD_R1024_D250_CE10/Plot_2.dat"    using 1:6  title "1024-2" with lines,\
"./dat_saves/Implicit_CD_R1024_D250_CE10/Plot_4.dat"    using 1:6  title "1024-4" with lines,\
"./dat_saves/Implicit_CD_R1024_D250_CE10/Plot_6.dat"    using 1:6  title "1024-6" with lines,\
"./dat_saves/Implicit_CD_R1024_D250_CE10/Plot_8.dat"    using 1:6  title "1024-8" with lines,\
"./dat_saves/Implicit_CD_R1024_D250_CE10/Plot_10.dat"   using 1:6  title "1024-10" with lines,\

#"./Movie/Plot_4th_05_256_250_6_0.dat"   using 1:4       title "Initial Solution " with lines,\
"./Movie/Plot_4th_05_256_250_6_5.dat"    using 1:4       title "Initial Solution " with lines,\
"./Movie/Plot_4th_05_256_250_6_10.dat"   using 1:4       title "Initial Solution " with lines,\
#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Movie/Plot_4th_2_256_250_7_100.dat"   using 1:5       title "Initial Solution " with lines,\
"./Movie/Plot_4th_2_256_250_7_110.dat"   using 1:5       title "Initial Solution " with lines,\
"./Movie/Plot_2nd_4096_250_10_10.dat"    using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_2th_4096_250_10_10.dat"    using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_2nd_512_250_8_10.dat"      using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_4th_512_250_8_10.dat"      using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_2nd_128_250_6_10.dat"      using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_4th_128_250_6_10.dat"      using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_2nd_32_250_5_10.dat"       using 1:3   title "Initial Solution " with lines,\
"./Movie/Plot_4th_32_250_5_10.dat"       using 1:3   title "Initial Solution " with lines,\
"./dat_saves/Implicit_R512_D500_VC/Plot_6_0_N.dat"             using 1:5 title "Solution" with lines,\
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
