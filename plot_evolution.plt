set title ""
set xlabel "x"
set ylabel ""
set key title ""
#set yrange [-0.1:1.1]
set xrange [-10:10]
plot \
"./Movie/Plot_HLLE_5_64_250_4_0.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_5_64_250_4_21.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_5_64_250_4_22.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_5_64_250_4_24.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_5_64_250_4_26.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_5_64_250_4_28.dat"   using 1:4 title "Initial Solution " with lines,\






#"./Tests/Problem4_Exact.dat"                          using 1:3  title " " with lines,\
"./Movie/Plot_CD_2000_16_250_4_10.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_CD_2000_16_500_4_10.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_CD_2000_16_1000_4_10.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_CD_2000_16_2000_4_10.dat"   using 1:4 title "Initial Solution " with lines,\
"./dat_saves/Implicit_CD_Resolution/Plot_32_3_0_N.dat"  using 1:4 title "Initial Condition " with lines,\
"./dat_saves/Implicit_CD_Resolution/Plot_4_2_5_N.dat"   using 1:4 title "Resolution 4 " with lines,\
"./dat_saves/Implicit_CD_Resolution/Plot_8_3_5_N.dat"   using 1:4 title "Resolution 8 " with lines,\
"./dat_saves/Implicit_CD_Resolution/Plot_16_3_5_N.dat"  using 1:4 title "Resolution 16" with lines,\
"./dat_saves/Implicit_CD_Resolution/Plot_32_5_3_N.dat"  using 1:4 title "Resolution 32" with lines,\
"./Tests/Problem4_HLLE.dat"                           using 1:3  title " " with lines,\
"./Tests/Problem4_S_HLLE.dat"                           using 1:3  title " " with lines,\


#"./Movie/Plot_HLLE_2000_16_100_4_0.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_16_100_4_1.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_32_100_5_1.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_64_100_5_5.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_256_100_5_5.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_256_100_5_2.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_CD_2000_64_100_5_7.dat"   using 1:4 title "Initial Solution " with lines,\
"./Movie/Plot_HLLE_2000_256_100_6_10.dat"   using 1:4 title "Initial Solution " with lines,\








pause -1
