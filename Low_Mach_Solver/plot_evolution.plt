# set title "Lambda Resolution Study "
set xlabel "Cells per Flame Length"
# set ylabel "Lambda"
# set key at 300, 124850
# set key title "Domain"
set yrange [0:20]
# set xrange [240:280]

plot \
"./RK4_Usual.dat"      using ($1+236.9):($3)    title "_D500" with lines ,\
"./Solution.dat"       using ($1):($3)    title "_D500" with lines ,\
# "./RK4.dat"      using ($1):($3)    title "_D500" with lines ,\
#////

pause -1
