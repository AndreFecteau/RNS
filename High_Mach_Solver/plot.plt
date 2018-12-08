# set title "Lambda Resolution Study "
set xlabel "Cells per Flame Length"
# set ylabel "Lambda"
# set key at 300, 124850
# set key title "Domain"
# set yrange [0.06:0.08]
set xrange [22:27]

plot \
"./RK4_Usual.dat"               using ($1):($3)    title "_D500" with lines ,\
"./RK4_Usual_Burned.dat"        using ($1):($3)    title "_D500" with lines ,\
# "./RK4_Usual_upstream.dat"      using ($1):($5)    title "_D500" with lines ,\
"./Le_1_64FL_182.dat"           using ($1-230.02):($5)    title "_D500" with lines ,\
# "./RK4_Usual_upstream.dat"      using ($1):($2)    title "_D500" with lines ,\
# "./RK4_Usual_upstream.dat"      using ($1):($3)    title "_D500" with lines ,\
# "./RK4_Usual.dat"      using ($1+37.7):($2)    title "_D500" with lines ,\
# "./RK4_Usual.dat"      using ($1+37.7):($3)    title "_D500" with lines ,\
# "./RK40.dat"      using ($1):($2)    title "_D500" with lines ,\
# "./RK40.dat"      using ($1):($3)    title "_D500" with lines ,\
# "./RK40.dat"      using ($1):($4)    title "_D500" with lines ,\

pause -1
