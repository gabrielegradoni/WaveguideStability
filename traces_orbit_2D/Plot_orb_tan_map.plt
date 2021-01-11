#!/gnuplot
FILE="fort.20"
FILE1="fort.15"
 set terminal  pngcairo  enhanced  font   "Helvetica,   20"
 set output 'Fig.png'
# set terminal postscript eps enhanced colour solid rounded "Helvetica" 25
# set output 'fig.ps'
#-------------------------------------------------------------
#  Plot delle orbite per waveguide 2D
#  In rosso i punti iniziali da FILE1
#------------------------------------------------------------
set size 0.75, 1
set title  "   " 
set xlabel "  x'=x/ (2{/Symbol p}) " 
set ylabel "  v_x  "
a=0.5
set xrange [-a:a]
set yrange [-1:1]
#  set yrange [0:a]
#  set arrow from  0,2 to a,2    lw 1   lc rgb "black"      nohead
#      set logscale y 
#      set format y '10^{%L}'
       set nokey     
      plot FILE     u  2:3         w   points     pt  1   lc   rgb "blue"     ps  0.1   t    '  orbits   '    #  ,\
#     FILE1    u  2:3         w   points     pt  5   lc   rgb "red"      ps  .5   t    ' initial conditions   '