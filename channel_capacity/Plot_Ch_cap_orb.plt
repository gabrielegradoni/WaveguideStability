#!/gnuplot
FILE="fort.21"
FILE1="fort.22"
FILE7="fort.27" 
 set terminal  pngcairo  enhanced  font   "Helvetica,   20"
 set output 'Fig1.png'
# set terminal postscript eps enhanced colour solid rounded "Helvetica" 25
# set output 'fig.ps'
set size 0.8,1
set title  "  " 
set xlabel "  n " 
set ylabel "  y  "

a1=0
a2=500
  set xrange [a1:a2]
#  set yrange [0:10]
#  set yrange [0:2.5]
#  set arrow from  a1,2 to a2,2    lw 2   lc rgb "grey"      nohead
#       set logscale y 
#       set format y '10^{%L}'
#       set logscale x 
#       set format x '10^{%L}'
       set nokey
f(x)=x*0.0015

# MEGNO for reg. orbits
  plot   FILE    u  1:4        w   lines      lt  1   lc   rgb "purple"       lw  2    t     ' log LE/n'  ,\
         FILE    u  1:5        w   lines      lt  1   lc   rgb "red"          lw  2    t     ' log RE/n'  ,\
	 FILE1   u  1:4        w   lines      lt  1   lc   rgb "pink"         lw  2    t     ' log LE/log n'  ,\
         FILE1   u  1:5        w   lines      lt  1   lc   rgb "orange"        lw  2    t     ' log RE/log n'  ,\
	 FILE7   u  1:4        w   lines      lt  1   lc   rgb "blue"         lw  2    t     ' log LE/log n'  ,\
         FILE7   u  1:5        w   lines      lt  1   lc   rgb "green"        lw  2    t     ' log RE/log n' 
	 