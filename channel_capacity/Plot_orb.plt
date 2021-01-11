#!/gnuplot
FILE="fort.20"
FILE1="fort.15"
FILE2="Ch_cap_in_points.dat"
 set terminal  pngcairo  enhanced  font   "Helvetica,   20"
 set output 'Fig.png'
# set terminal postscript eps enhanced colour solid rounded "Helvetica" 25
# set output 'fig.ps'
#-------------------------------------------------------------
#  Plot delle orbite per waveguide 2D
#  In rosso i punti iniziali da FILE1
#------------------------------------------------------------
set size 0.8, 1
set title  "   " 
set xlabel "  x/(2 {/Symbol p}) " 
set ylabel "  v_x  "
a=0.5
ay=0.15
set xrange [-a:a]
set yrange [-1:1]
eps=0.25
D=2*sqrt(eps)
 set arrow from   -a,D  to a,D    lw 1   lc rgb "black"      nohead
  set arrow from  -a,-D to a,-D    lw 1   lc rgb "black"      nohead
pi=3.1415926  
fp(x)= D*cos(pi*x)
fm(x)= -D*cos(pi*x)
#      set logscale y 
#      set format y '10^{%L}'
       set nokey
set xtics -.5,0.25,0.5       
plot FILE    u  2:3         w   points     pt  1   lc   rgb "blue"     ps  0.1   t    '    '     ,\
     FILE1   u  2:3         w   points     pt  5   lc   rgb "black"    ps  .5    t    '    '     ,\
     FILE2   u  2:3         w   points     pt  5   lc   rgb "red"      ps  .5    t    '    '    