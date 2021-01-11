#!/gnuplot
FILE_IN="fort.30"
   set terminal  pngcairo  enhanced  font   "Helvetica,   20"
   set output 'Fig_LE.png'
#    set size  0.94,1         #  con x,  v negli assi 
     set size  0.94,1.08      #  con x_0,  v_0 negli assi 
# set terminal postscript eps enhanced colour solid rounded "Helvetica" # set output FILE_OUT
set title '  LE  '
set xlabel ' x_0/(2{/Symbol p})'
set ylabel ' v_{x0}'
set cblabel ' '
a=0.5
set xrange [-a:a]
set yrange [-1:1]
# set cbrange [1:2000]
  set cbrange [1:1E10]  # E5
  set palette model RGB defined (0 "dark-blue", 1  "blue",2 "green",3 "yellow", 4 "red")
 set logscale cb
 set format cb '10^{%L}'
 set xtics -.5,0.25,0.5
plot FILE_IN u 3:4:6  w image t '    '  #
# set arrow from  -a,0 to a,0    lw 1   lc rgb "black"      nohead
# set arrow from   0,-a to 0,a    lw 1   lc rgb "black"      nohead


