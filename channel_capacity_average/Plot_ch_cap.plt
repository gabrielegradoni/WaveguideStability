#!/gnuplot
FILE="Ch_cap.dat"
  set terminal  pngcairo   enhanced  font "Helvetica, 16"
  set output 'Fig_Ch.png'
  set size 0.8,1 
set title  '  '         #         ' Channel capacity '
set ylabel ' C '
set xlabel '{/Symbol e} '
   set xrange [0:0.5]
   set yrange [0:1]
  
# set logscale y
# set format  y  '10^{%L}'
 set nokey
Nmax=25
Nmax=100

 f(x)= 1.2*(2*x)-0.4*(2*x)**2

#set ytics  0,0.05,0.15

plot   FILE    u  1:2         w   lines     lt  1   lc   rgb "blue"       lw  2   t    ' C-LE '   ,\
       FILE    u  1:3         w   lines     lt  1   lc   rgb "purple"     lw  2   t    ' C-RE '   ,\
       FILE    u  1:4         w   lines     lt  1   lc   rgb "red"        lw  2   t    ' C-RE '  ,\
       f(x)                   w   lines     lt  1   lc   rgb "black"      lw  3          