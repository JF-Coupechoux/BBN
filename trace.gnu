set logscale y
set logscale x
set yrange[0.001:1.e50]
set xrange[1.e-12:1]
plot  "homogene.out" using 2:5 title "rad" with lines, "homogene.out" using 2:6 title "b" with lines, "homogene.out" using 2:7 title "phi" with lines
