set logscale y
set logscale x
set yrange[0.001:1.e50]
set xrange[1.e-12:1]
plot  "ML.out" using 2:5 title "rad" with lines, "ML.out" using 2:6 title "b" with lines, "ML.out" using 2:8 title "de" with lines, "MQ.out" using 2:7 title "phi1" with lines, "MQ1.out" using 2:7 title "phi" with lines, "MQ3.out" using 2:7 title "phi" with lines


