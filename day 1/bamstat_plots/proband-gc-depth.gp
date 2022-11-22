
            set terminal png size 600,500 truecolor
            set output "./bamstat_plots/proband-gc-depth.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Mapped depth"
            set xlabel "Percentile of mapped sequence ordered by GC content"
            set x2label "GC Content [%]"
            set title "proband"
            set x2tics ("30" 4.681,"40" 45.936,"50" 76.183)
            set xtics nomirror
            set xrange [0.1:99.9]

            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#dedede" t '10-90th percentile' , \
                 '-' using 1:2:3 with filledcurve lt 1 lc rgb "#bbdeff" t '25-75th percentile' , \
                 '-' using 1:2 with lines lc rgb "#0084ff" t 'Median'
        0.309	0.000	0.005
0.412	0.005	0.010
0.463	0.005	0.005
0.617	0.005	0.010
0.772	0.005	0.015
0.823	0.005	0.005
0.926	0.005	0.005
1.080	0.005	0.010
1.286	0.005	0.005
1.543	0.005	0.020
1.646	0.005	0.010
1.800	0.005	0.005
1.852	0.005	0.005
1.903	0.005	0.005
2.160	0.005	0.417
2.212	0.010	0.010
2.623	0.005	3.131
3.035	0.005	28.474
3.498	0.005	10.849
4.681	0.005	23.456
6.019	0.005	28.734
8.642	0.010	24.005
13.683	1.882	27.587
18.056	1.303	24.132
22.531	1.539	26.999
27.623	1.220	28.435
32.047	1.882	33.487
36.523	0.005	22.555
41.101	0.608	31.154
45.936	0.073	33.124
50.463	0.005	23.549
54.115	1.548	24.108
57.665	0.906	25.073
61.780	0.005	23.309
65.741	0.010	24.721
67.953	0.005	21.462
70.422	0.005	20.717
71.914	0.005	21.883
73.971	0.010	21.379
76.183	0.005	17.184
78.652	0.005	24.010
80.710	0.005	13.181
82.202	0.010	32.546
83.488	0.005	56.355
84.979	2.097	44.850
86.060	0.274	26.107
86.831	0.005	27.308
87.963	0.010	33.541
89.198	0.617	47.883
90.226	0.010	33.869
91.615	0.005	37.544
92.901	0.402	26.636
94.290	0.304	29.919
95.267	0.265	36.358
96.245	0.216	43.130
96.965	0.020	13.450
97.737	0.157	24.049
98.045	0.039	3.146
98.405	0.049	0.510
98.817	0.010	1.294
99.074	0.010	0.559
99.383	0.039	0.333
99.588	0.064	0.392
99.794	0.020	0.147
99.897	0.078	0.127
99.949	0.010	0.010
100.000	0.059	0.059
end
0.309	0.000	0.005
0.412	0.005	0.010
0.463	0.005	0.005
0.617	0.005	0.010
0.772	0.005	0.015
0.823	0.005	0.005
0.926	0.005	0.005
1.080	0.005	0.010
1.286	0.005	0.005
1.543	0.005	0.005
1.646	0.005	0.010
1.800	0.005	0.005
1.852	0.005	0.005
1.903	0.005	0.005
2.160	0.005	0.098
2.212	0.010	0.010
2.623	0.005	1.382
3.035	0.005	5.419
3.498	0.005	7.909
4.681	1.764	9.418
6.019	1.490	14.852
8.642	3.802	16.297
13.683	3.430	15.165
18.056	2.254	17.258
22.531	3.773	17.900
27.623	3.170	14.710
32.047	3.773	16.533
36.523	2.386	16.092
41.101	2.435	16.332
45.936	1.940	17.307
50.463	2.097	15.509
54.115	3.293	15.386
57.665	1.872	16.165
61.780	1.715	11.368
65.741	2.533	12.593
67.953	1.725	14.475
70.422	2.999	13.901
71.914	1.352	13.470
73.971	1.588	13.299
76.183	0.990	12.676
78.652	1.029	9.687
80.710	1.245	10.976
82.202	3.606	13.558
83.488	1.617	20.070
84.979	3.734	13.152
86.060	1.372	10.506
86.831	1.352	18.003
87.963	0.029	9.859
89.198	1.225	32.761
90.226	1.666	15.802
91.615	0.598	17.586
92.901	1.421	9.535
94.290	1.303	16.944
95.267	2.259	28.435
96.245	0.960	9.457
96.965	0.412	4.704
97.737	0.319	1.225
98.045	0.039	2.097
98.405	0.294	0.392
98.817	0.049	0.265
99.074	0.010	0.176
99.383	0.039	0.304
99.588	0.064	0.348
99.794	0.020	0.039
99.897	0.078	0.127
99.949	0.010	0.010
100.000	0.059	0.059
end
0.309	0.005
0.412	0.005
0.463	0.005
0.617	0.005
0.772	0.010
0.823	0.005
0.926	0.005
1.080	0.010
1.286	0.005
1.543	0.005
1.646	0.005
1.800	0.005
1.852	0.005
1.903	0.005
2.160	0.010
2.212	0.010
2.623	1.107
3.035	1.798
3.498	2.538
4.681	4.175
6.019	3.626
8.642	7.438
13.683	7.673
18.056	7.551
22.531	9.653
27.623	8.134
32.047	8.644
36.523	7.654
41.101	7.409
45.936	7.683
50.463	5.454
54.115	7.967
57.665	7.605
61.780	5.096
65.741	6.086
67.953	4.410
70.422	7.664
71.914	5.067
73.971	7.120
76.183	2.773
78.652	2.347
80.710	3.464
82.202	8.693
83.488	6.880
84.979	7.732
86.060	5.743
86.831	5.199
87.963	3.616
89.198	8.644
90.226	11.486
91.615	5.370
92.901	4.665
94.290	6.047
95.267	11.932
96.245	5.478
96.965	1.186
97.737	0.823
98.045	0.470
98.405	0.353
98.817	0.167
99.074	0.167
99.383	0.108
99.588	0.265
99.794	0.020
99.897	0.078
99.949	0.010
100.000	0.059
end
