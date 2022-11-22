
        set terminal png size 600,400 truecolor
        set output "./bamstat_plots/mother-indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "mother"
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	5705
2	1471
3	878
4	860
5	316
6	158
7	65
8	3
9	42
10	51
11	0
12	53
13	45
14	4
15	0
16	2
17	0
18	0
19	0
20	0
21	19
22	12
24	0
25	0
27	0
30	0
end
1	10110
2	2799
3	1011
4	896
5	274
6	307
7	27
8	77
9	16
10	79
11	32
12	32
13	1
14	44
15	8
16	14
17	2
18	20
19	13
20	1
21	16
22	1
24	1
25	28
27	20
30	17
end
1	0.564293
2	0.525545
3	0.868447
4	0.959821
5	1.153285
6	0.514658
7	2.407407
8	0.038961
9	2.625000
10	0.645570
11	0.000000
12	1.656250
13	45.000000
14	0.090909
15	0.000000
16	0.142857
17	0.000000
18	0.000000
19	0.000000
20	0.000000
21	1.187500
22	12.000000
24	0.000000
25	0.000000
27	0.000000
30	0.000000
end
