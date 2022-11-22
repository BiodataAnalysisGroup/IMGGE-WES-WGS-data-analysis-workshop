
        set terminal png size 600,400 truecolor
        set output "./bamstat_plots/proband-indel-dist.png"
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
        set title "proband"
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	8559
2	2028
3	1361
4	901
5	348
6	257
7	97
8	33
9	44
10	23
11	1
12	105
13	68
14	5
15	0
16	0
17	0
18	12
19	0
20	1
21	45
22	8
24	8
25	0
27	0
30	0
end
1	14021
2	3785
3	991
4	1402
5	438
6	460
7	102
8	169
9	25
10	152
11	41
12	48
13	5
14	129
15	12
16	9
17	4
18	5
19	15
20	2
21	22
22	21
24	1
25	13
27	31
30	9
end
1	0.610441
2	0.535799
3	1.373360
4	0.642653
5	0.794521
6	0.558696
7	0.950980
8	0.195266
9	1.760000
10	0.151316
11	0.024390
12	2.187500
13	13.600000
14	0.038760
15	0.000000
16	0.000000
17	0.000000
18	2.400000
19	0.000000
20	0.500000
21	2.045455
22	0.380952
24	8.000000
25	0.000000
27	0.000000
30	0.000000
end
