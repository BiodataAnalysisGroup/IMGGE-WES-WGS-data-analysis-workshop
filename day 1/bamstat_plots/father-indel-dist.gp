
        set terminal png size 600,400 truecolor
        set output "./bamstat_plots/father-indel-dist.png"
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
        set title "father"
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	6893
2	1750
3	1140
4	846
5	246
6	274
7	60
8	30
9	23
10	36
11	4
12	43
13	45
14	1
15	9
16	0
17	0
18	5
19	1
20	0
21	17
22	4
24	2
25	0
27	0
29	0
30	0
end
1	11667
2	3090
3	1120
4	859
5	375
6	450
7	195
8	172
9	51
10	125
11	5
12	28
13	3
14	38
15	21
16	4
17	1
18	31
19	8
20	4
21	12
22	22
24	0
25	23
27	18
29	1
30	19
end
1	0.590812
2	0.566343
3	1.017857
4	0.984866
5	0.656000
6	0.608889
7	0.307692
8	0.174419
9	0.450980
10	0.288000
11	0.800000
12	1.535714
13	15.000000
14	0.026316
15	0.428571
16	0.000000
17	0.000000
18	0.161290
19	0.125000
20	0.000000
21	1.416667
22	0.181818
24	0.000000
25	0.000000
27	0.000000
29	0.000000
30	0.000000
end
