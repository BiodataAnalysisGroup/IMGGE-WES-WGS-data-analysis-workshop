
            set terminal png size 700,500 truecolor
            set output "./bamstat_plots/mother-quals2.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set multiplot
             set rmargin 0; set lmargin 0; set tmargin 0; set bmargin 0; set origin 0.1,0.1; set size 0.4,0.8
            set yrange [0:51]
            set ylabel "Quality"
            set xlabel "Cycle (fwd reads)"
            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#cccccc" t '25-75th percentile' , '-' using 1:2 with lines lc rgb "#000000" t 'Median', '-' using 1:2 with lines lt 1 t 'Mean'
        1	32	32
2	33	37
3	33	38
4	33	39
5	33	39
6	33	39
7	33	39
8	34	40
9	33	40
10	34	40
11	34	40
12	34	40
13	34	40
14	33	40
15	34	40
16	34	40
17	34	40
18	34	40
19	34	40
20	34	40
21	34	40
22	33	40
23	33	40
24	34	40
25	34	40
26	33	40
27	34	40
28	34	40
29	33	39
30	32	39
31	33	40
32	34	40
33	33	40
34	28	35
35	28	36
36	34	40
37	33	40
38	33	39
39	27	33
40	33	39
41	34	40
42	34	40
43	34	40
44	33	40
45	34	40
46	34	40
47	34	40
48	34	40
49	34	40
50	34	40
51	31	39
52	34	40
53	34	40
54	34	40
55	34	40
56	34	40
57	34	40
58	34	40
59	34	40
60	34	40
61	34	40
62	34	40
63	31	38
64	34	40
65	34	40
66	34	40
67	34	40
68	34	40
69	34	40
70	34	40
71	34	40
72	33	40
73	34	40
74	34	40
75	33	39
76	34	39
77	34	39
78	34	39
79	34	40
80	34	40
81	34	40
82	34	40
83	34	40
84	34	39
85	34	39
86	34	39
87	34	39
88	34	39
89	34	39
90	34	39
91	34	39
92	34	39
93	34	39
94	32	38
95	33	39
96	33	39
97	33	38
98	33	38
99	32	37
100	32	37
101	32	36
end
1	32
2	34
3	34
4	35
5	34
6	34
7	34
8	36
9	36
10	36
11	36
12	35
13	36
14	35
15	36
16	36
17	36
18	36
19	36
20	36
21	36
22	36
23	35
24	36
25	36
26	35
27	36
28	36
29	35
30	34
31	35
32	36
33	36
34	32
35	32
36	36
37	36
38	35
39	30
40	35
41	36
42	36
43	36
44	36
45	36
46	36
47	36
48	36
49	36
50	36
51	34
52	36
53	36
54	36
55	36
56	36
57	36
58	36
59	36
60	36
61	36
62	36
63	34
64	36
65	36
66	36
67	36
68	36
69	36
70	36
71	36
72	35
73	36
74	35
75	35
76	35
77	36
78	36
79	36
80	36
81	36
82	36
83	36
84	36
85	35
86	36
87	35
88	35
89	35
90	35
91	35
92	35
93	35
94	34
95	35
96	35
97	34
98	34
99	34
100	34
101	34
end
1	32.73
2	35.59
3	35.71
4	36.75
5	36.43
6	35.71
7	36.09
8	37.39
9	37.36
10	37.39
11	37.40
12	37.05
13	37.43
14	36.55
15	37.55
16	37.54
17	37.53
18	37.58
19	37.62
20	37.45
21	37.68
22	37.28
23	36.84
24	37.37
25	37.70
26	36.51
27	37.69
28	37.23
29	36.65
30	35.41
31	36.90
32	37.74
33	37.04
34	32.35
35	32.90
36	37.42
37	37.16
38	36.34
39	31.59
40	36.57
41	37.71
42	37.50
43	37.44
44	36.99
45	37.33
46	37.72
47	37.63
48	37.72
49	37.60
50	37.69
51	35.64
52	37.64
53	37.68
54	37.74
55	37.76
56	37.56
57	37.58
58	37.63
59	37.67
60	37.57
61	37.49
62	37.49
63	34.65
64	37.30
65	37.53
66	37.55
67	37.47
68	37.44
69	37.41
70	37.45
71	37.37
72	36.34
73	37.22
74	37.30
75	36.02
76	36.69
77	37.18
78	37.32
79	37.20
80	37.31
81	37.29
82	37.30
83	37.31
84	37.25
85	37.21
86	37.21
87	37.08
88	36.95
89	37.04
90	37.00
91	37.01
92	36.96
93	36.93
94	35.19
95	36.77
96	36.58
97	36.16
98	36.08
99	35.00
100	35.61
101	34.43
end

                set origin 0.55,0.1
                set size 0.4,0.8
                unset ytics
                set y2tics mirror
                set yrange [0:51]
                unset ylabel
                set xlabel "Cycle (rev reads)"
                set label "mother" at screen 0.5,0.95 center
                plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#cccccc" t '25-75th percentile' , '-' using 1:2 with lines lc rgb "#000000" t 'Median', '-' using 1:2 with lines lt 2 t 'Mean'
            1	32	32
2	32	36
3	33	37
4	32	39
5	32	39
6	32	39
7	33	39
8	33	39
9	32	39
10	32	39
11	32	39
12	33	39
13	33	39
14	33	39
15	33	39
16	33	39
17	33	39
18	33	40
19	34	40
20	32	39
21	33	40
22	33	40
23	34	40
24	32	39
25	33	40
26	34	40
27	34	40
28	34	40
29	34	40
30	34	40
31	34	40
32	34	40
33	34	40
34	32	39
35	34	40
36	33	39
37	34	40
38	34	40
39	34	40
40	34	40
41	34	40
42	34	40
43	33	39
44	34	40
45	34	40
46	34	40
47	34	40
48	33	39
49	34	40
50	31	38
51	34	39
52	34	40
53	34	40
54	34	40
55	34	40
56	34	40
57	34	40
58	34	40
59	34	40
60	34	40
61	34	40
62	34	40
63	34	40
64	34	40
65	34	40
66	34	40
67	34	40
68	34	40
69	34	40
70	34	40
71	34	40
72	34	40
73	34	40
74	34	39
75	34	40
76	33	39
77	34	40
78	34	40
79	34	40
80	34	39
81	34	39
82	34	39
83	34	39
84	34	39
85	34	39
86	34	39
87	33	38
88	34	39
89	34	39
90	32	38
91	34	39
92	33	39
93	32	38
94	34	39
95	33	39
96	33	38
97	32	38
98	30	36
99	32	37
100	32	37
101	31	36
end
1	32
2	33
3	35
4	34
5	33
6	34
7	35
8	35
9	34
10	34
11	34
12	35
13	35
14	35
15	36
16	36
17	35
18	36
19	36
20	34
21	35
22	36
23	36
24	34
25	35
26	36
27	36
28	36
29	36
30	36
31	35
32	35
33	36
34	34
35	35
36	35
37	35
38	35
39	35
40	35
41	35
42	35
43	35
44	35
45	36
46	35
47	35
48	35
49	35
50	34
51	35
52	36
53	35
54	36
55	36
56	36
57	36
58	35
59	36
60	35
61	35
62	35
63	36
64	35
65	35
66	36
67	35
68	35
69	36
70	35
71	36
72	35
73	36
74	35
75	36
76	35
77	35
78	35
79	35
80	35
81	35
82	35
83	35
84	35
85	35
86	35
87	35
88	35
89	35
90	35
91	35
92	35
93	34
94	35
95	35
96	35
97	34
98	33
99	34
100	34
101	33
end
1	32.60
2	34.20
3	35.81
4	35.51
5	35.44
6	35.71
7	36.39
8	36.66
9	35.74
10	35.54
11	35.55
12	36.90
13	37.00
14	37.04
15	37.04
16	37.06
17	36.81
18	37.16
19	37.20
20	35.38
21	37.04
22	37.25
23	37.32
24	35.58
25	37.12
26	37.27
27	37.26
28	37.14
29	37.26
30	37.18
31	37.27
32	37.18
33	37.26
34	35.55
35	37.18
36	36.29
37	37.14
38	37.18
39	37.16
40	37.14
41	37.11
42	37.09
43	36.76
44	37.06
45	37.12
46	37.14
47	37.13
48	36.25
49	37.16
50	34.77
51	36.81
52	37.14
53	37.20
54	37.25
55	37.28
56	37.27
57	37.30
58	37.21
59	37.42
60	37.24
61	37.28
62	37.23
63	37.31
64	37.30
65	37.35
66	37.22
67	37.23
68	37.21
69	37.27
70	37.19
71	37.19
72	37.20
73	37.21
74	37.08
75	37.25
76	36.39
77	37.03
78	37.13
79	37.14
80	37.16
81	37.07
82	36.91
83	36.85
84	37.02
85	37.01
86	37.03
87	35.92
88	36.95
89	36.87
90	35.73
91	36.88
92	36.85
93	35.08
94	36.73
95	36.73
96	36.40
97	35.63
98	33.97
99	35.92
100	35.66
101	33.93
end
