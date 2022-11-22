
            set terminal png size 600,400 truecolor
            set output "./bamstat_plots/father-gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "father"
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",38.94) at 38.94,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' , '-' smooth csplines with lines lc 2 title 'Last fragments'
        1	0.000100
3	0.000083
4	0.000100
5	0.000183
5	0.000249
6	0.000332
6	0.000415
7	0.000482
7	0.000498
8	0.000548
8	0.000598
9	0.000581
9	0.000631
10	0.000781
10	0.000997
11	0.001279
11	0.001578
12	0.001694
12	0.002126
13	0.002359
13	0.003870
14	0.004070
14	0.005482
15	0.005997
15	0.009518
16	0.010249
16	0.016180
17	0.016994
17	0.025565
18	0.026313
18	0.040150
19	0.041363
19	0.063107
20	0.064353
20	0.091480
21	0.093374
21	0.131896
22	0.133474
22	0.184388
23	0.186016
23	0.250835
24	0.254124
24	0.319690
25	0.323095
25	0.391717
26	0.393894
26	0.483231
27	0.486403
27	0.567318
28	0.570458
28	0.649795
29	0.650360
29	0.718534
30	0.723168
30	0.791957
31	0.792156
31	0.842257
32	0.844981
32	0.890330
33	0.921261
33	0.926577
34	0.945016
34	0.948155
35	0.960249
35	0.960730
36	0.983471
36	0.986412
37	0.990880
37	0.990598
38	0.994767
38	1.000000
39	0.999120
39	0.996894
40	0.974319
40	0.976196
41	0.948089
41	0.949069
42	0.896742
42	0.895879
43	0.869931
43	0.872622
44	0.841692
44	0.840429
45	0.804233
45	0.805993
46	0.783535
46	0.784864
47	0.738185
47	0.737521
48	0.691556
48	0.690958
49	0.659081
50	0.660094
50	0.616605
51	0.612419
51	0.600575
52	0.604993
52	0.598963
53	0.599146
53	0.580342
54	0.582634
54	0.566820
55	0.569710
55	0.550142
56	0.551936
56	0.543281
57	0.549311
57	0.545474
58	0.549145
58	0.575724
59	0.585525
59	0.612103
60	0.615608
60	0.586970
61	0.585957
61	0.563182
62	0.568681
62	0.564229
63	0.570259
63	0.531604
64	0.519294
64	0.461619
65	0.447964
65	0.409143
66	0.401319
66	0.350903
67	0.316417
67	0.309141
68	0.285835
68	0.279357
69	0.228957
69	0.222545
70	0.184787
70	0.177893
71	0.144620
71	0.135932
72	0.104586
72	0.099321
73	0.078224
73	0.072044
74	0.057160
74	0.054137
75	0.044569
75	0.040864
76	0.032293
76	0.028954
77	0.022276
77	0.020648
78	0.016828
78	0.015316
79	0.012575
79	0.011146
80	0.008273
80	0.007376
81	0.005880
81	0.005033
82	0.004153
82	0.003522
83	0.002425
83	0.002193
84	0.001512
84	0.001246
85	0.001096
85	0.001047
86	0.000880
86	0.000714
87	0.000548
87	0.000465
88	0.000365
89	0.000249
89	0.000166
90	0.000100
90	0.000083
91	0.000066
91	0.000050
92	0.000066
93	0.000050
94	0.000033
end
1	0.000084
2	0.000101
3	0.000051
4	0.000034
4	0.000084
5	0.000236
6	0.000270
6	0.000354
7	0.000388
8	0.000421
8	0.000691
9	0.000792
10	0.000843
10	0.001130
11	0.001399
11	0.001534
12	0.001804
12	0.003035
13	0.003372
13	0.004822
14	0.005041
14	0.006305
15	0.007115
15	0.012105
16	0.012577
16	0.018680
17	0.019489
17	0.030819
18	0.031595
18	0.048454
19	0.050224
19	0.074822
20	0.075092
20	0.104663
21	0.106248
21	0.149948
22	0.152426
22	0.212125
23	0.214401
23	0.279546
24	0.282463
24	0.351148
25	0.354722
25	0.425667
26	0.428786
26	0.513488
27	0.514803
27	0.606434
28	0.609772
28	0.677749
29	0.681896
29	0.751678
30	0.756364
30	0.826989
31	0.830681
31	0.863135
32	0.865951
32	0.919901
33	0.947314
33	0.952086
34	0.969097
34	0.972131
35	0.975992
35	0.977644
36	0.994268
36	0.998263
37	0.989884
37	0.990660
38	0.991435
38	0.996544
39	0.997960
39	1.000000
40	0.972637
40	0.971103
41	0.945325
41	0.947955
42	0.899012
42	0.900782
43	0.873217
43	0.875038
44	0.838672
44	0.836666
45	0.797839
45	0.800098
46	0.779024
46	0.781856
47	0.735324
47	0.737111
48	0.701032
48	0.700155
49	0.665020
50	0.662997
50	0.616549
51	0.616313
51	0.611053
52	0.615740
52	0.604562
53	0.605152
53	0.596031
54	0.598459
54	0.583538
55	0.581785
55	0.563577
56	0.569545
56	0.563223
57	0.570152
57	0.568045
58	0.572141
58	0.590940
59	0.599504
59	0.623462
60	0.628098
60	0.610530
61	0.614020
61	0.582780
62	0.584769
62	0.575783
63	0.581128
63	0.555130
64	0.546127
64	0.498786
65	0.486934
65	0.442256
66	0.432798
66	0.389959
67	0.351772
67	0.340358
68	0.312034
68	0.306572
69	0.257191
69	0.248441
70	0.207354
70	0.196884
71	0.159878
71	0.150504
72	0.121944
72	0.113818
73	0.088967
73	0.083286
74	0.069714
74	0.064572
75	0.051876
75	0.047813
76	0.038659
76	0.036231
77	0.030684
77	0.028020
78	0.022895
78	0.020956
79	0.017466
79	0.015443
80	0.012358
80	0.011077
81	0.008969
81	0.008396
82	0.006997
82	0.006069
83	0.004164
83	0.003372
84	0.002461
84	0.001871
85	0.001568
85	0.001416
86	0.001113
86	0.000944
87	0.000843
87	0.000607
88	0.000489
88	0.000337
89	0.000287
89	0.000236
90	0.000202
90	0.000118
91	0.000084
92	0.000051
92	0.000034
93	0.000017
end