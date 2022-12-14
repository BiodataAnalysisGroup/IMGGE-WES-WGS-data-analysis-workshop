
            set terminal png size 600,400 truecolor
            set output "./bamstat_plots/father-insert-size.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set rmargin 5
            set label sprintf("%d",253) at 253+10,31957
            set ylabel  "Number of pairs"
            set xlabel  "Insert Size"
            set title "father"
            plot \
                '-' with lines lc rgb 'black' title 'All pairs', \
                '-' with lines title 'Inward', \
                '-' with lines title 'Outward', \
                '-' with lines title 'Other'
        0	5
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	2
19	0
20	0
21	1
22	1
23	0
24	1
25	0
26	0
27	0
28	0
29	0
30	0
31	0
32	3
33	0
34	1
35	0
36	0
37	0
38	0
39	0
40	1
41	0
42	1
43	0
44	1
45	0
46	0
47	0
48	0
49	0
50	0
51	0
52	0
53	0
54	0
55	0
56	2
57	0
58	0
59	0
60	0
61	0
62	0
63	0
64	0
65	0
66	1
67	0
68	0
69	0
70	0
71	0
72	0
73	3
74	0
75	0
76	1
77	0
78	0
79	0
80	0
81	2
82	1
83	1
84	0
85	0
86	0
87	1
88	3
89	0
90	1
91	0
92	3
93	1
94	10
95	6
96	8
97	3
98	6
99	17
100	24
101	2
102	10
103	27
104	10
105	24
106	35
107	30
108	26
109	10
110	25
111	29
112	11
113	24
114	45
115	31
116	44
117	67
118	37
119	40
120	61
121	47
122	49
123	68
124	74
125	37
126	43
127	81
128	59
129	68
130	67
131	94
132	71
133	85
134	48
135	89
136	70
137	109
138	112
139	90
140	104
141	85
142	118
143	113
144	113
145	108
146	124
147	176
148	137
149	202
150	199
151	186
152	159
153	177
154	208
155	198
156	195
157	219
158	283
159	253
160	244
161	263
162	244
163	314
164	334
165	276
166	318
167	328
168	366
169	382
170	397
171	362
172	460
173	390
174	447
175	467
176	473
177	526
178	650
179	563
180	593
181	611
182	708
183	719
184	759
185	817
186	901
187	898
188	1054
189	1013
190	1090
191	1134
192	1179
193	1359
194	1369
195	1533
196	1622
197	1708
198	1864
199	1920
200	2234
201	2380
202	2505
203	2759
204	2929
205	3164
206	3620
207	3745
208	4197
209	4404
210	5071
211	5326
212	5666
213	6100
214	6856
215	7120
216	7925
217	8580
218	9242
219	9996
220	10276
221	11071
222	11948
223	13034
224	13932
225	14518
226	15291
227	16310
228	17107
229	17774
230	18669
231	19383
232	20285
233	21249
234	22184
235	22733
236	23851
237	24517
238	25247
239	26119
240	26840
241	27419
242	28223
243	28417
244	29337
245	29476
246	29800
247	30778
248	30806
249	30881
250	31673
251	31418
252	31502
253	31957
254	31796
255	31718
256	31351
257	31345
258	30319
259	30605
260	30284
261	30343
262	29643
263	29263
264	28812
265	28649
266	27674
267	27301
268	26687
269	25735
270	25700
271	24842
272	23846
273	23206
274	22373
275	21782
276	20796
277	20311
278	19655
279	18904
280	18088
281	17167
282	16334
283	15320
284	14838
285	13856
286	13699
287	12774
288	11584
289	10746
290	10280
291	9399
292	8739
293	8239
294	7658
295	7112
296	6450
297	5934
298	5357
299	4859
300	4243
301	4008
302	3717
303	3078
304	2900
305	2583
306	2268
307	1979
end
0	0
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
21	0
22	0
23	0
24	0
25	0
26	0
27	0
28	0
29	0
30	0
31	0
32	0
33	0
34	0
35	0
36	0
37	0
38	0
39	0
40	0
41	0
42	0
43	0
44	0
45	0
46	0
47	0
48	0
49	0
50	0
51	0
52	0
53	0
54	0
55	0
56	0
57	0
58	0
59	0
60	0
61	0
62	0
63	0
64	0
65	0
66	0
67	0
68	0
69	0
70	0
71	0
72	0
73	0
74	0
75	0
76	0
77	0
78	0
79	0
80	0
81	0
82	0
83	1
84	0
85	0
86	0
87	0
88	3
89	0
90	1
91	0
92	3
93	1
94	6
95	6
96	7
97	3
98	6
99	15
100	7
101	2
102	10
103	27
104	10
105	24
106	35
107	30
108	26
109	10
110	25
111	27
112	11
113	24
114	42
115	31
116	43
117	65
118	35
119	40
120	61
121	47
122	49
123	67
124	74
125	37
126	43
127	81
128	59
129	67
130	67
131	94
132	71
133	85
134	48
135	86
136	70
137	109
138	112
139	90
140	104
141	85
142	118
143	113
144	112
145	108
146	124
147	174
148	136
149	202
150	199
151	185
152	159
153	176
154	208
155	198
156	195
157	219
158	283
159	253
160	244
161	263
162	244
163	314
164	334
165	276
166	318
167	328
168	365
169	380
170	397
171	362
172	460
173	390
174	447
175	467
176	472
177	526
178	650
179	561
180	593
181	611
182	708
183	719
184	759
185	817
186	901
187	898
188	1054
189	1013
190	1090
191	1134
192	1179
193	1359
194	1367
195	1531
196	1622
197	1707
198	1864
199	1920
200	2234
201	2380
202	2505
203	2759
204	2929
205	3164
206	3620
207	3741
208	4196
209	4404
210	5071
211	5326
212	5665
213	6100
214	6856
215	7120
216	7924
217	8579
218	9242
219	9996
220	10276
221	11071
222	11948
223	13034
224	13932
225	14517
226	15291
227	16310
228	17107
229	17773
230	18669
231	19381
232	20283
233	21249
234	22184
235	22733
236	23851
237	24517
238	25247
239	26119
240	26840
241	27419
242	28223
243	28417
244	29337
245	29476
246	29800
247	30778
248	30806
249	30881
250	31673
251	31418
252	31502
253	31957
254	31796
255	31718
256	31351
257	31345
258	30319
259	30605
260	30284
261	30343
262	29643
263	29263
264	28811
265	28646
266	27674
267	27301
268	26682
269	25735
270	25700
271	24842
272	23846
273	23206
274	22373
275	21782
276	20796
277	20311
278	19655
279	18904
280	18088
281	17167
282	16334
283	15320
284	14838
285	13856
286	13699
287	12774
288	11584
289	10746
290	10280
291	9399
292	8738
293	8239
294	7658
295	7112
296	6450
297	5934
298	5357
299	4859
300	4243
301	4008
302	3717
303	3078
304	2900
305	2583
306	2268
307	1979
end
0	0
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	0
19	0
20	0
21	0
22	0
23	0
24	0
25	0
26	0
27	0
28	0
29	0
30	0
31	0
32	0
33	0
34	0
35	0
36	0
37	0
38	0
39	0
40	0
41	0
42	0
43	0
44	1
45	0
46	0
47	0
48	0
49	0
50	0
51	0
52	0
53	0
54	0
55	0
56	0
57	0
58	0
59	0
60	0
61	0
62	0
63	0
64	0
65	0
66	0
67	0
68	0
69	0
70	0
71	0
72	0
73	0
74	0
75	0
76	1
77	0
78	0
79	0
80	0
81	0
82	1
83	0
84	0
85	0
86	0
87	0
88	0
89	0
90	0
91	0
92	0
93	0
94	3
95	0
96	1
97	0
98	0
99	2
100	17
101	0
102	0
103	0
104	0
105	0
106	0
107	0
108	0
109	0
110	0
111	0
112	0
113	0
114	0
115	0
116	0
117	0
118	0
119	0
120	0
121	0
122	0
123	0
124	0
125	0
126	0
127	0
128	0
129	0
130	0
131	0
132	0
133	0
134	0
135	0
136	0
137	0
138	0
139	0
140	0
141	0
142	0
143	0
144	0
145	0
146	0
147	0
148	0
149	0
150	0
151	0
152	0
153	0
154	0
155	0
156	0
157	0
158	0
159	0
160	0
161	0
162	0
163	0
164	0
165	0
166	0
167	0
168	0
169	0
170	0
171	0
172	0
173	0
174	0
175	0
176	0
177	0
178	0
179	0
180	0
181	0
182	0
183	0
184	0
185	0
186	0
187	0
188	0
189	0
190	0
191	0
192	0
193	0
194	0
195	0
196	0
197	0
198	0
199	0
200	0
201	0
202	0
203	0
204	0
205	0
206	0
207	0
208	0
209	0
210	0
211	0
212	0
213	0
214	0
215	0
216	0
217	0
218	0
219	0
220	0
221	0
222	0
223	0
224	0
225	0
226	0
227	0
228	0
229	0
230	0
231	0
232	0
233	0
234	0
235	0
236	0
237	0
238	0
239	0
240	0
241	0
242	0
243	0
244	0
245	0
246	0
247	0
248	0
249	0
250	0
251	0
252	0
253	0
254	0
255	0
256	0
257	0
258	0
259	0
260	0
261	0
262	0
263	0
264	0
265	0
266	0
267	0
268	0
269	0
270	0
271	0
272	0
273	0
274	0
275	0
276	0
277	0
278	0
279	0
280	0
281	0
282	0
283	0
284	0
285	0
286	0
287	0
288	0
289	0
290	0
291	0
292	0
293	0
294	0
295	0
296	0
297	0
298	0
299	0
300	0
301	0
302	0
303	0
304	0
305	0
306	0
307	0
end
0	5
1	0
2	0
3	0
4	0
5	0
6	0
7	0
8	0
9	0
10	0
11	0
12	0
13	0
14	0
15	0
16	0
17	0
18	2
19	0
20	0
21	1
22	1
23	0
24	1
25	0
26	0
27	0
28	0
29	0
30	0
31	0
32	3
33	0
34	1
35	0
36	0
37	0
38	0
39	0
40	1
41	0
42	1
43	0
44	0
45	0
46	0
47	0
48	0
49	0
50	0
51	0
52	0
53	0
54	0
55	0
56	2
57	0
58	0
59	0
60	0
61	0
62	0
63	0
64	0
65	0
66	1
67	0
68	0
69	0
70	0
71	0
72	0
73	3
74	0
75	0
76	0
77	0
78	0
79	0
80	0
81	2
82	0
83	0
84	0
85	0
86	0
87	1
88	0
89	0
90	0
91	0
92	0
93	0
94	1
95	0
96	0
97	0
98	0
99	0
100	0
101	0
102	0
103	0
104	0
105	0
106	0
107	0
108	0
109	0
110	0
111	2
112	0
113	0
114	3
115	0
116	1
117	2
118	2
119	0
120	0
121	0
122	0
123	1
124	0
125	0
126	0
127	0
128	0
129	1
130	0
131	0
132	0
133	0
134	0
135	3
136	0
137	0
138	0
139	0
140	0
141	0
142	0
143	0
144	1
145	0
146	0
147	2
148	1
149	0
150	0
151	1
152	0
153	1
154	0
155	0
156	0
157	0
158	0
159	0
160	0
161	0
162	0
163	0
164	0
165	0
166	0
167	0
168	1
169	2
170	0
171	0
172	0
173	0
174	0
175	0
176	1
177	0
178	0
179	2
180	0
181	0
182	0
183	0
184	0
185	0
186	0
187	0
188	0
189	0
190	0
191	0
192	0
193	0
194	2
195	2
196	0
197	1
198	0
199	0
200	0
201	0
202	0
203	0
204	0
205	0
206	0
207	4
208	1
209	0
210	0
211	0
212	1
213	0
214	0
215	0
216	1
217	1
218	0
219	0
220	0
221	0
222	0
223	0
224	0
225	1
226	0
227	0
228	0
229	1
230	0
231	2
232	2
233	0
234	0
235	0
236	0
237	0
238	0
239	0
240	0
241	0
242	0
243	0
244	0
245	0
246	0
247	0
248	0
249	0
250	0
251	0
252	0
253	0
254	0
255	0
256	0
257	0
258	0
259	0
260	0
261	0
262	0
263	0
264	1
265	3
266	0
267	0
268	5
269	0
270	0
271	0
272	0
273	0
274	0
275	0
276	0
277	0
278	0
279	0
280	0
281	0
282	0
283	0
284	0
285	0
286	0
287	0
288	0
289	0
290	0
291	0
292	1
293	0
294	0
295	0
296	0
297	0
298	0
299	0
300	0
301	0
302	0
303	0
304	0
305	0
306	0
307	0
end
