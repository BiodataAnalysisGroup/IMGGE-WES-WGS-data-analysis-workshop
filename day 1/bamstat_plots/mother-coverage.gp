
            set terminal png size 600,400 truecolor
            set output "./bamstat_plots/mother-coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "mother"
            set xrange [:712]
            plot '-' with lines notitle
        1	213228
2	129749
3	97261
4	82769
5	73105
6	65291
7	60558
8	57032
9	53928
10	51983
11	49769
12	49165
13	46678
14	44696
15	44066
16	43732
17	42457
18	41462
19	41530
20	40181
21	39452
22	38813
23	38888
24	38844
25	39662
26	38733
27	38289
28	37995
29	37831
30	37393
31	37294
32	36211
33	37313
34	36391
35	37070
36	36532
37	36253
38	36431
39	36313
40	37085
41	36201
42	35613
43	35465
44	35353
45	35207
46	35577
47	35545
48	35228
49	35123
50	35383
51	35706
52	35736
53	35302
54	35076
55	34382
56	34013
57	33829
58	33784
59	33956
60	33423
61	33499
62	33584
63	33414
64	33099
65	32565
66	32190
67	31843
68	31234
69	31253
70	30296
71	30687
72	30018
73	29986
74	29129
75	29371
76	28774
77	28488
78	27872
79	27252
80	27234
81	26419
82	26344
83	26291
84	25702
85	25277
86	25188
87	24450
88	24363
89	23939
90	23690
91	23228
92	22850
93	22664
94	22309
95	21608
96	21370
97	21209
98	20436
99	19970
100	20343
101	19736
102	19021
103	18925
104	18807
105	18421
106	18199
107	18027
108	17890
109	17688
110	17184
111	16979
112	16305
113	16073
114	16043
115	15680
116	15232
117	14993
118	15100
119	14960
120	14556
121	14385
122	14037
123	13619
124	13664
125	13010
126	12922
127	12618
128	12247
129	12119
130	11969
131	11950
132	11709
133	11254
134	11060
135	10775
136	10543
137	10476
138	10310
139	9822
140	9531
141	9438
142	8973
143	8818
144	8687
145	8294
146	8255
147	7978
148	7828
149	7494
150	7413
151	7058
152	6992
153	6804
154	6791
155	6504
156	6543
157	6632
158	6265
159	6087
160	5727
161	5779
162	5629
163	5304
164	5137
165	4924
166	4812
167	4757
168	4573
169	4408
170	4258
171	4181
172	4070
173	4100
174	3973
175	3864
176	3857
177	3749
178	3756
179	3558
180	3492
181	3367
182	3311
183	3243
184	3238
185	3098
186	2989
187	2888
188	2829
189	2804
190	2716
191	2575
192	2444
193	2364
194	2389
195	2283
196	2247
197	2211
198	2221
199	2055
200	2069
201	2010
202	1952
203	1845
204	1832
205	1871
206	1758
207	1841
208	1799
209	1616
210	1707
211	1654
212	1705
213	1704
214	1564
215	1651
216	1571
217	1481
218	1391
219	1391
220	1398
221	1383
222	1303
223	1366
224	1225
225	1182
226	1158
227	1148
228	1130
229	1108
230	1098
231	1192
232	1140
233	1129
234	1072
235	1046
236	1035
237	1015
238	950
239	959
240	965
241	896
242	924
243	891
244	858
245	867
246	819
247	842
248	831
249	808
250	736
251	793
252	739
253	722
254	703
255	635
256	691
257	707
258	674
259	660
260	652
261	690
262	642
263	635
264	651
265	620
266	613
267	670
268	544
269	574
270	547
271	590
272	532
273	547
274	541
275	483
276	484
277	478
278	478
279	480
280	465
281	470
282	491
283	447
284	458
285	433
286	418
287	384
288	375
289	410
290	418
291	404
292	379
293	433
294	417
295	396
296	357
297	416
298	426
299	396
300	385
301	376
302	377
303	361
304	324
305	324
306	362
307	351
308	334
309	318
310	339
311	330
312	356
313	355
314	320
315	340
316	309
317	305
318	317
319	322
320	334
321	322
322	321
323	307
324	329
325	325
326	323
327	297
328	319
329	296
330	284
331	277
332	327
333	279
334	278
335	275
336	241
337	265
338	251
339	254
340	252
341	258
342	247
343	270
344	221
345	235
346	218
347	217
348	216
349	243
350	213
351	210
352	196
353	197
354	221
355	162
356	179
357	181
358	185
359	170
360	176
361	189
362	174
363	187
364	178
365	160
366	180
367	179
368	181
369	153
370	162
371	157
372	177
373	161
374	173
375	184
376	178
377	159
378	149
379	190
380	173
381	137
382	179
383	152
384	153
385	147
386	146
387	133
388	163
389	140
390	154
391	166
392	148
393	144
394	166
395	128
396	159
397	168
398	148
399	176
400	149
401	148
402	135
403	144
404	152
405	165
406	143
407	140
408	161
409	145
410	144
411	177
412	174
413	146
414	158
415	132
416	142
417	155
418	123
419	147
420	125
421	162
422	129
423	146
424	162
425	159
426	135
427	152
428	161
429	156
430	158
431	161
432	131
433	164
434	114
435	151
436	151
437	143
438	132
439	155
440	143
441	158
442	130
443	143
444	128
445	148
446	127
447	111
448	134
449	141
450	133
451	131
452	120
453	124
454	117
455	98
456	106
457	107
458	116
459	104
460	112
461	121
462	111
463	112
464	113
465	118
466	103
467	103
468	115
469	120
470	85
471	96
472	95
473	114
474	95
475	83
476	98
477	79
478	82
479	84
480	105
481	93
482	80
483	84
484	92
485	85
486	106
487	103
488	103
489	68
490	89
491	83
492	70
493	88
494	105
495	95
496	92
497	84
498	99
499	86
500	96
501	90
502	70
503	64
504	62
505	55
506	83
507	79
508	74
509	68
510	73
511	79
512	49
513	69
514	59
515	64
516	62
517	72
518	66
519	63
520	61
521	64
522	61
523	48
524	76
525	70
526	57
527	71
528	59
529	69
530	66
531	67
532	66
533	89
534	61
535	67
536	65
537	77
538	63
539	73
540	72
541	79
542	70
543	76
544	73
545	76
546	65
547	70
548	71
549	51
550	53
551	64
552	54
553	61
554	69
555	68
556	54
557	52
558	66
559	54
560	52
561	60
562	51
563	53
564	50
565	56
566	57
567	45
568	70
569	61
570	71
571	58
572	47
573	40
574	59
575	61
576	46
577	65
578	37
579	42
580	47
581	47
582	59
583	56
584	62
585	57
586	63
587	53
588	67
589	53
590	54
591	48
592	55
593	57
594	46
595	53
596	53
597	55
598	48
599	40
600	59
601	64
602	54
603	54
604	52
605	62
606	49
607	42
608	52
609	46
610	44
611	50
612	46
613	64
614	54
615	46
616	52
617	47
618	54
619	49
620	60
621	49
622	55
623	46
624	39
625	37
626	53
627	45
628	51
629	45
630	40
631	39
632	33
633	27
634	34
635	29
636	32
637	33
638	26
639	38
640	35
641	45
642	28
643	28
644	25
645	44
646	48
647	33
648	34
649	49
650	50
651	30
652	41
653	42
654	34
655	40
656	26
657	46
658	41
659	38
660	53
661	45
662	39
663	36
664	36
665	39
666	43
667	37
668	27
669	38
670	41
671	29
672	28
673	31
674	37
675	36
676	31
677	32
678	26
679	41
680	24
681	44
682	42
683	27
684	31
685	28
686	29
687	29
688	32
689	31
690	35
691	25
692	24
693	25
694	20
695	28
696	31
697	33
698	23
699	21
700	26
701	28
702	24
703	20
704	32
705	42
706	22
707	30
708	36
709	34
710	27
711	25
712	19
713	30
714	31
715	33
716	37
717	21
718	31
719	31
720	28
721	37
722	37
723	37
724	46
725	49
726	51
727	48
728	37
729	40
730	52
731	54
732	50
733	43
734	37
735	45
736	50
737	63
738	40
739	36
740	32
741	39
742	44
743	42
744	34
745	49
746	44
747	57
748	34
749	33
750	40
751	46
752	42
753	43
754	52
755	37
756	53
757	35
758	56
759	48
760	44
761	45
762	40
763	38
764	52
765	41
766	49
767	39
768	40
769	33
770	33
771	27
772	33
773	29
774	40
775	27
776	34
777	20
778	38
779	23
780	31
781	24
782	21
783	17
784	24
785	28
786	23
787	19
788	23
789	26
790	25
791	23
792	24
793	34
794	16
795	18
796	23
797	21
798	23
799	25
800	17
801	16
802	12
803	18
804	19
805	11
806	18
807	20
808	19
809	14
810	9
811	16
812	14
813	13
814	16
815	24
816	17
817	26
818	15
819	16
820	17
821	21
822	12
823	17
824	8
825	13
826	14
827	10
828	17
829	15
830	14
831	10
832	14
833	15
834	17
835	10
836	15
837	17
838	14
839	12
840	18
841	8
842	7
843	20
844	16
845	15
846	11
847	11
848	16
849	13
850	17
851	15
852	9
853	10
854	21
855	20
856	9
857	10
858	11
859	10
860	9
861	17
862	20
863	22
864	19
865	10
866	9
867	20
868	19
869	18
870	15
871	17
872	14
873	13
874	12
875	9
876	11
877	13
878	9
879	12
880	11
881	17
882	8
883	21
884	13
885	12
886	8
887	18
888	8
889	11
890	9
891	13
892	5
893	17
894	8
895	9
896	13
897	12
898	12
899	11
900	8
901	15
902	8
903	16
904	12
905	13
906	8
907	8
908	9
909	10
910	16
911	14
912	7
913	12
914	16
915	8
916	21
917	19
918	6
919	9
920	11
921	13
922	10
923	14
924	12
925	17
926	10
927	9
928	10
929	9
930	11
931	10
932	11
933	7
934	13
935	9
936	9
937	9
938	3
939	4
940	15
941	10
942	15
943	9
944	11
945	8
946	18
947	14
948	13
949	16
950	12
951	19
952	16
953	6
954	13
955	15
956	11
957	8
958	14
959	10
960	5
961	3
962	8
963	5
964	13
965	5
966	9
967	7
968	10
969	6
970	6
971	6
972	6
973	11
974	5
975	7
976	9
977	4
978	5
979	4
980	10
981	7
982	3
983	12
984	10
985	9
986	8
987	12
988	5
989	5
990	7
991	7
992	10
993	4
994	8
995	7
996	6
997	6
998	7
999	8
1000	5
1000	4219
end