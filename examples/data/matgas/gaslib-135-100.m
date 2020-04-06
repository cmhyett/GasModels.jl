function mgc = gaslib-135-100

%% required global data
mgc.gas_specific_gravity         = 0.6;
mgc.specific_heat_capacity_ratio = 1.4;  % unitless
mgc.temperature                  = 273.15;  % K
mgc.compressibility_factor       = 0.8;  % unitless
mgc.units                        = 'si';

%% optional global data (that was either provided or computed based on required global data)
mgc.gas_molar_mass              = 0.01857; % kg/mol
mgc.R                           = 8.314;  % J/(mol K)
mgc.base_pressure               = 8101325;  % Pa
mgc.base_flow                   = 2200.0;
mgc.base_length                 = 5000;  % m
mgc.is_per_unit                 = 0;
mgc.sound_speed                 = 312.8060

%% junction data
% id	p_min	p_max	p_nominal	junction_type	status	pipeline_name	edi_id	lat	lon
mgc.junction = [
0	      101325	8101325	101325	0	1	'gaslib-135'	0	      49.9631	6.3113
1	      101325	8101325	101325	0	1	'gaslib-135'	1	      51.5156	17.2343
2	      101325	8101325	101325	0	1	'gaslib-135'	2	      52.1648	19.4013
3	      101325	8101325	101325	0	1	'gaslib-135'	3	      48.4385	7.2498
4	      101325	8101325	101325	0	1	'gaslib-135'	4	      51.2692	16.9216
5	      3101325	8101325	3101325	0	1	'gaslib-135'	5	      49.4445	8.8146
6	      101325	8101325	101325	0	1	'gaslib-135'	6	      52.5495	16.6284
7	      101325	8101325	101325	0	1	'gaslib-135'	7	      54.7328	18.1588
8	      101325	8101325	101325	0	1	'gaslib-135'	8	      52.7131	10.6960
9	      101325	8101325	101325	0	1	'gaslib-135'	9	      50.3393	6.4216
10	    101325	8101325	101325	0	1	'gaslib-135'	10	    54.3424	12.3038
11	    101325	8101325	101325	0	1	'gaslib-135'	11	    50.2538	6.7533
12	    101325	8101325	101325	0	1	'gaslib-135'	12	    51.9783	11.0518
13	    101325	8101325	101325	0	1	'gaslib-135'	13	    53.7471	9.8835
14	    101325	8101325	101325	0	1	'gaslib-135'	14	    52.9274	9.8013
15	    101325	8101325	101325	0	1	'gaslib-135'	15	    49.5908	14.0671
16	    101325	8101325	101325	0	1	'gaslib-135'	16	    54.8084	18.1004
17	    3101325	8101325	3101325	0	1	'gaslib-135'	17	    51.0905	16.5184
18	    101325	8101325	101325	0	1	'gaslib-135'	18	    49.8335	5.9731
19	    101325	8101325	101325	0	1	'gaslib-135'	19	    51.1473	13.8830
20	    101325	8101325	101325	0	1	'gaslib-135'	20	    51.6988	15.2097
21	    101325	8101325	101325	0	1	'gaslib-135'	21	    51.5947	14.1580
22	    101325	8101325	101325	0	1	'gaslib-135'	22	    51.0828	7.6024
23	    101325	8101325	101325	0	1	'gaslib-135'	23	    51.2520	16.7159
24	    101325	8101325	101325	0	1	'gaslib-135'	24	    50.8516	15.1007
25	    101325	8101325	101325	0	1	'gaslib-135'	25	    51.7403	16.7419
26	    101325	8101325	101325	0	1	'gaslib-135'	26	    50.1445	6.9712
27	    101325	8101325	101325	0	1	'gaslib-135'	27	    52.2903	18.6940
28	    101325	8101325	101325	0	1	'gaslib-135'	28	    49.7007	7.4678
29	    101325	8101325	101325	0	1	'gaslib-135'	29	    54.6856	16.9072
30	    3101325	8101325	3101325	0	1	'gaslib-135'	30	    49.6554	11.6318
31	    101325	8101325	101325	0	1	'gaslib-135'	31	    49.7056	11.4570
32	    101325	8101325	101325	0	1	'gaslib-135'	32	    50.6319	8.4142
33	    101325	8101325	101325	0	1	'gaslib-135'	33	    49.9115	10.8765
34	    3101325	8101325	3101325	0	1	'gaslib-135'	34	    47.2275	12.0628
35	    101325	8101325	101325	0	1	'gaslib-135'	35	    50.9196	7.4417
36	    101325	8101325	101325	0	1	'gaslib-135'	36	    55.1028	18.2162
37	    101325	8101325	101325	0	1	'gaslib-135'	37	    52.3959	18.0084
38	    101325	8101325	101325	0	1	'gaslib-135'	38	    50.9866	11.4186
39	    101325	8101325	101325	0	1	'gaslib-135'	39	    49.5533	11.5163
40	    101325	8101325	101325	0	1	'gaslib-135'	40	    0.0     0.0
41	    101325	8101325	101325	0	1	'gaslib-135'	41	    0.0     0.0
42	    101325	8101325	101325	0	1	'gaslib-135'	42	    0.0     0.0
43	    101325	8101325	101325	0	1	'gaslib-135'	43	    0.0     0.0
44	    3101325	8101325	3101325	0	1	'gaslib-135'	44	    0.0     0.0
45	    101325	8101325	101325	0	1	'gaslib-135'	45	    0.0     0.0
46	    101325	8101325	101325	0	1	'gaslib-135'	46	    0.0     0.0
47	    101325	8101325	101325	0	1	'gaslib-135'	47	    0.0     0.0
48	    101325	8101325	101325	0	1	'gaslib-135'	48	    0.0     0.0
49	    3101325	8101325	3101325	0	1	'gaslib-135'	49	    0.0     0.0
50	    101325	8101325	101325	0	1	'gaslib-135'	50	    0.0     0.0
51	    101325	8101325	101325	0	1	'gaslib-135'	51	    0.0     0.0
52	    101325	8101325	101325	0	1	'gaslib-135'	52	    0.0     0.0
53	    101325	8101325	101325	0	1	'gaslib-135'	53	    0.0     0.0
54	    101325	8101325	101325	0	1	'gaslib-135'	54	    0.0     0.0
55	    101325	8101325	101325	0	1	'gaslib-135'	55	    0.0     0.0
56	    3101325	8101325	3101325	0	1	'gaslib-135'	56	    0.0     0.0
57	    3101325	8101325	3101325	0	1	'gaslib-135'	57	    0.0     0.0
58	    101325	8101325	101325	0	1	'gaslib-135'	58	    0.0     0.0
59	    3101325	8101325	3101325	0	1	'gaslib-135'	59	    0.0     0.0
60	    101325	8101325	101325	0	1	'gaslib-135'	60	    0.0     0.0
61	    101325	8101325	101325	0	1	'gaslib-135'	61	    0.0     0.0
62	    101325	8101325	101325	0	1	'gaslib-135'	62	    0.0     0.0
63	    3101325	8101325	3101325	0	1	'gaslib-135'	63	    0.0     0.0
64	    101325	8101325	101325	0	1	'gaslib-135'	64	    0.0     0.0
65	    101325	8101325	101325	0	1	'gaslib-135'	65	    0.0     0.0
66	    3101325	8101325	3101325	0	1	'gaslib-135'	66	    0.0     0.0
67	    101325	8101325	101325	0	1	'gaslib-135'	67	    0.0     0.0
68	    101325	8101325	101325	0	1	'gaslib-135'	68	    0.0     0.0
69	    101325	8101325	101325	0	1	'gaslib-135'	69	    0.0     0.0
70	    101325	8101325	101325	0	1	'gaslib-135'	70	    0.0     0.0
71	    3101325	8101325	3101325	0	1	'gaslib-135'	71	    0.0     0.0
72	    101325	8101325	101325	0	1	'gaslib-135'	72	    0.0     0.0
73    	101325	8101325	101325	0	1	'gaslib-135'	73	    0.0     0.0
74	    101325	8101325	101325	0	1	'gaslib-135'	74	    0.0     0.0
75	    101325	8101325	101325	0	1	'gaslib-135'	75	    0.0     0.0
76	    3101325	8101325	3101325	0	1	'gaslib-135'	76	    0.0     0.0
77	    101325	8101325	101325	0	1	'gaslib-135'	77	    0.0     0.0
78	    101325	8101325	101325	0	1	'gaslib-135'	78	    0.0     0.0
79	    101325	8101325	101325	0	1	'gaslib-135'	79	    0.0     0.0
80	    101325	8101325	101325	0	1	'gaslib-135'	80	    0.0     0.0
81	    101325	8101325	101325	0	1	'gaslib-135'	81	    0.0     0.0
82	    101325	8101325	101325	0	1	'gaslib-135'	82	    0.0     0.0
83	    101325	8101325	101325	0	1	'gaslib-135'	83	    0.0     0.0
84	    101325	8101325	101325	0	1	'gaslib-135'	84	    0.0     0.0
85	    101325	8101325	101325	0	1	'gaslib-135'	85	    0.0     0.0
86	    101325	8101325	101325	0	1	'gaslib-135'	86	    0.0     0.0
87	    101325	8101325	101325	0	1	'gaslib-135'	87	    0.0     0.0
88	    101325	8101325	101325	0	1	'gaslib-135'	88	    0.0     0.0
89	    101325	8101325	101325	0	1	'gaslib-135'	89	    0.0     0.0
90	    101325	8101325	101325	0	1	'gaslib-135'	90	    0.0     0.0
91     	101325	8101325	101325	0	1	'gaslib-135'	91	    0.0     0.0
92	    101325	8101325	101325	0	1	'gaslib-135'	92	    0.0     0.0
93	    101325	8101325	101325	0	1	'gaslib-135'	93	    0.0     0.0
94	    101325	8101325	101325	0	1	'gaslib-135'	94	    0.0     0.0
95	    101325	8101325	101325	0	1	'gaslib-135'	95	    0.0     0.0
96	    101325	8101325	101325	0	1	'gaslib-135'	96	    0.0     0.0
97	    101325	8101325	101325	0	1	'gaslib-135'	97	    0.0     0.0
98	    101325	8101325	101325	0	1	'gaslib-135'	98	    0.0     0.0
99	    101325	8101325	101325	0	1	'gaslib-135'	99	    0.0     0.0
100	    101325	8101325	101325	0	1	'gaslib-135'	100   	0.0     0.0
101	    101325	8101325	101325	0	1	'gaslib-135'	101	    0.0     0.0
102	    101325	8101325	101325	0	1	'gaslib-135'	102	    0.0     0.0
103	    101325	8101325	101325	0	1	'gaslib-135'	103	    0.0     0.0
104	    101325	8101325	101325	0	1	'gaslib-135'	104	    0.0     0.0
105	    101325	7101325	101325	0	1	'gaslib-135'	105	    0.0     0.0
106	    101325	7101325	101325	0	1	'gaslib-135'	106	    0.0     0.0
107	    101325	7101325	101325	0	1	'gaslib-135'	107	    0.0     0.0
108	    101325	7101325	101325	0	1	'gaslib-135'	108	    0.0     0.0
109	    101325	7101325	101325	0	1	'gaslib-135'	109	    0.0     0.0
110	    101325	7101325	101325	0	1	'gaslib-135'	110	    0.0     0.0
111	    101325	7101325	101325	0	1	'gaslib-135'	111	    0.0     0.0
112	    101325	7101325	101325	0	1	'gaslib-135'	112	    0.0     0.0
113	    3101325	8101325	3101325	0	1	'gaslib-135'	113	    0.0     0.0
114	    3101325	8101325	3101325	0	1	'gaslib-135'	114	    0.0     0.0
115	    101325	7101325	101325	0	1	'gaslib-135'	115	    0.0     0.0
116	    101325	7101325	101325	0	1	'gaslib-135'	116	    0.0     0.0
117	    101325	7101325	101325	0	1	'gaslib-135'	117	    0.0     0.0
118	    101325	7101325	101325	0	1	'gaslib-135'	118	    0.0     0.0
119	    101325	7101325	101325	0	1	'gaslib-135'	119	    0.0     0.0
120	    101325	7101325	101325	0	1	'gaslib-135'	120	    0.0     0.0
121	    101325	7101325	101325	0	1	'gaslib-135'	121	    0.0     0.0
122	    101325	7101325	101325	0	1	'gaslib-135'	122	    0.0     0.0
123	    101325	7101325	101325	0	1	'gaslib-135'	123	    0.0     0.0
124	    101325	7101325	101325	0	1	'gaslib-135'	124	    0.0     0.0
125	    101325	7101325	101325	0	1	'gaslib-135'	125	    0.0     0.0
126	    101325	7101325	101325	0	1	'gaslib-135'	126	    0.0     0.0
127	    101325	7101325	101325	0	1	'gaslib-135'	127	    0.0     0.0
128	    101325	7101325	101325	0	1	'gaslib-135'	128	    0.0     0.0
129	    101325	7101325	101325	0	1	'gaslib-135'	129	    0.0     0.0
130	    101325	7101325	101325	0	1	'gaslib-135'	130	    0.0     0.0
131	    101325	7101325	101325	0	1	'gaslib-135'	131	    0.0     0.0
132	    101325	7101325	101325	0	1	'gaslib-135'	132	    0.0     0.0
133	    101325	7101325	101325	0	1	'gaslib-135'	133	    0.0     0.0
134	    101325	7101325	101325	0	1	'gaslib-135'	134	    0.0     0.0
100057	101325	8101325	101325	0	1	'gaslib-135'	100057	0.0     0.0
200044	101325	8101325	101325	0	1	'gaslib-135'	200044	0.0     0.0
300071	101325	8101325	101325	0	1	'gaslib-135'	300071	0.0     0.0
400044	101325	8101325	101325	0	1	'gaslib-135'	400044	0.0     0.0
500017	101325	8101325	101325	0	1	'gaslib-135'	500017	0.0     0.0
600030	101325	8101325	101325	0	1	'gaslib-135'	600030	0.0     0.0
700044	101325	8101325	101325	0	1	'gaslib-135'	700044	0.0     0.0
800057	101325	8101325	101325	0	1	'gaslib-135'	800057	0.0     0.0
900005	101325	8101325	101325	0	1	'gaslib-135'	900005	0.0     0.0
1000059	101325	8101325	101325	0	1	'gaslib-135'	1000059	0.0     0.0
1100066	101325	8101325	101325	0	1	'gaslib-135'	1100066	0.0     0.0
1200071	101325	8101325	101325	0	1	'gaslib-135'	1200071	0.0     0.0
1300005	101325	8101325	101325	0	1	'gaslib-135'	1300005	0.0     0.0
1400030	101325	8101325	101325	0	1	'gaslib-135'	1400030	0.0     0.0
1500044	101325	8101325	101325	0	1	'gaslib-135'	1500044	0.0     0.0
1600113	101325	8101325	101325	0	1	'gaslib-135'	1600113	0.0     0.0
1700063	101325	8101325	101325	0	1	'gaslib-135'	1700063	0.0     0.0
1800076	101325	8101325	101325	0	1	'gaslib-135'	1800076	0.0     0.0
1900049	101325	8101325	101325	0	1	'gaslib-135'	1900049	0.0     0.0
2000056	101325	8101325	101325	0	1	'gaslib-135'	2000056	0.0     0.0
2100059	101325	8101325	101325	0	1	'gaslib-135'	2100059	0.0     0.0
2200056	101325	8101325	101325	0	1	'gaslib-135'	2200056	0.0     0.0
2300057	101325	8101325	101325	0	1	'gaslib-135'	2300057	0.0     0.0
2400066	101325	8101325	101325	0	1	'gaslib-135'	2400066	0.0     0.0
2500034	101325	8101325	101325	0	1	'gaslib-135'	2500034	0.0     0.0
2600005	101325	8101325	101325	0	1	'gaslib-135'	2600005	0.0     0.0
2700063	101325	8101325	101325	0	1	'gaslib-135'	2700063	0.0     0.0
2800034	101325	8101325	101325	0	1	'gaslib-135'	2800034	0.0     0.0
2900114	101325	8101325	101325	0	1	'gaslib-135'	2900114	0.0     0.0
];

%% pipe data
% id	fr_junction	to_junction	diameter	length	friction_factor	p_min	p_max	status
mgc.pipe = [
0	  2	  27	0.6	50304.7921	0.0078	101325	8101325	1
1	  130	24	1.0	100340.6882	0.0071	101325	8101325	1
2 	89	90	0.4	17722.0417	0.0085	101325	8101325	1
3	  89	90	0.4	17722.0417	0.0085	101325	8101325	1
4	  90	37	1.0	27575.9859	0.0071	101325	8101325	1
5	  89	6	  1.0	50039.0727	0.0071	101325	8101325	1
6	  4	  1	  1.0	35005.3067	0.0071	101325	8101325	1
7	  27	50	0.4	35640.0908	0.0085	101325	8101325	1
8	  50	37	0.6	12536.2361	0.0078	101325	8101325	1
9	  100	104	0.6	60939.7150	0.0078	101325	8101325	1
10	104	97	0.6	64597.2239	0.0078	101325	8101325	1
11	97	98	0.4	23354.6171	0.0085	101325	8101325	1
12	70	6	  1.0	61293.0070	0.0071	101325	8101325	1
13	98	29	0.6	43259.2640	0.0078	101325	8101325	1
14	29	85	0.6	39290.8095	0.0078	101325	8101325	1
15	85	86	0.6	38503.3883	0.0078	101325	8101325	1
16	86	79	1.0	65028.3898	0.0071	101325	8101325	1
17	95	16	1.0	162363.1287	0.0071	101325	8101325	1
18	81	80	0.6	6510.8232	  0.0078	101325	8101325	1
19	80	88	0.6	23219.4718	0.0078	101325	8101325	1
20	88	93	0.6	69525.6851	0.0078	101325	8101325	1
21	88	93	1.0	69525.6851	0.0071	101325	8101325	1
22	93	92	0.6	29067.2163	0.0078	101325	8101325	1
23	70	72	1.0	35220.8674	0.0071	101325	8101325	1
24	92	79	1.0	73550.1565	0.0071	101325	8101325	1
25	93	103	0.6	46081.5965	0.0078	101325	8101325	1
26	103	87	0.6	49928.7695	0.0078	101325	8101325	1
27	93	101	0.6	85857.9821	0.0078	101325	8101325	1
28	93	101	1.0	85857.9821	0.0071	101325	8101325	1
29	101	10	0.6	10634.4893	0.0078	101325	8101325	1
30	101	10	1.0	10634.4893	0.0071	101325	8101325	1
31	79	65	1.0	137584.5291	0.0071	101325	8101325	1
32	101	84	1.0	37563.1488	0.0071	101325	8101325	1
33	84	82	1.0	74204.5482	0.0071	101325	8101325	1
34	72	20	1.0	44307.5784	0.0071	101325	8101325	1
35	82	78	0.6	34864.4565	0.0078	101325	8101325	1
36	96	78	1.0	10819.0309	0.0071	101325	8101325	1
37	96	102	0.6	13882.5906	0.0078	101325	8101325	1
38	102	13	1.0	59093.6132	0.0071	101325	8101325	1
39	78	91	0.6	74213.0364	0.0078	101325	8101325	1
40	91	8	  1.0	10380.0380	0.0071	101325	8101325	1
41	8	  99	1.0	39573.0643	0.0071	101325	8101325	1
42	99	14	1.0	25443.4906	0.0071	101325	8101325	1
43	91	12	1.0	81583.8293	0.0071	101325	8101325	1
44	12	67	1.0	49535.2363	0.0071	101325	8101325	1
45	6	  64	1.0	60573.4639	0.0071	101325	8101325	1
46	11	26	1.0	19740.0494	0.0071	101325	8101325	1
47	9	  11	0.6	25479.1293	0.0078	101325	8101325	1
48	64	65	1.0	94839.1666	0.0071	101325	8101325	1
49	65	29	1.0	89395.5362	0.0071	101325	8101325	1
50	29	83	1.0	67770.9172	0.0071	101325	8101325	1
51	83	16	1.0	10478.3818	0.0071	101325	8101325	1
52	16	7	  0.6	9209.8789	  0.0078	101325	8101325	1
53	27	50	1.0	35640.0908	0.0071	101325	8101325	1
54	16	36	1.0	33603.2761	0.0071	101325	8101325	1
55	71	20	1.0	96065.0297	0.0071	101325	8101325	1
56	71	24	1.0	61286.1713	0.0071	101325	8101325	1
57	71	21	1.0	60740.4522	0.0071	101325	8101325	1
58	128	74	1.0	32476.9494	0.0071	101325	8101325	1
59	129	19	1.0	26313.6685	0.0071	101325	8101325	1
60	74	75	1.0	3274.8442	  0.0071	101325	8101325	1
61	75	48	1.0	54490.9543	0.0071	101325	8101325	1
62	48	49	1.0	17417.5481	0.0071	101325	8101325	1
63	49	66	1.0	77016.0840	0.0071	3101325	8101325	1
64	50	37	1.0	12536.2361	0.0071	101325	8101325	1
65	66	15	1.0	9158.2217	  0.0071	101325	8101325	1
66	66	15	1.0	9158.2217	  0.0071	101325	8101325	1
67	66	51	1.0	8869.5424	  0.0071	101325	8101325	1
68	66	51	1.0	8869.5424	  0.0071	101325	8101325	1
69	19	38	1.0	173659.7030	0.0071	101325	8101325	1
70	127	59	1.0	85824.9497	0.0071	101325	8101325	1
71	38	59	1.0	112865.0129	0.0071	101325	8101325	1
72	107	62	1.0	79066.7411	0.0071	101325	8101325	1
73	106	44	1.0	169915.7859	0.0071	101325	8101325	1
74	62	44	1.0	99696.9688	0.0071	101325	8101325	1
75	1	  25	1.0	42281.4742	0.0071	101325	8101325	1
76	44	55	0.6	92797.1469	0.0078	101325	8101325	1
77	119	41	1.0	24981.3111	0.0071	101325	8101325	1
78	120	41	1.0	25743.4219	0.0071	101325	8101325	1
79	44	123	1.0	84458.1686	0.0071	101325	8101325	1
80	44	124	1.0	84946.6519	0.0071	101325	8101325	1
81	125	30	1.0	80621.8096	0.0071	101325	8101325	1
82	126	30	1.0	80904.3281	0.0071	101325	8101325	1
83	121	45	1.0	83082.9611	0.0071	101325	8101325	1
84	122	45	1.0	83957.6804	0.0071	101325	8101325	1
85	45	63	1.0	13264.2145	0.0071	101325	8101325	1
86	25	77	1.0	49223.1218	0.0071	101325	8101325	1
87	45	63	1.0	13264.2145	0.0071	101325	8101325	1
88	109	34	1.0	87864.0574	0.0071	101325	8101325	1
89	108	34	1.0	87454.8381	0.0071	101325	8101325	1
90	131	69	1.0	33438.2062	0.0071	101325	8101325	1
91	132	69	1.0	33473.8697	0.0071	101325	8101325	1
92	38	52	1.0	42519.4409	0.0071	101325	8101325	1
93	52	67	1.0	21281.5244	0.0071	101325	8101325	1
94	38	58	1.0	37246.4955	0.0071	101325	8101325	1
95	58	57	1.0	63409.1675	0.0071	101325	8101325	1
96	30	31	1.0	13796.8759	0.0071	101325	8101325	1
97	4	  77	1.0	10397.5683	0.0071	101325	8101325	1
98	30	31	1.0	13796.8759	0.0071	101325	8101325	1
99	30	31	0.6	13796.8759	0.0078	101325	8101325	1
100	30	39	1.0	14097.1093	0.0071	101325	8101325	1
101	30	39	0.6	14097.1093	0.0078	101325	8101325	1
102	31	115	1.0	81461.5321	0.0071	101325	8101325	1
103	31	116	1.0	82408.7726	0.0071	101325	8101325	1
104	31	33	0.6	47648.2478	0.0078	101325	8101325	1
105	33	133	1.0	37524.1310	0.0071	101325	8101325	1
106	57	54	1.0	55637.0878	0.0071	101325	8101325	1
107	57	117	1.0	89446.2275	0.0071	101325	8101325	1
108	23	77	1.0	6844.7490	  0.0071	101325	8101325	1
109	57	118	1.0	89179.2180	0.0071	101325	8101325	1
110	54	42	0.6	48033.1621	0.0078	101325	8101325	1
111	54	46	1.0	19465.5730	0.0071	101325	8101325	1
112	46	32	1.0	93300.8611	0.0071	101325	8101325	1
113	46	56	1.0	31776.6029	0.0071	101325	8101325	1
114	56	110	1.0	75196.4294	0.0071	101325	8101325	1
115	56	134	1.0	61179.7197	0.0071	101325	8101325	1
116	56	134	1.0	61179.7197	0.0071	101325	8101325	1
117	134	111	1.0	27318.1496	0.0071	101325	7101325	1
118	134	112	1.0	26466.9067	0.0071	101325	7101325	1
119	17	23	1.0	22659.1404	0.0071	101325	8101325	1
120	113	61	1.0	37045.9904	0.0071	101325	8101325	1
121	114	61	1.0	37887.8337	0.0071	101325	8101325	1
122	61	53	1.0	51079.7763	0.0071	101325	8101325	1
123	61	105	1.0	135035.4911	0.0071	101325	8101325	1
124	76	3	  1.0	3000.0085	  0.0071	101325	8101325	1
125	32	73	1.0	28202.1519	0.0071	101325	8101325	1
126	73	22	0.4	56443.9175	0.0085	101325	8101325	1
127	22	40	0.6	78489.7647	0.0078	101325	8101325	1
128	73	35	0.4	50533.8662	0.0085	101325	8101325	1
129	68	35	0.6	12358.4449	0.0078	101325	8101325	1
130	25	70	1.0	52638.2656	0.0071	101325	8101325	1
131	73	43	1.0	22865.4790	0.0071	101325	8101325	1
132	43	26	1.0	70219.5178	0.0071	101325	8101325	1
133	53	28	1.0	6799.7886	  0.0071	101325	8101325	1
134	28	26	1.0	60895.6319	0.0071	101325	8101325	1
135	53	94	1.0	73579.2704	0.0071	101325	8101325	1
136	94	60	1.0	7053.1298	  0.0071	101325	8101325	1
137	26	60	1.0	44718.9773	0.0071	101325	8101325	1
138	60	18	1.0	38054.0887	0.0071	101325	8101325	1
139	60	0	  1.0	15990.7836	0.0071	101325	8101325	1
140	80	47	0.6	9159.0147	  0.0078	101325	8101325	1
];

%% compressor data
% id	fr_junction	to_junction	c_ratio_min	c_ratio_max	power_max	flow_min	flow_max	inlet_p_min	inlet_p_max	outlet_p_min	outlet_p_max	status	operating_cost	directionality
mgc.compressor = [
141	    17	500017	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
142	    76	1800076	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
143	    66	1100066	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
144	    66	2400066	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
145	    63	2700063	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
146	    63	1700063	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
147	    34	2800034	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
148	    34	2500034	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
149	    5	  900005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
150	    5	  1300005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
151	    5	  2600005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
152	    113	1600113	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
153	    114	2900114	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
154	    57	800057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
155	    57	100057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
156	    57	2300057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
157	    56	2000056	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
158	    56	2200056	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
159	    44	200044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
160	    44	700044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
161	    44	1500044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
162	    44	400044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
163	    30	600030	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
164	    30	1400030	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
165	    59	2100059	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
166	    59	1000059	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
167	    49	1900049	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
168	    71	300071	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
169	    71	1200071	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100000	116	100057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100001	119	200044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100002	128	300071	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100003	122	400044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100004	130	500017	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100005	123	600030	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100006	120	700044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100007	115	800057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100008	110	900005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100009	126	1000059	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100010	106	1100066	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100011	129	1200071	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100012	111	1300005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100013	124	1400030	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100014	121	1500044	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100015	134	1600113	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100016	109	1700063	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100017	105	1800076	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100018	127	1900049	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100019	117	2000056	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100020	125	2100059	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100021	118	2200056	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100022	133	2300057	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100023	107	2400066	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100024	132	2500034	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100025	112	2600005	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100026	108	2700063	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100027	131	2800034	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
100028	134	2900114	1.0	5.0	1000000000	-3300 3300 101325	8101325	101325	8101325	1	10	2
];

%% receipt data
% id	junction_id	injection_min	injection_max	injection_nominal	is_dispatchable	status
mgc.receipt = [
0	0	0.0	390.0	389.4	1	1
1	1	0.0	389.4	389.4	0	1
2	2	0.0	389.4	389.4	0	1
3	3	0.0	389.4	389.4	0	1
4	4	0.0	389.4	389.4	0	1
5	5	0.0	389.4	389.4	0	1
];

%% delivery data
% id	junction_id	withdrawal_min	withdrawal_max	withdrawal_nominal	is_dispatchable	status
mgc.delivery = [
6	  6	   0	23.6	23.6	0	1
7	  7	   0	23.6	23.6	0	1
8	  8	   0	23.6	23.6	0	1
9	  9	   0	23.6	23.6	0	1
10	10	 0	23.6	23.6	0	1
11	11	 0	23.6	23.6	0	1
12	12	 0	23.6	23.6	0	1
13	13	 0	23.6	23.6	0	1
14	14	 0	23.6	23.6	0	1
15	15	 0	23.6	23.6	0	1
16	16	 0	23.6	23.6	0	1
17	17	 0	23.6	23.6	0	1
18	18	 0	23.6	23.6	0	1
19	19	 0	23.6	23.6	0	1
20	20	 0	23.6	23.6	0	1
21	21	 0	23.6	23.6	0	1
22	22	 0	23.6	23.6	0	1
23	23	 0	23.6	23.6	0	1
24	24	 0	23.6	23.6	0	1
25	25	 0	23.6	23.6	0	1
26	26	 0	23.6	23.6	0	1
27	27	 0	23.6	23.6	0	1
28	28	 0	23.6	23.6	0	1
29	29	 0	23.6	23.6	0	1
30	30	 0	23.6	23.6	0	1
31	31	 0	23.6	23.6	0	1
32	32	 0	23.6	23.6	0	1
33	33	 0	23.6	23.6	0	1
34	34	 0	23.6	23.6	0	1
35	35	 0	23.6	23.6	0	1
36	36	 0	23.6	23.6	0	1
37	37	 0	23.6	23.6	0	1
38	38	 0	23.6	23.6	0	1
39	39	 0	23.6	23.6	0	1
40	40	 0	23.6	23.6	0	1
41	41	 0	23.6	23.6	0	1
42	42	 0	23.6	23.6	0	1
43	43	 0	23.6	23.6	0	1
44	44	 0	23.6	23.6	0	1
45	45	 0	23.6	23.6	0	1
46	46	 0	23.6	23.6	0	1
47	47	 0	23.6	23.6	0	1
48	48	 0	23.6	23.6	0	1
49	49	 0	23.6	23.6	0	1
50	50	 0	23.6	23.6	0	1
51	51	 0	23.6	23.6	0	1
52	52	 0	23.6	23.6	0	1
53	53	 0	23.6	23.6	0	1
54	54	 0	23.6	23.6	0	1
55	55	 0	23.6	23.6	0	1
56	56	 0	23.6	23.6	0	1
57	57	 0	23.6	23.6	0	1
58	58	 0	23.6	23.6	0	1
59	59	 0	23.6	23.6	0	1
60	60	 0	23.6	23.6	0	1
61	61	 0	23.6	23.6	0	1
62	62	 0	23.6	23.6	0	1
63	63	 0	23.6	23.6	0	1
64	64	 0	23.6	23.6	0	1
65	65	 0	23.6	23.6	0	1
66	66	 0	23.6	23.6	0	1
67	67	 0	23.6	23.6	0	1
68	68	 0	23.6	23.6	0	1
69	69	 0	23.6	23.6	0	1
70	70	 0	23.6	23.6	0	1
71	71	 0	23.6	23.6	0	1
72	72	 0	23.6	23.6	0	1
73	73	 0	23.6	23.6	0	1
74	74	 0	23.6	23.6	0	1
75	75	 0	23.6	23.6	0	1
76	76	 0	23.6	23.6	0	1
77	77	 0	23.6	23.6	0	1
78	78	 0	23.6	23.6	0	1
79	79	 0	23.6	23.6	0	1
80	80	 0	23.6	23.6	0	1
81	81	 0	23.6	23.6	0	1
82	82	 0	23.6	23.6	0	1
83	83	 0	23.6	23.6	0	1
84	84	 0	23.6	23.6	0	1
85	85	 0	23.6	23.6	0	1
86	86	 0	23.6	23.6	0	1
87	87	 0	23.6	23.6	0	1
88	88	 0	23.6	23.6	0	1
89	89	 0	23.6	23.6	0	1
90	90	 0	23.6	23.6	0	1
91	91	 0	23.6	23.6	0	1
92	92	 0	23.6	23.6	0	1
93	93	 0	23.6	23.6	0	1
94	94	 0	23.6	23.6	0	1
95	95	 0	23.6	23.6	0	1
96	96	 0	23.6	23.6	0	1
97	97	 0	23.6	23.6	0	1
98	98	 0	23.6	23.6	0	1
99	99	 0	23.6	23.6	0	1
100	100	 0	23.6	23.6	0	1
101	101	 0	23.6	23.6	0	1
102	102	 0	23.6	23.6	0	1
103	103	 0	23.6	23.6	0	1
104	104	 0	23.6	23.6	0	1
];

%% ne_pipe data
% id	fr_junction	to_junction	diameter	length	friction_factor	p_min	p_max	status	construction_cost
mgc.ne_pipe = [
171	2	  27	0.6	50304.7921	0.0078	101325	8101325	1	60.440
172	130	24	1.0	100340.6882	0.0071	101325	8101325	1	207.476
173	89	90	0.4	17722.0417	0.0085	101325	8101325	1	15.049
174	89	90	0.4	17722.0417	0.0085	101325	8101325	1	15.049
175	90	37	1.0	27575.9859	0.0071	101325	8101325	1	57.019
176	89	6	  1.0	50039.0727	0.0071	101325	8101325	1	103.466
177	4	  1	  1.0	35005.3067	0.0071	101325	8101325	1	72.381
178	27	50	0.4	35640.0908	0.0085	101325	8101325	1	30.264
179	50	37	0.6	12536.2361	0.0078	101325	8101325	1	15.062
180	100	104	0.6	60939.7150	0.0078	101325	8101325	1	73.217
181	104	97	0.6	64597.2239	0.0078	101325	8101325	1	77.612
182	97	98	0.4	23354.6171	0.0085	101325	8101325	1	19.832
183	70	6	  1.0	61293.0070	0.0071	101325	8101325	1	126.736
184	98	29	0.6	43259.2640	0.0078	101325	8101325	1	51.975
185	29	85	0.6	39290.8095	0.0078	101325	8101325	1	47.207
186	85	86	0.6	38503.3883	0.0078	101325	8101325	1	46.261
187	86	79	1.0	65028.3898	0.0071	101325	8101325	1	134.460
188	95	16	1.0	162363.1287	0.0071	101325	8101325	1	335.720
189	81	80	0.6	6510.8232	  0.0078	101325	8101325	1	7.823
190	80	88	0.6	23219.4718	0.0078	101325	8101325	1	27.897
191	88	93	0.6	69525.6851	0.0078	101325	8101325	1	83.533
192	88	93	1.0	69525.6851	0.0071	101325	8101325	1	143.759
193	93	92	0.6	29067.2163	0.0078	101325	8101325	1	34.923
194	70	72	1.0	35220.8674	0.0071	101325	8101325	1	72.827
195	92	79	1.0	73550.1565	0.0071	101325	8101325	1	152.080
196	93	103	0.6	46081.5965	0.0078	101325	8101325	1	55.366
197	103	87	0.6	49928.7695	0.0078	101325	8101325	1	59.988
198	93	101	0.6	85857.9821	0.0078	101325	8101325	1	103.156
199	93	101	1.0	85857.9821	0.0071	101325	8101325	1	177.530
200	101	10	0.6	10634.4893	0.0078	101325	8101325	1	12.777
201	101	10	1.0	10634.4893	0.0071	101325	8101325	1	21.989
202	79	65	1.0	137584.5291	0.0071	101325	8101325	1	284.485
203	101	84	1.0	37563.1488	0.0071	101325	8101325	1	77.679
204	84	82	1.0	74204.5482	0.0071	101325	8101325	1	153.434
205	72	20	1.0	44307.5784	0.0071	101325	8101325	1	91.615
206	82	78	0.6	34864.4565	0.0078	101325	8101325	1	41.889
207	96	78	1.0	10819.0309	0.0071	101325	8101325	1	22.371
208	96	102	0.6	13882.5906	0.0078	101325	8101325	1	16.680
209	102	13	1.0	59093.6132	0.0071	101325	8101325	1	122.189
210	78	91	0.6	74213.0364	0.0078	101325	8101325	1	89.165
211	91	8	  1.0	10380.0380	0.0071	101325	8101325	1	21.463
212	8	  99	1.0	39573.0643	0.0071	101325	8101325	1	81.826
213	99	14	1.0	25443.4906	0.0071	101325	8101325	1	52.610
214	91	12	1.0	81583.8293	0.0071	101325	8101325	1	168.692
215	12	67	1.0	49535.2363	0.0071	101325	8101325	1	102.425
216	6	  64	1.0	60573.4639	0.0071	101325	8101325	1	125.248
217	11	26	1.0	19740.0494	0.0071	101325	8101325	1	40.817
218	9	  11	0.6	25479.1293	0.0078	101325	8101325	1	30.612
219	64	65	1.0	94839.1666	0.0071	101325	8101325	1	196.100
220	65	29	1.0	89395.5362	0.0071	101325	8101325	1	184.844
221	29	83	1.0	67770.9172	0.0071	101325	8101325	1	140.131
222	83	16	1.0	10478.3818	0.0071	101325	8101325	1	21.666
223	16	7	  0.6	9209.8789	  0.0078	101325	8101325	1	11.065
224	27	50	1.0	35640.0908	0.0071	101325	8101325	1	73.693
225	16	36	1.0	33603.2761	0.0071	101325	8101325	1	69.482
226	71	20	1.0	96065.0297	0.0071	101325	8101325	1	198.635
227	71	24	1.0	61286.1713	0.0071	101325	8101325	1	126.722
228	71	21	1.0	60740.4522	0.0071	101325	8101325	1	125.594
229	128	74	1.0	32476.9494	0.0071	101325	8101325	1	67.153
230	129	19	1.0	26313.6685	0.0071	101325	8101325	1	54.409
231	74	75	1.0	3274.8442	  0.0071	101325	8101325	1	6.771
232	75	48	1.0	54490.9543	0.0071	101325	8101325	1	112.672
233	48	49	1.0	17417.5481	0.0071	101325	8101325	1	36.0145
234	49	66	1.0	77016.0840	0.0071	3101325	8101325	1	159.247
235	50	37	1.0	12536.2361	0.0071	101325	8101325	1	25.921
236	66	15	1.0	9158.2217	  0.0071	101325	8101325	1	18.937
237	66	15	1.0	9158.2217  	0.0071	101325	8101325	1	18.937
238	66	51	1.0	8869.5424	  0.0071	101325	8101325	1	18.340
239	66	51	1.0	8869.5424	  0.0071	101325	8101325	1	18.340
240	19	38	1.0	173659.7030	0.0071	101325	8101325	1	359.078
241	127	59	1.0	85824.9497	0.0071	101325	8101325	1	177.461
242	38	59	1.0	112865.0129	0.0071	101325	8101325	1	233.372
243	107	62	1.0	79066.7411	0.0071	101325	8101325	1	163.487
244	106	44	1.0	169915.7859	0.0071	101325	8101325	1	351.337
245	62	44	1.0	99696.9688	0.0071	101325	8101325	1	206.145
246	1	  25	1.0	42281.4742	0.0071	101325	8101325	1	87.426
247	44	55	0.6	92797.1469	0.0078	101325	8101325	1	111.493
248	119	41	1.0	24981.3111	0.0071	101325	8101325	1	051.654
249	120	41	1.0	25743.4219	0.0071	101325	8101325	1	53.230
250	44	123	1.0	84458.1686	0.0071	101325	8101325	1	174.635
251	44	124	1.0	84946.6519	0.0071	101325	8101325	1	175.645
252	125	30	1.0	80621.8096	0.0071	101325	8101325	1	166.703
253	126	30	1.0	80904.3281	0.0071	101325	8101325	1	167.287
254	121	45	1.0	83082.9611	0.0071	101325	8101325	1	171.792
255	122	45	1.0	83957.6804	0.0071	101325	8101325	1	173.600
256	45	63	1.0	13264.2145	0.0071	101325	8101325	1	27.427
257	25	77	1.0	49223.1218	0.0071	101325	8101325	1	101.779
258	45	63	1.0	13264.2145	0.0071	101325	8101325	1	27.427
259	109	34	1.0	87864.0574	0.0071	101325	8101325	1	181.677
260	108	34	1.0	87454.8381	0.0071	101325	8101325	1	180.831
261	131	69	1.0	33438.2062	0.0071	101325	8101325	1	69.141
262	132	69	1.0	33473.8697	0.0071	101325	8101325	1	69.214
263	38	52	1.0	42519.4409	0.0071	101325	8101325	1	87.918
264	52	67	1.0	21281.5244	0.0071	101325	8101325	1	44.004
265	38	58	1.0	37246.4955	0.0071	101325	8101325	1	77.015
266	58	57	1.0	63409.1675	0.0071	101325	8101325	1	131.112
267	30	31	1.0	13796.8759	0.0071	101325	8101325	1	28.528
268	4	  77	1.0	10397.5683	0.0071	101325	8101325	1	21.499
269	30	31	1.0	13796.8759	0.0071	101325	8101325	1	28.528
270	30	31	0.6	13796.8759	0.0078	101325	8101325	1	16.577
271	30	39	1.0	14097.1093	0.0071	101325	8101325	1	29.149
272	30	39	0.6	14097.1093	0.0078	101325	8101325	1	16.937
273	31	115	1.0	81461.5321	0.0071	101325	8101325	1	168.439
274	31	116	1.0	82408.7726	0.0071	101325	8101325	1	170.398
275	31	33	0.6	47648.2478	0.0078	101325	8101325	1	57.248
276	33	133	1.0	37524.1310	0.0071	101325	8101325	1	77.589
277	57	54	1.0	55637.0878	0.0071	101325	8101325	1	115.041
278	57	117	1.0	89446.2275	0.0071	101325	8101325	1	184.949
279	23	77	1.0	6844.7490	  0.0071	101325	8101325	1	14.153
280	57	118	1.0	89179.2180	0.0071	101325	8101325	1	184.397
281	54	42	0.6	48033.1621	0.0078	101325	8101325	1	57.710
282	54	46	1.0	19465.5730	0.0071	101325	8101325	1	40.249
283	46	32	1.0	93300.8611	0.0071	101325	8101325	1	192.919
284	46	56	1.0	31776.6029	0.0071	101325	8101325	1	65.705
285	56	110	1.0	75196.4294	0.0071	101325	8101325	1	155.485
286	56	134	1.0	61179.7197	0.0071	101325	8101325	1	126.502
287	56	134	1.0	61179.7197	0.0071	101325	8101325	1	126.502
288	134	111	1.0	27318.1496	0.0071	101325	7101325	1	56.486
289	134	112	1.0	26466.9067	0.0071	101325	7101325	1	54.726
290	17	23	1.0	22659.1404	0.0071	101325	8101325	1	46.853
291	113	61	1.0	37045.9904	0.0071	101325	8101325	1	76.600
292	114	61	1.0	37887.8337	0.0071	101325	8101325	1	78.341
293	61	53	1.0	51079.7763	0.0071	101325	8101325	1	105.618
294	61	105	1.0	135035.4911	0.0071	101325	8101325	1	279.214
295	76	3	  1.0	3000.0085	  0.0071	101325	8101325	1	6.203
296	32	73	1.0	28202.1519	0.0071	101325	8101325	1	58.314
297	73	22	0.4	56443.9175	0.0085	101325	8101325	1	47.929
298	22	40	0.6	78489.7647	0.0078	101325	8101325	1	94.303
299	73	35	0.4	50533.8662	0.0085	101325	8101325	1	42.911
300	68	35	0.6	12358.4449	0.0078	101325	8101325	1	14.848
301	25	70	1.0	52638.2656	0.0071	101325	8101325	1	108.841
302	73	43	1.0	22865.4790	0.0071	101325	8101325	1	47.279
303	43	26	1.0	70219.5178	0.0071	101325	8101325	1	145.194
304	53	28	1.0	6799.7886	  0.0071	101325	8101325	1	14.060
305	28	26	1.0	60895.6319	0.0071	101325	8101325	1	125.915
306	53	94	1.0	73579.2704	0.0071	101325	8101325	1	152.141
307	94	60	1.0	7053.1298	  0.0071	101325	8101325	1	14.584
308	26	60	1.0	44718.9773	0.0071	101325	8101325	1	92.466
309	60	18	1.0	38054.0887	0.0071	101325	8101325	1	78.685
310	60	0	  1.0	15990.7836	0.0071	101325	8101325	1	33.064
311	80	47	0.6	9159.0147	  0.0078	101325	8101325	1	11.004
];


end
