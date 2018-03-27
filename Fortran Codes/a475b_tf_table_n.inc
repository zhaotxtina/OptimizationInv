	include 'a475b_tf_table_def_n.inc'
c 
	data (p1_acal_475b(i),i=1,ntr)  / 100.4444401 ,100.3040603,
	1 100.3040532 ,100.2648325 ,100.2602424/
c 
	data (p2_acal_475b(i),i=1,ntr)  / 100.3497593 ,100.3552781,
	1 100.2738923 ,100.2856821 ,100.2466489/
c 
	data (ps_acal_475b(i),i=1,ntr)  /  -0.0946808 , -0.0512178,
	1  -0.0301608 , -0.0208496 , -0.0135935/
c 
	data (ad_acal_475b(i),i=1,ntr)  /  13.8232019 ,  8.9242878,
	1   6.6372094 ,  5.2981108 ,  4.4053324/
c 
	data (ps_offset_475b(i),i=1,ntr)  /   0.0932081 ,  0.0489266,
	1   0.0270390 ,  0.0168855 ,  0.0087769/
c 
	data (ru475b(i),i=1,ntf4)      /   0.0500000 ,  0.0700000,
	1   0.1000000 ,  0.1500000 ,  0.2000000 ,  0.3000000 ,  0.5000000,
	2   0.7000000 ,  1.0000000 ,  1.5000000 ,  2.0000000 ,  3.0000000,
	3   5.0000000 ,  7.0000000 , 10.0000000 , 15.0000000 , 20.0000000,
	4  30.0000000 , 50.0000000 , 70.0000000 ,100.0000000 ,150.0000000,
	5 200.0000000 ,300.0000000 ,500.0000000 ,700.0000000 ,1000.000000,
	6 1500.0000000 ,2000.0000000 ,3000.0000000 ,5000.0000000/
c 
c j =  1, TR spacing =     10 
c 
	data (psu475b(i,1),i=1,ntf4)  /  43.0124819 , 35.1401360,
	1  28.1423343 , 21.6245216 , 17.7986940 , 13.3691411 ,  9.1334173,
	2   7.0165597 ,  5.2476913 ,  3.7209334 ,  2.8901161 ,  1.9989129,
	3   1.2281140 ,  0.8767070 ,  0.6023230 ,  0.3809191 ,  0.2669905,
	4   0.1504171 ,  0.0548087 ,  0.0130315 , -0.0186945 , -0.0436518,
	5  -0.0562415 , -0.0689209 , -0.0791434 , -0.0835513 , -0.0868704,
	6  -0.0894615 , -0.0907608 , -0.0920633 , -0.0931081/
c 
	data (adu475b(i,1),i=1,ntf4)  /  17.8175242 , 16.8774913,
	1  16.0870161 , 15.4013435 , 15.0280641 , 14.6302412 , 14.2937273,
	2  14.1460980 , 14.0359419 , 13.9526568 , 13.9128989 , 13.8755717,
	3  13.8487793 , 13.8387016 , 13.8320085 , 13.8275361 , 13.8256350,
	4  13.8240445 , 13.8230874 , 13.8228029 , 13.8226623 , 13.8226149,
	5  13.8226208 , 13.8226565 , 13.8227195 , 13.8227627 , 13.8228061,
	6  13.8228507 , 13.8228790 , 13.8229141 , 13.8229510/
c 
	data (p1u475b(i,1),i=1,ntf4)  / 137.1365801 ,129.2335312,
	1 122.5106556 ,116.5728939 ,113.2654443 ,109.6332912 ,106.3943505,
	2 104.8764190 ,103.6673341 ,102.6725386 ,102.1524177 ,101.6131213,
	3 101.1642414 ,100.9657914 ,100.8139369 ,100.6936509 ,100.6326447,
	4 100.5709404 ,100.5209638 ,100.4993377 ,100.4830171 ,100.4702513,
	5 100.4638400 ,100.4574060 ,100.4522388 ,100.4500175 ,100.4483482,
	6 100.4470475 ,100.4463962 ,100.4457442 ,100.4452219/
c 
	data (p2u475b(i,1),i=1,ntf4)  / 180.1490620 ,164.3736672,
	1 150.6529899 ,138.1974155 ,131.0641383 ,123.0024323 ,115.5277678,
	2 111.8929787 ,108.9150254 ,106.3934720 ,105.0425338 ,103.6120342,
	3 102.3923554 ,101.8424984 ,101.4162599 ,101.0745700 ,100.8996352,
	4 100.7213575 ,100.5757725 ,100.5123692 ,100.4643226 ,100.4265995,
	5 100.4075985 ,100.3884851 ,100.3730954 ,100.3664662 ,100.3614778,
	6 100.3575860 ,100.3556354 ,100.3536809 ,100.3521138/
c 
c j =  2, TR spacing =     16 
c 
	data (psu475b(i,2),i=1,ntf4)  /  45.6568925 , 37.8443400,
	1  30.8532950 , 24.2661314 , 20.3415300 , 15.7100822 , 11.1414212,
	2   8.7787302 ,  6.7450587 ,  4.9307072 ,  3.9132074 ,  2.7910924,
	3   1.7876503 ,  1.3171207 ,  0.9425634 ,  0.6347893 ,  0.4741116,
	4   0.3077631 ,  0.1695721 ,  0.1085803 ,  0.0619628 ,  0.0250758,
	5   0.0063828 , -0.0125121 , -0.0278064 , -0.0344218 , -0.0394133,
	6  -0.0433174 , -0.0452780 , -0.0472458 , -0.0488266/
c 
	data (adu475b(i,2),i=1,ntf4)  /  13.8044079 , 12.7679781,
	1  11.8729010 , 11.0690758 , 10.6149859 , 10.1111346 ,  9.6594424,
	2   9.4492386 ,  9.2847260 ,  9.1537891 ,  9.0883202 ,  9.0242269,
	3   8.9758056 ,  8.9567377 ,  8.9436418 ,  8.9345758 ,  8.9305885,
	4   8.9271319 ,  8.9249235 ,  8.9242031 ,  8.9237944 ,  8.9235876,
	5   8.9235365 ,  8.9235372 ,  8.9235959 ,  8.9236476 ,  8.9237042,
	6   8.9237656 ,  8.9238056 ,  8.9238561 ,  8.9239100/
c 
	data (p1u475b(i,2),i=1,ntf4)  / 225.8475465 ,202.2566768,
	1 181.5413463 ,162.4941370 ,151.4330518 ,138.7355943 ,126.6875330,
	2 120.6874531 ,115.6735290 ,111.3355348 ,108.9659122 ,106.4119827,
	3 104.1876807 ,103.1667399 ,102.3655032 ,101.7156776 ,101.3798768,
	4 101.0350561 ,100.7511189 ,100.6266530 ,100.5319348 ,100.4572819,
	5 100.4195663 ,100.3815361 ,100.3508345 ,100.3375823 ,100.3275966,
	6 100.3197963 ,100.3158828 ,100.3119582 ,100.3088086/
c 
	data (p2u475b(i,2),i=1,ntf4)  / 180.1906540 ,164.4123368,
	1 150.6880513 ,138.2280056 ,131.0915218 ,123.0255121 ,115.5461118,
	2 111.9087229 ,108.9284703 ,106.4048276 ,105.0527048 ,103.6208903,
	3 102.4000304 ,101.8496192 ,101.4229398 ,101.0808883 ,100.9057652,
	4 100.7272930 ,100.5815468 ,100.5180727 ,100.4699720 ,100.4322061,
	5 100.4131835 ,100.3940482 ,100.3786409 ,100.3720041 ,100.3670099,
	6 100.3631137 ,100.3611608 ,100.3592040 ,100.3576352/
c 
c j =  3, TR spacing =     22 
c 
	data (psu475b(i,3),i=1,ntf4)  /  46.8592767 , 39.1073124,
	1  32.1654786 , 25.6085715 , 21.6849311 , 17.0220365 , 12.3582772,
	2   9.9030814 ,  7.7522651 ,  5.7912016 ,  4.6675074 ,  3.4014567,
	3   2.2376352 ,  1.6782899 ,  1.2250672 ,  0.8461828 ,  0.6455845,
	4   0.4354765 ,  0.2586677 ,  0.1798338 ,  0.1191796 ,  0.0708958,
	5   0.0463122 ,  0.0213694 ,  0.0010971 , -0.0076999 , -0.0143513,
	6  -0.0195635 , -0.0221851 , -0.0248197 , -0.0269390/
c 
	data (adu475b(i,3),i=1,ntf4)  /  12.0489856 , 10.9702010,
	1  10.0273642 ,  9.1664805 ,  8.6709169 ,  8.1087849 ,  7.5872582,
	2   7.3352843 ,  7.1315295 ,  6.9631863 ,  6.8759842 ,  6.7876987,
	3   6.7180999 ,  6.6896065 ,  6.6694715 ,  6.6551180 ,  6.6486346,
	4   6.6428692 ,  6.6390456 ,  6.6377380 ,  6.6369534 ,  6.6365096,
	5   6.6363660 ,  6.6362991 ,  6.6363309 ,  6.6363828 ,  6.6364470,
	6   6.6365210 ,  6.6365709 ,  6.6366354 ,  6.6367053/
c 
	data (p1u475b(i,3),i=1,ntf4)  / 225.8359913 ,202.2475519,
	1 181.5342435 ,162.4887671 ,151.4286279 ,138.7322215 ,126.6851995,
	2 120.6856876 ,115.6722639 ,111.3347021 ,108.9653034 ,106.4115946,
	3 104.1874593 ,103.1665850 ,102.3653956 ,101.7156052 ,101.3798213,
	4 101.0350172 ,100.7510928 ,100.6266324 ,100.5319183 ,100.4572685,
	5 100.4195545 ,100.3815258 ,100.3508255 ,100.3375738 ,100.3275885,
	6 100.3197886 ,100.3158752 ,100.3119508 ,100.3088013/
c 
	data (p2u475b(i,3),i=1,ntf4)  / 272.6952680 ,241.3548643,
	1 213.6997221 ,188.0973386 ,173.1135590 ,155.7542580 ,139.0434767,
	2 130.5887690 ,123.4245290 ,117.1259037 ,113.6328108 ,109.8130513,
	3 106.4250945 ,104.8448749 ,103.5904628 ,102.5617880 ,102.0254058,
	4 101.4704937 ,101.0097605 ,100.8064662 ,100.6510979 ,100.5281643,
	5 100.4658667 ,100.4028952 ,100.3519226 ,100.3298739 ,100.3132372,
	6 100.3002251 ,100.2936901 ,100.2871311 ,100.2818623/
c 
c j =  4, TR spacing =     28 
c 
	data (psu475b(i,4),i=1,ntf4)  /  47.4974724 , 39.7911085,
	1  32.8929329 , 26.3759267 , 22.4718031 , 17.8204290 , 13.1395774,
	2  10.6525928 ,  8.4514168 ,  6.4161933 ,  5.2324532 ,  3.8772177,
	3   2.6036379 ,  1.9786428 ,  1.4641741 ,  1.0272097 ,  0.7927686,
	4   0.5444290 ,  0.3327725 ,  0.2374359 ,  0.1635932 ,  0.1044498,
	5   0.0741922 ,  0.0433743 ,  0.0182220 ,  0.0072715 , -0.0010259,
	6  -0.0075407 , -0.0108228 , -0.0141251 , -0.0167855/
c 
	data (adu475b(i,4),i=1,ntf4)  /  11.0624415 ,  9.9616558,
	1   8.9933440 ,  8.1010635 ,  7.5819183 ,  6.9853212 ,  6.4199029,
	2   6.1399383 ,  5.9083652 ,  5.7117646 ,  5.6071532 ,  5.4983753,
	3   5.4095482 ,  5.3719588 ,  5.3447276 ,  5.3248106 ,  5.3156051,
	4   5.3072433 ,  5.3015351 ,  5.2995180 ,  5.2982664 ,  5.2975175,
	5   5.2972492 ,  5.2970848 ,  5.2970682 ,  5.2971120 ,  5.2971782,
	6   5.2972611 ,  5.2973193 ,  5.2973960 ,  5.2974808/
c 
	data (p1u475b(i,4),i=1,ntf4)  / 320.2167384 ,281.1705037,
	1 246.6168853 ,214.4963323 ,195.6073208 ,173.5950361 ,152.2016652,
	2 141.2590379 ,131.8927488 ,123.5579932 ,118.8805616 ,113.7047988,
	3 109.0424517 ,106.8368076 ,105.0675585 ,103.6015965 ,102.8305965,
	4 102.0271550 ,101.3546017 ,101.0558963 ,100.8266274 ,100.6445038,
	5 100.5519247 ,100.4581108 ,100.3819659 ,100.3489580 ,100.3240173,
	6 100.3044852 ,100.2946654 ,100.2848014 ,100.2768701/
c 
	data (p2u475b(i,4),i=1,ntf4)  / 272.7192660 ,241.3793952,
	1 213.7239524 ,188.1204056 ,173.1355177 ,155.7746071 ,139.0620878,
	2 130.6064451 ,123.4413320 ,117.1417999 ,113.6481084 ,109.8275811,
	3 106.4388138 ,104.8581648 ,103.6033844 ,102.5743868 ,102.0378279,
	4 101.4827260 ,101.0218292 ,100.8184604 ,100.6630342 ,100.5400540,
	5 100.4777325 ,100.4147365 ,100.3637439 ,100.3416865 ,100.3250432,
	6 100.3120259 ,100.3054882 ,100.2989265 ,100.2936556/
c 
c j =  5, TR spacing =     34 
c 
	data (psu475b(i,5),i=1,ntf4)  /  47.9040449 , 40.2233970,
	1  33.3527226 , 26.8652505 , 22.9792615 , 18.3466318 , 13.6729954,
	2  11.1782557 ,  8.9570536 ,  6.8849873 ,  5.6674220 ,  4.2569616,
	3   2.9083912 ,  2.2350643 ,  1.6731540 ,  1.1890200 ,  0.9260554,
	4   0.6444991 ,  0.4015539 ,  0.2910212 ,  0.2048362 ,  0.1353809,
	5   0.0996755 ,  0.0631666 ,  0.0332425 ,  0.0201711 ,  0.0102452,
	6   0.0024360 , -0.0015044 , -0.0054743 , -0.0086769/
c 
	data (adu475b(i,5),i=1,ntf4)  /  10.4208064 ,  9.3065678,
	1   8.3228647 ,  7.4115165 ,  6.8778274 ,  6.2594680 ,  5.6651979,
	2   5.3659964 ,  5.1145047 ,  4.8966700 ,  4.7783593 ,  4.6526987,
	3   4.5470651 ,  4.5010881 ,  4.4670496 ,  4.4415807 ,  4.4295664,
	4   4.4184474 ,  4.4106688 ,  4.4078484 ,  4.4060554 ,  4.4049430,
	5   4.4045214 ,  4.4042326 ,  4.4041476 ,  4.4041757 ,  4.4042386,
	6   4.4043266 ,  4.4043914 ,  4.4044788 ,  4.4045774/
c 
	data (p1u475b(i,5),i=1,ntf4)  / 320.2047556 ,281.1604555,
	1 246.6082266 ,214.4886954 ,195.6003164 ,173.5890800 ,152.1972394,
	2 141.2554944 ,131.8898450 ,123.5554005 ,118.8779670 ,113.7020030,
	3 109.0392611 ,106.8333515 ,105.0638503 ,103.5976518 ,102.8265168,
	4 102.0229264 ,101.3502414 ,101.0514753 ,100.8221587 ,100.6399965,
	5 100.5473975 ,100.4535633 ,100.3774016 ,100.3443865 ,100.3194402,
	6 100.2999037 ,100.2900818 ,100.2802156 ,100.2722825/
c 
	data (p2u475b(i,5),i=1,ntf4)  / 368.1088005 ,321.3838525,
	1 279.9609492 ,241.3539459 ,218.5795779 ,191.9357118 ,165.8702348,
	2 152.4337501 ,140.8468986 ,130.4403878 ,124.5453890 ,117.9589646,
	3 111.9476523 ,109.0684158 ,106.7370043 ,104.7866718 ,103.7525722,
	4 102.6674255 ,101.7517953 ,101.3424965 ,101.0269949 ,100.7753774,
	5 100.6470730 ,100.5167299 ,100.4106441 ,100.3645576 ,100.3296854,
	6 100.3023397 ,100.2885774 ,100.2747413 ,100.2636056/
c 