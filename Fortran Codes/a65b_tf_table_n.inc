	include 'a65b_tf_table_def_n.inc'
c 
	data (p1_acal_65b(i),i=1,ntr)  / 100.7477495 ,100.6820463,
	1 100.6822570 ,100.6570689 ,100.6543658/
c 
	data (p2_acal_65b(i),i=1,ntr)  / 100.7056114 ,100.7087859,
	1 100.6629471 ,100.6702640 ,100.6434197/
c 
	data (ps_acal_65b(i),i=1,ntr)  /  -0.0421381 , -0.0267396,
	1  -0.0193099 , -0.0131951 , -0.0109460/
c 
	data (ad_acal_65b(i),i=1,ntr)  /   8.6511970 ,  6.5001650,
	1   5.1791923 ,  4.3537245 ,  3.7121489/
c 
	data (ps_offset_65b(i),i=1,ntr)  /   0.0397272 ,  0.0234981,
	1   0.0152245 ,  0.0082557 ,  0.0051415/
c 
	data (ru65b(i),i=1,ntf4)      /   0.0500000 ,  0.0700000,
	1   0.1000000 ,  0.1500000 ,  0.2000000 ,  0.3000000 ,  0.5000000,
	2   0.7000000 ,  1.0000000 ,  1.5000000 ,  2.0000000 ,  3.0000000,
	3   5.0000000 ,  7.0000000 , 10.0000000 , 15.0000000 , 20.0000000,
	4  30.0000000 , 50.0000000 , 70.0000000 ,100.0000000 ,150.0000000,
	5 200.0000000 ,300.0000000 ,500.0000000 ,700.0000000 ,1000.000000,
	6 1500.0000000 ,2000.0000000 ,3000.0000000 ,5000.0000000/
c 
c j =  1, TR spacing =     16 
c 
	data (psu65b(i,1),i=1,ntf4)  /  46.0869225 , 38.2505960,
	1  31.2332682 , 24.6149753 , 20.6674630 , 16.0025243 , 11.3899690,
	2   8.9979176 ,  6.9335632 ,  5.0860198 ,  4.0467243 ,  2.8971154,
	3   1.8651202 ,  1.3795403 ,  0.9920510 ,  0.6728901 ,  0.5059433,
	4   0.3328262 ,  0.1887550 ,  0.1250782 ,  0.0763639 ,  0.0377853,
	5   0.0182225 , -0.0015622 , -0.0175856 , -0.0245196 , -0.0297531,
	6  -0.0338473 , -0.0359039 , -0.0379685 , -0.0396272/
c 
	data (adu65b(i,1),i=1,ntf4)  /  13.6402134 , 12.5902823,
	1  11.6817567 , 10.8637728 , 10.4003954 ,  9.8845847 ,  9.4198648,
	2   9.2024162 ,  9.0314216 ,  8.8945795 ,  8.8258008 ,  8.7581310,
	3   8.7066811 ,  8.6863010 ,  8.6722428 ,  8.6624665 ,  8.6581487,
	4   8.6543904 ,  8.6519742 ,  8.6511796 ,  8.6507241 ,  8.6504885,
	5   8.6504266 ,  8.6504197 ,  8.6504762 ,  8.6505286 ,  8.6505867,
	6   8.6506502 ,  8.6506918 ,  8.6507444 ,  8.6508008/
c 
	data (p1u65b(i,1),i=1,ntf4)  / 182.3516363 ,166.4055586,
	1 152.5034980 ,139.8424443 ,132.5662521 ,124.3118423 ,116.6172660,
	2 112.8556642 ,109.7608349 ,107.1286898 ,105.7130283 ,104.2088228,
	3 102.9210844 ,102.3385780 ,101.8859919 ,101.5223981 ,101.3359296,
	4 101.1456327 ,100.9899952 ,100.9221328 ,100.8706672 ,100.8302312,
	5 100.8098524 ,100.7893440 ,100.7728230 ,100.7657038 ,100.7603454,
	6 100.7561639 ,100.7540677 ,100.7519670 ,100.7502824/
c 
	data (p2u65b(i,1),i=1,ntf4)  / 228.4385588 ,204.6561546,
	1 183.7367662 ,164.4574196 ,153.2337151 ,140.3143666 ,128.0072350,
	2 121.8535818 ,116.6943981 ,112.2147096 ,109.7597526 ,107.1059382,
	3 104.7862046 ,103.7181183 ,102.8780429 ,102.1952882 ,101.8418729,
	4 101.4784589 ,101.1787502 ,101.0472110 ,100.9470311 ,100.8680165,
	5 100.8280749 ,100.7877818 ,100.7552374 ,100.7411842 ,100.7305923,
	6 100.7223166 ,100.7181638 ,100.7139985 ,100.7106552/
c 
c j =  2, TR spacing =     22 
c 
	data (psu65b(i,2),i=1,ntf4)  /  47.1726651 , 39.3975897,
	1  32.4317300 , 25.8487003 , 21.9075145 , 17.2212032 , 12.5294739,
	2  10.0563799 ,  7.8869059 ,  5.9051739 ,  4.7673924 ,  3.4827461,
	3   2.2983587 ,  1.7275391 ,  1.2640408 ,  0.8757414 ,  0.6697920,
	4   0.4537531 ,  0.2716434 ,  0.1903352 ,  0.1277215 ,  0.0778370,
	5   0.0524220 ,  0.0266224 ,  0.0056420 , -0.0034662 , -0.0103549,
	6  -0.0157545 , -0.0184710 , -0.0212012 , -0.0233981/
c 
	data (adu65b(i,2),i=1,ntf4)  /  11.9895851 , 10.9013589,
	1   9.9493177 ,  9.0789856 ,  8.5773243 ,  8.0073607 ,  7.4771370,
	2   7.2201364 ,  7.0116830 ,  6.8388136 ,  6.7489301 ,  6.6575851,
	3   6.5852106 ,  6.5554389 ,  6.5343246 ,  6.5192168 ,  6.5123697,
	4   6.5062620 ,  6.5021940 ,  6.5007957 ,  6.4999523 ,  6.4994708,
	5   6.4993122 ,  6.4992340 ,  6.4992607 ,  6.4993122 ,  6.4993771,
	6   6.4994529 ,  6.4995042 ,  6.4995704 ,  6.4996427/
c 
	data (p1u65b(i,2),i=1,ntf4)  / 275.3088072 ,243.7845463,
	1 215.9341228 ,190.1098706 ,174.9704728 ,157.3978256 ,140.4349221,
	2 131.8282082 ,124.5176230 ,118.0729916 ,114.4899534 ,110.5624454,
	3 107.0683760 ,105.4342834 ,104.1345965 ,103.0667798 ,102.5091238,
	4 101.9314557 ,101.4511348 ,101.2389535 ,101.0766715 ,100.9481791,
	5 100.8830293 ,100.8171463 ,100.7637918 ,100.7407042 ,100.7232794,
	6 100.7096478 ,100.7028004 ,100.6959270 ,100.6904046/
c 
	data (p2u65b(i,2),i=1,ntf4)  / 228.1361421 ,204.3869566,
	1 183.5023928 ,164.2611703 ,153.0629583 ,140.1766224 ,127.9054482,
	2 121.7718283 ,116.6307171 ,112.1678177 ,109.7225610 ,107.0796993,
	3 104.7700173 ,103.7067443 ,102.8705557 ,102.1910384 ,101.8393318,
	4 101.4777026 ,101.1794914 ,101.0486183 ,100.9489500 ,100.8703421,
	5 100.8306073 ,100.7905239 ,100.7581498 ,100.7441704 ,100.7336343,
	6 100.7254023 ,100.7212714 ,100.7171282 ,100.7138027/
c 
c j =  3, TR spacing =     28 
c 
	data (psu65b(i,3),i=1,ntf4)  /  47.7498395 , 40.0230911,
	1  33.1032746 , 26.5627673 , 22.6431563 , 17.9719067 , 13.2690123,
	2  10.7688890 ,  8.5545647 ,  6.5050379 ,  5.3114920 ,  3.9430121,
	3   2.6541312 ,  2.0202001 ,  1.4974336 ,  1.0525831 ,  0.8135228,
	4   0.5599318 ,  0.3434480 ,  0.2458081 ,  0.1701152 ,  0.1094407,
	5   0.0783802 ,  0.0467283 ,  0.0208807 ,  0.0096228 ,  0.0010900,
	6  -0.0056115 , -0.0089883 , -0.0123864 , -0.0151245/
c 
	data (adu65b(i,3),i=1,ntf4)  /  11.0051060 ,  9.8964995,
	1   8.9207723 ,  8.0211056 ,  7.4973216 ,  6.8949017 ,  6.3231175,
	2   6.0394541 ,  5.8043544 ,  5.6042390 ,  5.4974631 ,  5.3861096,
	3   5.2948099 ,  5.2560206 ,  5.2278340 ,  5.2071520 ,  5.1975653,
	4   5.1888345 ,  5.1828539 ,  5.1807327 ,  5.1794120 ,  5.1786174,
	5   5.1783301 ,  5.1781509 ,  5.1781267 ,  5.1781692 ,  5.1782354,
	6   5.1783195 ,  5.1783789 ,  5.1784574 ,  5.1785442/
c 
	data (p1u65b(i,3),i=1,ntf4)  / 275.1116001 ,243.6065048,
	1 215.7782060 ,189.9792137 ,174.8568447 ,157.3062092 ,140.3672383,
	2 131.7738228 ,124.4751501 ,118.0414461 ,114.4646375 ,110.5440835,
	3 107.0563347 ,105.4252593 ,104.1280121 ,103.0622352 ,102.5056610,
	4 101.9291286 ,101.4497666 ,101.2380143 ,101.0760630 ,100.9478346,
	5 100.8828194 ,100.8170732 ,100.7638301 ,100.7407909 ,100.7234027,
	6 100.7097999 ,100.7029670 ,100.6961080 ,100.6905974/
c 
	data (p2u65b(i,3),i=1,ntf4)  / 322.8614396 ,283.6295959,
	1 248.8814806 ,216.5419810 ,197.5000010 ,175.2781159 ,153.6362506,
	2 142.5427118 ,133.0297148 ,124.5464840 ,119.7761295 ,114.4870956,
	3 109.7104659 ,107.4454594 ,105.6254457 ,104.1148183 ,103.3191838,
	4 102.4890604 ,101.7932146 ,101.4838224 ,101.2461782 ,101.0572753,
	5 100.9611996 ,100.8638015 ,100.7847108 ,100.7504137 ,100.7244927,
	6 100.7041884 ,100.6939787 ,100.6837216 ,100.6754729/
c 
c j =  4, TR spacing =     34 
c 
	data (psu65b(i,4),i=1,ntf4)  /  48.1324166 , 40.4311792,
	1  33.5385729 , 27.0272814 , 23.1257702 , 18.4737073 , 13.7796098,
	2  11.2734149 ,  9.0413251 ,  6.9579262 ,  5.7327653 ,  4.3120798,
	3   2.9515007 ,  2.2709660 ,  1.7021858 ,  1.2113324 ,  0.9443319,
	4   0.6580847 ,  0.4107168 ,  0.2980310 ,  0.2100931 ,  0.1391697,
	5   0.1026868 ,  0.0653645 ,  0.0347570 ,  0.0213812 ,  0.0112214,
	6   0.0032261 , -0.0008089 , -0.0048749 , -0.0081557/
c 
	data (adu65b(i,4),i=1,ntf4)  /  10.4217615 ,  9.3005683,
	1   8.3103013 ,  7.3925454 ,  6.8549342 ,  6.2317787 ,  5.6324053,
	2   5.3302741 ,  5.0759813 ,  4.8553132 ,  4.7352184 ,  4.6073748,
	3   4.4995567 ,  4.4524753 ,  4.4175286 ,  4.3913078 ,  4.3789080,
	4   4.3674061 ,  4.3593358 ,  4.3564009 ,  4.3545302 ,  4.3533653,
	5   4.3529215 ,  4.3526145 ,  4.3525193 ,  4.3525451 ,  4.3526074,
	6   4.3526962 ,  4.3527618 ,  4.3528509 ,  4.3529514/
c 
	data (p1u65b(i,4),i=1,ntf4)  / 370.7158980 ,323.8113984,
	1 282.2017033 ,243.3854523 ,220.4653358 ,193.6221695 ,167.3204770,
	2 153.7400442 ,142.0123688 ,131.4618543 ,125.4756325 ,118.7764418,
	3 112.6490260 ,109.7082704 ,107.3234176 ,105.3252945 ,104.2644683,
	4 103.1500223 ,102.2084586 ,101.7871276 ,101.4621261 ,101.2027664,
	5 101.0704473 ,100.9359712 ,100.8264729 ,100.7788872 ,100.7428725,
	6 100.7146249 ,100.7004064 ,100.6861097 ,100.6746014/
c 
	data (p2u65b(i,4),i=1,ntf4)  / 322.5834814 ,283.3802192,
	1 248.6631304 ,216.3581709 ,197.3395656 ,175.1484622 ,153.5408672,
	2 142.4666293 ,132.9710437 ,124.5039281 ,119.7428672 ,114.4643620,
	3 109.6975253 ,107.4373044 ,105.6212318 ,104.1139621 ,103.3201364,
	4 102.4919376 ,101.7977418 ,101.4890966 ,101.2520330 ,101.0635967,
	5 100.9677605 ,100.8706067 ,100.7917159 ,100.7575060 ,100.7316511,
	6 100.7113988 ,100.7012153 ,100.6909846 ,100.6827571/
c 
c j =  5, TR spacing =     40 
c 
	data (psu65b(i,5),i=1,ntf4)  /  48.3682673 , 40.6906641,
	1  33.8222952 , 27.3368772 , 23.4520607 , 18.8199247 , 14.1419128,
	2  11.6388272 ,  9.4022349 ,  7.3033673 ,  6.0608156 ,  4.6079926,
	3   3.1984401 ,  2.4836207 ,  1.8792509 ,  1.3511156 ,  1.0606257,
	4   0.7460755 ,  0.4710359 ,  0.3445187 ,  0.2451355 ,  0.1644883,
	5   0.1228015 ,  0.0799874 ,  0.0447247 ,  0.0292623 ,  0.0174915,
	6   0.0082098 ,  0.0035180 , -0.0012160 , -0.0050415/
c 
	data (adu65b(i,5),i=1,ntf4)  /   9.9630521 ,  8.8336572,
	1   7.8337928 ,  6.9039778 ,  6.3570699 ,  5.7198071 ,  5.1012313,
	2   4.7858821 ,  4.5174626 ,  4.2811294 ,  4.1505231 ,  4.0091823,
	3   3.8871645 ,  3.8326209 ,  3.7913729 ,  3.7598006 ,  3.7445977,
	4   3.7302590 ,  3.7199807 ,  3.7161614 ,  3.7136807 ,  3.7120955,
	5   3.7114699 ,  3.7110100 ,  3.7108255 ,  3.7108270 ,  3.7108802,
	6   3.7109700 ,  3.7110403 ,  3.7111383 ,  3.7112514/
c 
	data (p1u65b(i,5),i=1,ntf4)  / 370.4867709 ,323.6027002,
	1 282.0166777 ,243.2289075 ,220.3294479 ,193.5146316 ,167.2444192,
	2 153.6808451 ,141.9674151 ,131.4290208 ,125.4492466 ,118.7568473,
	3 112.6352702 ,109.6972184 ,107.3145322 ,105.3182220 ,104.2583636,
	4 103.1449435 ,102.2042586 ,101.7833262 ,101.4586351 ,101.1995255,
	5 101.0673351 ,100.9329905 ,100.8236001 ,100.7760616 ,100.7400827,
	6 100.7118634 ,100.6976591 ,100.6833768 ,100.6718801/
c 
	data (p2u65b(i,5),i=1,ntf4)  / 418.8550382 ,364.2933643,
	1 315.8389729 ,270.5657847 ,243.7815086 ,212.3345563 ,181.3863320,
	2 165.3196723 ,151.3696500 ,138.7323881 ,131.5100622 ,123.3648399,
	3 115.8337103 ,112.1808391 ,109.1937831 ,106.6693376 ,105.3189893,
	4 103.8910190 ,102.6752945 ,102.1278449 ,101.7037706 ,101.3640138,
	5 101.1901366 ,101.0129779 ,100.8683248 ,100.8053239 ,100.7575742,
	6 100.7200732 ,100.7011771 ,100.6821608 ,100.6668386/
c 