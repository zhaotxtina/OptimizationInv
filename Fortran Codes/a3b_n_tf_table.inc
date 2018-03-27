	include 'a3b_tf_table_def.inc'
c 
	data (p1_acal_3b(i),i=1,ntr)  / 102.5646853 ,102.3804252,
	1 102.3855812 ,102.3458361 ,102.3407849/
c 
	data (p2_acal_3b(i),i=1,ntr)  / 102.4318169 ,102.4470837,
	1 102.3491911 ,102.3703758 ,102.3254190/
c 
	data (ps_acal_3b(i),i=1,ntr)  /  -0.1328684 , -0.0666585,
	1  -0.0363901 , -0.0245397 , -0.0153659/
c 
	data (ad_acal_3b(i),i=1,ntr)  /  15.8049063 , 10.1872340,
	1   7.5554761 ,  6.1593852 ,  5.0919664/
c 
	data (ps_offset_3b(i),i=1,ntr)  /   0.1312910 ,  0.0644251,
	1   0.0333856 ,  0.0209706 ,  0.0113048/
c 
	data (ru3b(i),i=1,ntf4)      /   0.0500000 ,  0.0700000,
	1   0.1000000 ,  0.1500000 ,  0.2000000 ,  0.3000000 ,  0.5000000,
	2   0.7000000 ,  1.0000000 ,  1.5000000 ,  2.0000000 ,  3.0000000,
	3   5.0000000 ,  7.0000000 , 10.0000000 , 15.0000000 , 20.0000000,
	4  30.0000000 , 50.0000000 , 70.0000000 ,100.0000000 ,150.0000000,
	5 200.0000000 ,300.0000000 ,500.0000000 ,700.0000000 ,1000.0000000,
	6 1500.0000000 ,2000.0000000 ,3000.0000000 ,5000.0000000/
c 
c j =  1, TR spacing =     10 
c 
	data (psu3b(i,1),i=1,ntf4)  /  48.7744647 , 39.8380130,
	1  31.8936518 , 24.4943195 , 20.1516053 , 15.1247700 , 10.3199512,
	2   7.9198220 ,  5.9150368 ,  4.1853662 ,  3.2444638 ,  2.2354880,
	3   1.3631446 ,  0.9655635 ,  0.6551902 ,  0.4047947 ,  0.2759681,
	4   0.1441674 ,  0.0360856 , -0.0111369 , -0.0469956 , -0.0752020,
	5  -0.0894302 , -0.1037591 , -0.1153110 , -0.1202920 , -0.1240426,
	6  -0.1269705 , -0.1284387 , -0.1299104 , -0.1310910/
c 
	data (adu3b(i,1),i=1,ntf4)  /  20.3294811 , 19.2631151,
	1  18.3667115 , 17.5895691 , 17.1667515 , 16.7164327 , 16.3358614,
	2  16.1690483 , 16.0446582 , 15.9506742 , 15.9058363 , 15.8637642,
	3  15.8335907 , 15.8222517 , 15.8147272 , 15.8097054 , 15.8075739,
	4  15.8057945 , 15.8047286 , 15.8044146 , 15.8042621 , 15.8042142,
	5  15.8042242 , 15.8042684 , 15.8043434 , 15.8043942 , 15.8044451,
	6  15.8044973 , 15.8045303 , 15.8045712 , 15.8046142/
c 
	data (p1u3b(i,1),i=1,ntf4)  / 138.3160419 ,130.5008212,
	1 123.8831198 ,118.0697234 ,114.8483029 ,111.3285007 ,108.2100638,
	2 106.7570996 ,105.6046458 ,104.6603738 ,104.1683573 ,103.6596657,
	3 103.2376244 ,103.0515214 ,102.9093561 ,102.7969205 ,102.7399655,
	4 102.6824153 ,102.6358529 ,102.6157211 ,102.6005364 ,102.5886648,
	5 102.5827050 ,102.5767258 ,102.5719255 ,102.5698625 ,102.5683125,
	6 102.5671049 ,102.5665003 ,102.5658951 ,102.5654104/
c 
	data (p2u3b(i,1),i=1,ntf4)  / 187.0905066 ,170.3388342,
	1 155.7767716 ,142.5640429 ,134.9999082 ,126.4532707 ,118.5300150,
	2 114.6769216 ,111.5196826 ,108.8457400 ,107.4128211 ,105.8951537,
	3 104.6007690 ,104.0170849 ,103.5645463 ,103.2017152 ,103.0159336,
	4 102.8265827 ,102.6719385 ,102.6045842 ,102.5535408 ,102.5134628,
	5 102.4932748 ,102.4729667 ,102.4566145 ,102.4495705 ,102.4442699,
	6 102.4401344 ,102.4380616 ,102.4359847 ,102.4343194/
c 
c j =  2, TR spacing =     16 
c 
	data (psu3b(i,2),i=1,ntf4)  /  51.9765779 , 43.0967785,
	1  35.1467362 , 27.6517265 , 23.1839447 , 17.9090184 , 12.7027544,
	2  10.0091283 ,  7.6898962 ,  5.6201746 ,  4.4591876 ,  3.1785812,
	3   2.0331631 ,  1.4959695 ,  1.0682999 ,  0.7168499 ,  0.5333569,
	4   0.3433772 ,  0.1855445 ,  0.1158803 ,  0.0626324 ,  0.0204975,
	5  -0.0008555 , -0.0224396 , -0.0399111 , -0.0474684 , -0.0531709,
	6  -0.0576309 , -0.0598709 , -0.0621191 , -0.0639251/
c 
	data (adu3b(i,2),i=1,ntf4)  /  15.7584281 , 14.5763957,
	1  13.5550999 , 12.6375258 , 12.1190054 , 11.5435123 , 11.0274514,
	2  10.7872392 , 10.5992118 , 10.4495374 , 10.3746917 , 10.3014136,
	3  10.2460522 , 10.2242527 , 10.2092829 , 10.1989228 , 10.1943685,
	4  10.1904232 , 10.1879068 , 10.1870885 , 10.1866265 , 10.1863956,
	5  10.1863407 , 10.1863458 , 10.1864177 , 10.1864795 , 10.1865467,
	6  10.1866192 , 10.1866664 , 10.1867260 , 10.1867895/
c 
	data (p1u3b(i,2),i=1,ntf4)  / 233.6708420 ,208.9889900,
	1 187.3250845 ,167.4136790 ,155.8538566 ,142.5856297 ,129.9954942,
	2 123.7238563 ,118.4811474 ,113.9429354 ,111.4626266 ,108.7879476,
	3 106.4568287 ,105.3861793 ,104.5455496 ,103.8634761 ,103.5108852,
	4 103.1487170 ,102.8503977 ,102.7195935 ,102.6200350 ,102.5415549,
	5 102.5019009 ,102.4619121 ,102.4296258 ,102.4156883 ,102.4051855,
	6 102.3969810 ,102.3928644 ,102.3887360 ,102.3854226/
c 
	data (p2u3b(i,2),i=1,ntf4)  / 181.6942641 ,165.8922115,
	1 152.1783483 ,139.7619525 ,132.6699119 ,124.6766113 ,117.2927398,
	2 113.7147280 ,110.7912512 ,108.3227608 ,107.0034390 ,105.6093664,
	3 104.4236656 ,103.8902098 ,103.4772497 ,103.1466262 ,102.9775283,
	4 102.8053398 ,102.6648532 ,102.6037132 ,102.5574026 ,102.5210574,
	5 102.5027564 ,102.4843517 ,102.4695369 ,102.4631567 ,102.4583564,
	6 102.4546119 ,102.4527353 ,102.4508551 ,102.4493477/
c 
c j =  3, TR spacing =     22 
c 
	data (psu3b(i,3),i=1,ntf4)  /  53.3967728 , 44.5842306,
	1  36.6888329 , 29.2267960 , 24.7589208 , 19.4461019 , 14.1280324,
	2  11.3262145 ,  8.8702066 ,  6.6293456 ,  5.3445353 ,  3.8961288,
	3   2.5637447 ,  1.9230015 ,  1.4036040 ,  0.9692211 ,  0.7391629,
	4   0.4981315 ,  0.2952375 ,  0.2047504 ,  0.1351192 ,  0.0796808,
	5   0.0514509 ,  0.0228058 , -0.0004782 , -0.0105830 , -0.0182237,
	6  -0.0242116 , -0.0272236 , -0.0302506 , -0.0326856/
c 
	data (adu3b(i,3),i=1,ntf4)  /  13.7465258 , 12.5150329,
	1  11.4380015 , 10.4538828 ,  9.8870147 ,  9.2435998 ,  8.6461661,
	2   8.3572843 ,  8.1235313 ,  7.9302660 ,  7.8300898 ,  7.7286097,
	3   7.6485526 ,  7.6157584 ,  7.5925755 ,  7.5760445 ,  7.5685767,
	4   7.5619366 ,  7.5575355 ,  7.5560323 ,  7.5551325 ,  7.5546261,
	5   7.5544640 ,  7.5543915 ,  7.5544333 ,  7.5544961 ,  7.5545728,
	6   7.5546610 ,  7.5547203 ,  7.5547967 ,  7.5548795/
c 
	data (p1u3b(i,3),i=1,ntf4)  / 227.9986520 ,204.2653938,
	1 183.4515344 ,164.3432073 ,153.2636672 ,140.5652704 ,128.5427235,
	2 122.5680704 ,117.5838046 ,113.2791698 ,110.9314794 ,108.4048013,
	3 106.2080097 ,105.2011665 ,104.4117921 ,103.7722025 ,103.4419465,
	4 103.1030343 ,102.8241566 ,102.7019756 ,102.6090290 ,102.5357960,
	5 102.4988071 ,102.4615170 ,102.4314194 ,102.4184301 ,102.4086435,
	6 102.4009996 ,102.3971648 ,102.3933194 ,102.3902335/
c 
	data (p2u3b(i,3),i=1,ntf4)  / 281.3954248 ,248.8496244,
	1 220.1403673 ,193.5700033 ,178.0225880 ,160.0113723 ,142.6707559,
	2 133.8942849 ,126.4540112 ,119.9085154 ,116.2760147 ,112.3009301,
	3 108.7717544 ,107.1241680 ,105.8153961 ,104.7414236 ,104.1811094,
	4 103.6011658 ,103.1193941 ,102.9067260 ,102.7441482 ,102.6154768,
	5 102.5502580 ,102.4843228 ,102.4309412 ,102.4078471 ,102.3904198,
	6 102.3767880 ,102.3699412 ,102.3630688 ,102.3575479/
c 
c j =  4, TR spacing =     28 
c 
	data (psu3b(i,4),i=1,ntf4)  /  54.0622589 , 45.2957074,
	1  37.4447324 , 30.0236182 , 25.5758551 , 20.2748910 , 14.9389067,
	2  12.1039119 ,  9.5954165 ,  7.2773783 ,  5.9301702 ,  4.3892269,
	3   2.9431831 ,  2.2345734 ,  1.6519415 ,  1.1576664 ,  0.8927437,
	4   0.6123588 ,  0.3736263 ,  0.2661790 ,  0.1829994 ,  0.1164093,
	5   0.0823548 ,  0.0476796 ,  0.0193882 ,  0.0070740 , -0.0022553,
	6  -0.0095794 , -0.0132686 , -0.0169805 , -0.0199706/
c 
	data (adu3b(i,4),i=1,ntf4)  /  12.7145572 , 11.4603189,
	1  10.3569791 ,  9.3404234 ,  8.7491711 ,  8.0700898 ,  7.4272123,
	2   7.1093554 ,  6.8468237 ,  6.6243556 ,  6.5062088 ,  6.3836025,
	3   6.2837574 ,  6.2416178 ,  6.2111536 ,  6.1889216 ,  6.1786676,
	4   6.1693728 ,  6.1630477 ,  6.1608216 ,  6.1594465 ,  6.1586305,
	5   6.1583423 ,  6.1581719 ,  6.1581659 ,  6.1582215 ,  6.1583011,
	6   6.1583992 ,  6.1584674 ,  6.1585570 ,  6.1586556/
c 
	data (p1u3b(i,4),i=1,ntf4)  / 323.9623695 ,284.5313008,
	1 249.6559331 ,217.2568073 ,198.2155971 ,176.0394461 ,154.5040343,
	2 143.4968996 ,134.0807935 ,125.7067742 ,121.0098357 ,115.8150759,
	3 111.1384948 ,108.9273074 ,107.1542944 ,105.6857637 ,104.9136491,
	4 104.1092507 ,103.4360862 ,103.1371757 ,102.9077821 ,102.7255828,
	5 102.6329744 ,102.5391379 ,102.4629807 ,102.4299694 ,102.4050270,
	6 102.3854940 ,102.3756741 ,102.3658100 ,102.3578786/
c 
	data (p2u3b(i,4),i=1,ntf4)  / 269.9001106 ,239.2355934,
	1 212.2112007 ,187.2331891 ,172.6397420 ,155.7645551 ,139.5651276,
	2 131.3929877 ,124.4853770 ,118.4293959 ,115.0796655 ,111.4258490,
	3 108.1953117 ,106.6927340 ,105.5023529 ,104.5280973 ,104.0209054,
	4 103.4968919 ,103.0624599 ,102.8709967 ,102.7247827 ,102.6091735,
	5 102.5506196 ,102.4914583 ,102.4435925 ,102.4228954 ,102.4072823,
	6 102.3950734 ,102.3889427 ,102.3827905 ,102.3778492/
c 
c j =  5, TR spacing =     34 
c 
	data (psu3b(i,5),i=1,ntf4)  /  54.5424955 , 45.8070698,
	1  37.9895474 , 30.6043665 , 26.1786180 , 20.9002860 , 15.5727442,
	2  12.7281291 , 10.1952275 ,  7.8326464 ,  6.4447552 ,  4.8376902,
	3   3.3023285 ,  2.5364000 ,  1.8976791 ,  1.3478022 ,  1.0493367,
	4   0.7299670 ,  0.4545915 ,  0.3293772 ,  0.2317823 ,  0.1531601,
	5   0.1127533 ,  0.0714464 ,  0.0375975 ,  0.0228142 ,  0.0115895,
	6   0.0027595 , -0.0016956 , -0.0061840 , -0.0098048/
c 
	data (adu3b(i,5),i=1,ntf4)  /  11.9437614 , 10.6737824,
	1   9.5523124 ,  8.5131867 ,  7.9046697 ,  7.1997123 ,  6.5225118,
	2   6.1817893 ,  5.8956163 ,  5.6480013 ,  5.5136695 ,  5.3711697,
	3   5.2515920 ,  5.1996383 ,  5.1612296 ,  5.1325356 ,  5.1190200,
	4   5.1065302 ,  5.0978121 ,  5.0946598 ,  5.0926621 ,  5.0914292,
	5   5.0909660 ,  5.0906541 ,  5.0905717 ,  5.0906103 ,  5.0906873,
	6   5.0907924 ,  5.0908688 ,  5.0909716 ,  5.0910870/
c 
	data (p1u3b(i,5),i=1,ntf4)  / 318.1304649 ,279.6369320,
	1 245.6009204 ,213.9950389 ,195.4293575 ,173.8204545 ,152.8573336,
	2 142.1552545 ,133.0100445 ,124.8876097 ,120.3377069 ,115.3121933,
	3 110.7958999 ,108.6639759 ,106.9566013 ,105.5441660 ,104.8023065,
	4 104.0301026 ,103.3845196 ,103.0980847 ,102.8783808 ,102.7039622,
	5 102.6153423 ,102.5255748 ,102.4527444 ,102.4211836 ,102.3973412,
	6 102.3786726 ,102.3692884 ,102.3598631 ,102.3522854/
c 
	data (p2u3b(i,5),i=1,ntf4)  / 372.6729604 ,325.4440018,
	1 283.5904678 ,244.5994054 ,221.6079755 ,194.7207405 ,168.4300778,
	2 154.8833836 ,143.2052720 ,132.7202561 ,126.7824621 ,120.1498835,
	3 114.0982284 ,111.2003759 ,108.8542804 ,106.8919682 ,105.8516432,
	4 104.7600696 ,103.8391111 ,103.4274619 ,103.1101631 ,102.8571223,
	5 102.7280956 ,102.5970212 ,102.4903419 ,102.4439978 ,102.4089307,
	6 102.3814321 ,102.3675928 ,102.3536791 ,102.3424806/
c 