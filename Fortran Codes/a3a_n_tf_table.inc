	include 'a3a_tf_table_def.inc'
c 
	data (p1_acal_3a(i),i=1,ntr)  /  96.1320426 , 96.0505248,
	1  96.0528365 , 96.0486394 , 96.0379504/
c 
	data (p2_acal_3a(i),i=1,ntr)  /  96.0735117 , 96.0802105,
	1  96.0367360 , 96.0584912 , 96.0318596/
c 
	data (ps_acal_3a(i),i=1,ntr)  /  -0.0585309 , -0.0296857,
	1  -0.0161005 , -0.0098518 , -0.0060908/
c 
	data (ad_acal_3a(i),i=1,ntr)  /  15.7937038 , 10.1807502,
	1   7.5510470 ,  6.1556709 ,  5.0885746/
c 
	data (ps_offset_3a(i),i=1,ntr)  /   0.0490714 ,  0.0152386,
	1   0.0000000 ,  0.0000000 ,  0.0000000/
c 
	data (ru3a(i),i=1,ntf2)      /   0.1000000 ,  0.1500000,
	1   0.2000000 ,  0.3000000 ,  0.5000000 ,  0.7000000 ,  1.0000000,
	2   1.5000000 ,  2.0000000 ,  3.0000000 ,  5.0000000 ,  7.0000000,
	3  10.0000000 , 15.0000000 , 20.0000000 , 30.0000000 , 50.0000000,
	4  70.0000000 ,100.0000000 ,150.0000000 ,200.0000000 ,300.0000000,
	5 500.0000000 ,700.0000000 ,1000.0000000 ,1500.0000000 ,2000.0000000,
	6 3000.0000000 ,5000.0000000/
c 
c j =  1, TR spacing =     10 
c 
	data (psu3a(i,1),i=1,ntf2)  /  81.0137409 , 64.6048541,
	1  54.7897972 , 43.1271829 , 31.4564373 , 25.3030132 , 19.9001107,
	2  14.9558386 , 12.1103409 ,  8.8878212 ,  5.9025704 ,  4.4564339,
	3   3.2771067 ,  2.2843860 ,  1.7555478 ,  1.1985234 ,  0.7265493,
	4   0.5148298 ,  0.3512109 ,  0.2203707 ,  0.1534875 ,  0.0853786,
	5   0.0297643 ,  0.0055231 , -0.0128733 , -0.0273531 , -0.0346712,
	6  -0.0420664 , -0.0480714/
c 
	data (adu3a(i,1),i=1,ntf2)  /  24.4558264 , 22.2843799,
	1  21.0314072 , 19.6054880 , 18.2746565 , 17.6267663 , 17.0988672,
	2  16.6584714 , 16.4280482 , 16.1923501 , 16.0039883 , 15.9259269,
	3  15.8703393 , 15.8305233 , 15.8125550 , 15.7967317 , 15.7866134,
	4  15.7834426 , 15.7818274 , 15.7812887 , 15.7813948 , 15.7819161,
	5  15.7828589 , 15.7835328 , 15.7842389 , 15.7850000 , 15.7855032,
	6  15.7861548 , 15.7868784/
c 
	data (p1u3a(i,1),i=1,ntf2)  / 162.1564638 ,146.3606290,
	1 137.2347751 ,126.8182547 ,117.0164643 ,112.1766469 ,108.1609469,
	2 104.7135599 ,102.8436278 ,100.8411788 , 99.1106060 , 98.3214232,
	3  97.7047873 , 97.2067026 , 96.9501205 , 96.6872858 , 96.4713916,
	4  96.3769142 , 96.3050826 , 96.2485035 , 96.2199276 , 96.1911142,
	5  96.1678471 , 96.1577981 , 96.1502206 , 96.1442947 , 96.1413172,
	6  96.1383252 , 96.1359152/
c 
	data (p2u3a(i,1),i=1,ntf2)  / 243.1702047 ,210.9654831,
	1 192.0245723 ,169.9454376 ,148.4729016 ,137.4796601 ,128.0610576,
	2 119.6693985 ,114.9539687 ,109.7290000 ,105.0131764 ,102.7778571,
	3 100.9818940 , 99.4910886 , 98.7056683 , 97.8858092 , 97.1979408,
	4  96.8917440 , 96.6562936 , 96.4688742 , 96.3734150 , 96.2764928,
	5  96.1976113 , 96.1633212 , 96.1373473 , 96.1169416 , 96.1066460,
	6  96.0962589 , 96.0878438/
c 
c j =  2, TR spacing =     16 
c 
	data (psu3a(i,2),i=1,ntf2)  /  83.8186741 , 67.6096562,
	1  57.9064480 , 46.3443699 , 34.6819234 , 28.4528682 , 22.8993068,
	2  17.7052416 , 14.6425290 , 11.0783951 ,  7.6448281 ,  5.9159427,
	3   4.4623602 ,  3.1992979 ,  2.5077712 ,  1.7617189 ,  1.1117470,
	4   0.8134323 ,  0.5793017 ,  0.3893373 ,  0.2910924 ,  0.1900796,
	5   0.1066976 ,  0.0700275 ,  0.0420259 ,  0.0198486 ,  0.0085780,
	6  -0.0028719 , -0.0122386/
c 
	data (adu3a(i,2),i=1,ntf2)  /  20.1661879 , 17.8734652,
	1  16.5282399 , 14.9654708 , 13.4562936 , 12.6917221 , 12.0449455,
	2  11.4799629 , 11.1702910 , 10.8380776 , 10.5547768 , 10.4296974,
	3  10.3360870 , 10.2653258 , 10.2317196 , 10.2005678 , 10.1789999,
	4  10.1714465 , 10.1669484 , 10.1645746 , 10.1639823 , 10.1640478,
	5  10.1649290 , 10.1657287 , 10.1666404 , 10.1676774 , 10.1683851,
	6  10.1693211 , 10.1703806/
c 
	data (p1u3a(i,2),i=1,ntf2)  / 318.1774907 ,271.5001798,
	1 243.8954952 ,211.4918529 ,179.6019192 ,163.0416752 ,148.6565797,
	2 135.6155291 ,128.1555711 ,119.7326973 ,111.9302787 ,108.1377789,
	3 105.0305437 ,102.3984405 ,100.9872545 , 99.4914242 , 98.2137591,
	4  97.6365325 , 97.1881886 , 96.8279008 , 96.6429822 , 96.4540327,
	5  96.2991424 , 96.2314089 , 96.1798891 , 96.1392443 , 96.1186601,
	6  96.0978178 , 96.0808464/
c 
	data (p2u3a(i,2),i=1,ntf2)  / 234.3588166 ,203.8905236,
	1 185.9890472 ,165.1474830 ,144.9199958 ,134.5888070 ,125.7572729,
	2 117.9102875 ,113.5130421 ,108.6543022 ,104.2854506 ,102.2218362,
	3 100.5681835 , 99.1991426 , 98.4794833 , 97.7297053 , 97.1020121,
	4  96.8231002 , 96.6088869 , 96.4385635 , 96.3518898 , 96.2639531,
	5  96.1924448 , 96.1613814 , 96.1378631 , 96.1193957 , 96.1100820,
	6  96.1006897 , 96.0930851/
c 
c j =  3, TR spacing =     22 
c 
	data (psu3a(i,3),i=1,ntf2)  /  85.0103444 , 68.9045278,
	1  59.2708137 , 47.7969912 , 36.2179877 , 30.0184389 , 24.4673903,
	2  19.2351435 , 16.1179639 , 12.4413256 ,  8.8196635 ,  6.9503336,
	3   5.3442407 ,  3.9142222 ,  3.1135945 ,  2.2317442 ,  1.4438358,
	4   1.0743293 ,  0.7799289 ,  0.5375942 ,  0.4107853 ,  0.2791193,
	5   0.1692183 ,  0.1204402 ,  0.0829550 ,  0.0530769 ,  0.0378071,
	6   0.0222107 ,  0.0093572/
c 
	data (adu3a(i,3),i=1,ntf2)  /  18.2703624 , 15.9283309,
	1  14.5451563 , 12.9247182 , 11.3364171 , 10.5166319 ,  9.8098898,
	2   9.1769228 ,  8.8205357 ,  8.4267072 ,  8.0760045 ,  7.9140765,
	3   7.7883105 ,  7.6892396 ,  7.6403124 ,  7.5931887 ,  7.5587407,
	4   7.5458827 ,  7.5376725 ,  7.5327344 ,  7.5310521 ,  7.5302544,
	5   7.5307387 ,  7.5315193 ,  7.5325278 ,  7.5337565 ,  7.5346271,
	6   7.5358063 ,  7.5371694/
c 
	data (p1u3a(i,3),i=1,ntf2)  / 309.1054755 ,264.1593490,
	1 237.5908599 ,206.4220568 ,175.7788081 ,159.8858583 ,146.0975891,
	2 133.6177221 ,126.4909154 ,118.4590344 ,111.0383687 ,107.4408550,
	3 104.4995370 ,102.0135044 ,100.6832332 , 99.2756279 , 98.0758024,
	4  97.5346820 , 97.1148829 , 96.7779155 , 96.6051250 , 96.4287031,
	5  96.2842078 , 96.2210656 , 96.1730620 , 96.1352104 , 96.1160495,
	6  96.0966568 , 96.0808753/
c 
	data (p2u3a(i,3),i=1,ntf2)  / 394.1158199 ,333.0638768,
	1 296.8616736 ,254.2190480 ,211.9967958 ,189.9042972 ,170.5649794,
	2 152.8528656 ,142.6088793 ,130.9003600 ,119.8580322 ,114.3911886,
	3 109.8437777 ,105.9277266 ,103.7968277 ,101.5073721 , 99.5196382,
	4  98.6090112 , 97.8948118 , 97.3155096 , 97.0159103 , 96.7078223,
	5  96.4534261 , 96.3415057 , 96.2560170 , 96.1882874 , 96.1538565,
	6  96.1188675 , 96.0902325/
c 
c j =  4, TR spacing =     28 
c 
	data (psu3a(i,4),i=1,ntf2)  /  85.5491502 , 69.4978945,
	1  59.9048257 , 48.4881195 , 36.9747457 , 30.8103523 , 25.2857452,
	2  20.0657775 , 16.9437057 , 13.2393414 ,  9.5492605 ,  7.6179652,
	3   5.9363739 ,  4.4148473 ,  3.5495707 ,  2.5818050 ,  1.7000826,
	4   1.2793367 ,  0.9398865 ,  0.6570173 ,  0.5074873 ,  0.3508884,
	5   0.2188818 ,  0.1598119 ,  0.1141585 ,  0.0775623 ,  0.0587648,
	6   0.0394740 ,  0.0234725/
c 
	data (adu3a(i,4),i=1,ntf2)  /  17.2894735 , 14.9260923,
	1  13.5261179 , 11.8794856 , 10.2539861 ,  9.4073227 ,  8.6703765,
	2   8.0016190 ,  7.6194821 ,  7.1899168 ,  6.7971908 ,  6.6106057,
	3   6.4620476 ,  6.3416322 ,  6.2804993 ,  6.2199885 ,  6.1740430,
	4   6.1561673 ,  6.1442816 ,  6.1366650 ,  6.1337737 ,  6.1319378,
	5   6.1318276 ,  6.1324796 ,  6.1334869 ,  6.1348155 ,  6.1357944,
	6   6.1371523 ,  6.1387539/
c 
	data (p1u3a(i,4),i=1,ntf2)  / 461.4105341 ,387.7498733,
	1 344.0133926 ,292.4052088 ,241.1420227 ,214.2088059 ,190.5307601,
	2 168.7182124 ,156.0208122 ,141.3995127 ,127.4520603 ,120.4623019,
	3 114.5867599 ,109.4666264 ,106.6496466 ,103.5913655 ,100.9014220,
	4  99.6549659 , 98.6693890 , 97.8635857 , 97.4440966 , 97.0103113,
	5  96.6498220 , 96.4903771 , 96.3681310 , 96.2709143 , 96.2213277,
	6  96.1707755 , 96.1292201/
c 
	data (p2u3a(i,4),i=1,ntf2)  / 375.8613839 ,318.2519788,
	1 284.1085669 ,243.9170893 ,204.1672770 ,183.3984536 ,165.2450149,
	2 148.6524349 ,139.0771065 ,128.1601713 ,117.9027998 ,112.8443367,
	3 108.6503860 ,105.0517791 ,103.1000759 ,101.0095605 , 99.2013394,
	4  98.3756292 , 97.7295025 , 97.2065684 , 96.9366093 , 96.6594229,
	5  96.4309402 , 96.3305652 , 96.2539725 , 96.1933520 , 96.1625629,
	6  96.1313015 , 96.1057476/
c 
c j =  5, TR spacing =     34 
c 
	data (psu3a(i,5),i=1,ntf2)  /  85.9498770 , 69.9321751,
	1  60.3640124 , 48.9844110 , 37.5195614 , 31.3857636 , 25.8894931,
	2  20.6928396 , 17.5792822 , 13.8728221 , 10.1539554 ,  8.1880970,
	3   6.4584762 ,  4.8722823 ,  3.9576860 ,  2.9200613 ,  1.9565036,
	4   1.4885270 ,  1.1059682 ,  0.7829347 ,  0.6102633 ,  0.4276892,
	5   0.2720698 ,  0.2017855 ,  0.1471110 ,  0.1029977 ,  0.0802093,
	6   0.0566962 ,  0.0370500/
c 
	data (adu3a(i,5),i=1,ntf2)  /  16.5559031 , 14.1766099,
	1  12.7645669 , 11.0996343 ,  9.4484637 ,  8.5830980 ,  7.8247703,
	2   7.1299973 ,  6.7285836 ,  6.2713262 ,  5.8443183 ,  5.6365353,
	3   5.4674828 ,  5.3268895 ,  5.2536751 ,  5.1793214 ,  5.1208113,
	4   5.0971669 ,  5.0808870 ,  5.0699257 ,  5.0654559 ,  5.0621914,
	5   5.0611694 ,  5.0615619 ,  5.0624772 ,  5.0638450 ,  5.0649083,
	6   5.0664293 ,  5.0682687/
c 
	data (p1u3a(i,5),i=1,ntf2)  / 452.2119879 ,380.2681449,
	1 337.5573942 ,287.1694916 ,237.1362786 ,210.8613751 ,187.7734144,
	2 166.5189281 ,154.1558396 ,139.9321359 ,126.3824515 ,119.6020030,
	3 113.9097035 ,108.9565275 ,106.2351716 ,103.2846243 ,100.6937910,
	4  99.4950510 , 98.5482245 , 97.7749246 , 97.3727118 , 96.9571046,
	5  96.6120218 , 96.4595023 , 96.3426255 , 96.2497267 , 96.2023642,
	6  96.1541006 , 96.1144505/
c 
	data (p2u3a(i,5),i=1,ntf2)  / 538.1618649 ,450.2003199,
	1 397.9214066 ,336.1539026 ,274.6558400 ,242.2471387 ,213.6629075,
	2 187.2117677 ,171.7351218 ,153.8049580 ,136.5364069 ,127.7901000,
	3 120.3681797 ,113.8288098 ,110.1928576 ,106.2046856 ,102.6502946,
	4 100.9835780 , 99.6541927 , 98.5578593 , 97.9829751 , 97.3847938,
	5  96.8840916 , 96.6612877 , 96.4897365 , 96.3527244 , 96.2825735,
	6  96.2107968 , 96.1515005/
c 