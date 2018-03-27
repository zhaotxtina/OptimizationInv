	  parameter (nrm=19,nrt=25,nhd=28,ntr=5,ntrm=7,necc=67)
	  integer necc_hd(nhd)    ! number of ecc for a given hd
	  real*4 rt(nrt),hd(nhd),rm(nrm),ecc(nhd,necc)  !grid for rt, hd, rm and ecc, ecc grid are different for different HD
	  real*4 ps_a(nrt,nhd,nrm,ntr,necc),ad_a(nrt,nhd,nrm,ntr,necc)
	  real*4 ps_b(nrt,nhd,nrm,ntr,necc),ad_b(nrt,nhd,nrm,ntr,necc)
	  integer nrt_t,nrm_t,nhd_t,necc_t
	  integer nrm,nrt,nhd,ntr,ntrm,necc
      common /arc_bhctab_grid/ rt,hd,rm,ecc,nrm_t,nrt_t,nhd_t,necc_t,necc_hd
	  common /arc_bhctab_ps_a/ ps_a   !2 MHZ unmixed attenuation
	  common /arc_bhctab_ad_a/ ad_a   !2 MHZ unmixed phaseshift
	  common /arc_bhctab_ps_b/ ps_b   !400k MHZ unmixed attenuation
	  common /arc_bhctab_ad_b/ ad_b   !400k MHZ unmixed phaseshift

