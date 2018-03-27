	  parameter (nrm=19,nrt=26,nhd=19,ntr=5,ntrm=7,nepsr=27)
	  real*4 rt(nrt),hd(nhd),rm(nrm),epsr(nepsr)  !grid for rt, hd, rm and epsr
	  real*4 ps_a(nrt,nhd,nrm,ntr,nepsr),ad_a(nrt,nhd,nrm,ntr,nepsr)
	  real*4 ps_b(nrt,nhd,nrm,ntr,nepsr),ad_b(nrt,nhd,nrm,ntr,nepsr)
	  integer nrt_t,nrm_t,nhd_t,nepsr_t    !  what are those variables for?
	  integer nrm,nrt,nhd,ntr,ntrm,nepsr
      common /arc_bhctab_grid/ rt,hd,rm,epsr,nrm_t,nrt_t,nhd_t,nepsr_t
	  common /arc_bhctab_ps_a/ ps_a   !2 MHZ unmixed attenuation
	  common /arc_bhctab_ad_a/ ad_a   !2 MHZ unmixed phaseshift
	  common /arc_bhctab_ps_b/ ps_b   !400k MHZ unmixed attenuation
	  common /arc_bhctab_ad_b/ ad_b   !400k MHZ unmixed phaseshift

