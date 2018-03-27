
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	subroutine arc_bhcinv2_sub(rto,hdo,erro,id_tool,ibadtx1,
c	1	i_init,bsi,rmi,psm2i,adm2i,wps2,wad2,ndepth)
C
C To Do 2-parameter HD, Rt inversion for ARC 2MHz data using TTM model
c
C Inputs:	
c	id_tool	- integer switch =5, 6, 8, 9 for ARC5, 6, 8, 9 tool, respectively 
c	ibadtx1	- integer flag =1,2,3,4,5,6 for NG tx #, 6=all good 
c	i_init	- integer switch =0 initialize the routine for a new set of input para
c				=1 subsequent run through different depth frame 
c	ndepth	- integer for no. of depth frames in the input arrays
c				for level-by-level onperation, set ndepth=1
c	bsi		- 1xndepth array for bit size (in) 
c	
c	psm2i	- 5xndepth array for 2MHz air_cal_corrected mbhc phase-shift (deg) 
c	adm2i	- 5xndepth array for 2MHz air_cal_corrected mbhc attenuation (dB)
c	wps2	- 5x1 array for the weighting factor for psm2 inputs
c			default should be (1 1 1 1 1) for selecting all
c			with equal weights of one.
c	wad2	- 7x1 array for the weighting factor for adm2 inputs
c			default should be (1 1 1 1 1) for selecting all
c			with equal weights of one.
c			To invert for 3-parameter, at least 3 of the wps2 and wad2
c			need to be non-zero.
c
c Outputs:
c	  rto     1xndepth array for inverted Rt resistivity (ohm-m)
c	  hdo	  1xndepth array for inverted HD (inch)
c	  erro	  1xndepth array for relative error in inversion
c Cindy 06/11 add: QC  1xndepth array inversion quality indicator
c   	  The relative amplitude of erro can be used as a 
c   	  quality indicator of the inversion.  The better the
c   	  data fit the TTM uniform formation with a borehole model, the 
c	  smaller the value of err_inv. Usually, near bed boundaries 
c   	  or over the zones where invasion, anisotropy, or other effects
c   	  contribute significantly to the
c   	  separation of the resistivity curves of various 
c   	  spacings, the erro will increase significantly
c   	  compared with zones where the data fit the uniform 
c   	  formation with a borehole model
c
c Author:Peter Wu/SPC/Nov 19, 2001
c Revised by Dong/SPC/Aug,2010 for not inverting fixed rt_in or rm_in
c revised by T zhao/Feb 2015, from Keli's ecaliper inversion, change to dielectric inversion
c call Ven_diinv_sub(dh_bhc(i),rm_bhc(i),rt_bhc(i),epsr_bhc(i),err_bhc(i),QC(i),
c	1             id_tool,icontig,i-1,nframe,bs,nch_select,ind_forinv, weight_forinv,psm2i,adm2i,psm4i,adm4i,
c	2             is_OBM,inv_switch,hdi,rmi,rti,epsri)
C ----------------------------------------------------------------------
	subroutine Ven_diinv_sub(hdo2,rmo2,rto2,epsro2,erro2,QC2,
     &     hdo4,rmo4,rto4, epsro4,erro4,QC4,id_tool,ibadtx1,
     &     i_init,ndepth,bsi,nch,index_channel,weight_channel,psm2i,adm2i,psm4i,adm4i,
     &     is_OBM,inv_switch,dh_in,rm_in,rt_in,epsr_in,tool_mode)

c
	parameter (max_m=20) !max. number of measurements
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=26) !5 mbhc spacings,25-point transform
	parameter (ng=100) !no. of initial guess candidates
	parameter (rm_OBM=1000.0) ! OBM mud resistivity
      parameter (max_ms=10) !max. number of measurements for splitting 2M and 400kHz, each has 10 channels
c
     
! Input
      integer tool_mode ! inversion mode just added it for initilization purpose
	integer i_init,ndepth,id_tool,ibadtx1,ibadtx
	real*4 bsi(ndepth),rm_t,hdm
      integer nch
	integer index_channel(nch)
	real*4 weight_channel(nch),w1_in(max_m),w1_in_point(max_m)
	real*4 psm2i(ntr,ndepth),adm2i(ntr,ndepth),psm4i(ntr,ndepth),adm4i(ntr,ndepth)  !air-cal corrected mbhc measurements in db/deg
!	real*4 measi_raw(max_m,ndepth) !air-cal corrected mbhc measurements in db/deg, channels defined by index_channel(1:nch)
	integer is_OBM ! OBM indicator
	integer Nx, inv_switch(4) ! switch for whether inverting dh, rm, rt and ecc, 1 for invert, 0 for not inverting.
	real*4 dh_in,rm_in,rt_in,epsr_in ! input for dh,rm,rt,ecc in cases they are not inverted 

! Output
	real*4 rto2(ndepth),hdo2(ndepth),erro2(ndepth),rmo2(ndepth),epsro2(ndepth),QC2(ndepth)
      real*4 rto4(ndepth),hdo4(ndepth),erro4(ndepth),rmo4(ndepth),epsro4(ndepth),QC4(ndepth)
!      real*4 std_hd2(ndepth), std_rm2(ndepth),std_rt2(ndepth),std_epsr2(ndepth) ! added aug 05,2015
!      real*4 std_hd4(ndepth), std_rm4(ndepth),std_rt4(ndepth),std_epsr4(ndepth)
      real*4 std_hd2, std_rm2,std_rt2,std_epsr2 ! added aug 05,2015
      real*4 std_hd4, std_rm4,std_rt4,std_epsr4
      real    bhc_para(4), bhc_para2(4),bhc_para4(4)   ! borehole model parameter (dh,rm,rt,ecc)

	INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,     ! for SNSL1
	1	NWRITE
	INTEGER IPVT(max_n)
	
	REAL FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
c	REAL FVEC(max_m),FJAC(max_m*max_n),DIAG(max_n),        ! those are old ones for 20 channels
c	1 QTF(max_n),WA1(max_n),WA2(max_n),WA3(max_n),WA4(max_m)
      
      REAL FVEC(max_ms),FJAC(max_ms*max_n),DIAG(max_n),      ! change them for 10 channel measurement inversion
	1 QTF(max_n),WA1(max_n),WA2(max_n),WA3(max_n),WA4(max_ms)
	
	real, allocatable:: X(:),X2(:),X4(:)
	REAL ENORM,R1MACH,ENORM2,ENORM4
	EXTERNAL FCN_bhcmodel2, FCN_bhcmodel2_2M, FCN_bhcmodel2_400k  
      real*4 decadem(7),decaderm(8)
	data (decadem(i),i=1,7) /5.00, 2.5,  1.5, 1, 0.1,  0.01, 0.001 /
      
      data (decaderm(i),i=1,8) /1000.0, 250.000, 100.0,  50.0,
	1	10.0, 1.0, 0.1, 0.01 /
	DATA NWRITE /6/
!	data absent /-999.25/ 
c wps_factor=[20*log10(exp(1))*(pi/180)].^2 
!	data wps_factor /0.0229818/    ! how it is obtained?
	data wps_factor /0.3/

c     adding split input and intermediate result
	real*4 y_in(max_m),ry_in(max_m),w_in(max_m),w_uncty(max_m), y_in2(max_ms),y_in4(max_ms)
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng), el2(ng)
      real*4 rt_ig2(ng),hd_ig2(ng),rm_ig2(ng),epsr_ig2(ng), el22(ng)   ! adding split initialized inversion guess start
      real*4 rt_ig4(ng),hd_ig4(ng),rm_ig4(ng),epsr_ig4(ng), el24(ng)   ! adding split initialized inversion guess start
	real*4 rt_candi2(ng),hd_candi2(ng),error_candi2(ng),rm_candi2(ng),epsr_candi2(ng)
      real*4 rt_candi4(ng),hd_candi4(ng),error_candi4(ng),rm_candi4(ng),epsr_candi4(ng)
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer i1(1),i2(1)
	integer tool_type
	integer No_x, inv_switch_x(4)
c
	real*4 ru(ntf),xu(ntf,max_m),x_nuncty(ntf,max_m)
!	real*4 psm2(ntr,ndepth),adm2(ntr,ndepth),psm4(ntr,ndepth),adm4(ntr,ndepth)  !air-cal corrected mbhc channels in db/deg --- modeled
	real*4 psm2_vec(ntr), adm2_vec(ntr),psm4_vec(ntr), adm4_vec(ntr)
	real*4 sensit
	real*4 rto1,rmo1,hdo1,epsro1, erro1,QC1
	real*4 t1_wm,t2_wm,t3_wm,t4_wm,t5_wm
	real*4 tmp(1)
	real*4 rt_x,rm_x,hd_x,epsr_x  !input model parameters in case it is not inverted
	integer is_OBM_x
	integer kk2, kk4
c      add split weighting for 2M and 400K separately
      real*4 w_in2(10), w_in4(10)
      
	common /data_in_bhc/ y_in,n,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,ibadtx,rm_t,hdm,
     1	tool_type,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in
c
	save
c


      Nx = sum(inv_switch)
      allocate(x2(Nx))
      allocate(x4(Nx))
!      write(150,*) rto,hdo,erro,QC,id_tool,ibadtx1,
!     1	i_init,bsi,rmi,psm2i,adm2i,t1_wm,t2_wm,
!     2     t3_wm,t4_wm,t5_wm,ndepth,psm2,adm2,rmo,Nx,rt_in,rm_in

!     assigning values in the common block data_in_bhc
      tool_type = id_tool
      No_x = Nx
      rt_x = rt_in
      rm_x = rm_in
      hd_x = dh_in
!      ecc_x= 0.0
      epsr_x = epsr_in
      inv_switch_x=inv_switch
      is_OBM_x=is_OBM
      rm_t=rm_in


	 k=nch
	 ldfjac=k
	 m=k  !number of measurements
	 n=Nx !Number of unknowns
c
	 if(m .lt. n) return 

      
	if (i_init .eq. 0) then
	   ibadtx=ibadtx1
	   IOPT = 1
C C
C C     Set FTOL and XTOL to the square root of the machine precision
C C     and GTOL to zero.  Unless high precision solutions are
C C     required, these are the recommended settings.
C C
!	   FTOL = SQRT(R1MACH(4))
!	   XTOL = SQRT(R1MACH(4))
	   FTOL = 1.0e-20
	   XTOL = 1.0e-20

	   GTOL = 0.E0
C C
	   MAXFEV = 400
	   EPSFCN = 0.0
	   MODE = 1
	   FACTOR = 1.E2
	   NPRINT = 0



         do i=1,5
           idx_ps2(i)=i
           idx_ad2(i)=i+ntr
           idx_ps4(i)=i+2*ntr
           idx_ad4(i)=i+3*ntr
         end do

! fetch ps_a,ad_a, ps_b, ad_b values from the table, generate the mixed values and save to xu(ntf,max_m), and save uncertainty too 
! order of xu: rps2,rad2,rps4,rad2 (i.e. following order of idx_ps2,idx_ad2,idx_ps4,idx_ad4)
	   call get_r_tf_table_di(ru,xu,x_nuncty,idx_ps2,idx_ad2,idx_ps4,idx_ad4,  
	1	   id_tool,ibadtx)   
         
         

      end if

c weights builds the weight matrix w1_in(4*ntr)
! index values:  idx_ps2=1:5, idx_ad2(i)=1+ntr:2ntr, idx_ps4=11:15, idx_ad4=16:20

!        write(1001,*) '***************index_channel,weight_channel***********'
!      
!        write(1001,*) index_channel,weight_channel
!        write(1001,*) idx_ps2,idx_ad2,idx_ps4,idx_ad4

       call weights(w1_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,nps,nad,
     &      nch,index_channel,weight_channel,wps_factor)

!        write(1001,*) index_channel,weight_channel
!        write(1001,*) idx_ps2,idx_ad2,idx_ps4,idx_ad4
!
!
!         write(1001,*) 'w1_in', w1_in  
 !!!!!!!!!!!!!!
      i_init=1

	do jj=1,1 ! ndepth=1, point by point processing

!  Initialize model parameters, user input values kept in cases not inverted
        bhc_para(1)= dh_in
        bhc_para(2)= rm_in
        bhc_para(3)= rt_in
        bhc_para(4)= epsr_in
        
 
	 y_in=-999.25  !measurement data in db/deg, in order of ps2,ad2,ps4,ad4
       
       do j=1,nch
	        if (index_channel(j).le.5) y_in(index_channel(j))=psm2i(index_channel(j),jj)
	        if (index_channel(j).gt.5.and.index_channel(j).le.10) y_in(index_channel(j))=adm2i(index_channel(j)-5,jj) 
	        if (index_channel(j).gt.10.and.index_channel(j).le.15) y_in(index_channel(j))=psm4i(index_channel(j)-10,jj) 
	        if (index_channel(j).gt.15) y_in(index_channel(j))=adm4i(index_channel(j)-15,jj) 
       enddo
      
       
       ry_in=-999.25
	  do i=1,20
          if(abs(y_in(i)+999.25).gt.1.0e-5)   call linter_log_v(ry_in(i),i1,i2,    !out put, interpolated apparent resistivity, index before and after point y_in(i)
     &      xu(1,i),ru,y_in(i),0,ntf,1)        !input --- xu,ru: table,   y_in:raw signal at one point 
        end do

        !  add the split part for measurement y_in, 

        do i=1,10
            y_in2(i) = y_in(i)
            y_in4(i) = y_in(i+10)
        end do
        
        
        
      ! get weights assigned

       w_sum=0.0
       w_sum2=0.0
       w_sum4=0.0
       
       do i=1,10
	     if (abs(y_in(i)+999.25).lt.1.0e-5) then
	         w1_in(i)=0.0    !assign weight zero for absent channels
	     else 
	         w_sum=w_sum+w1_in(i)
               w_sum2= w_sum2+w1_in(i)               
	     endif    
       end do
       
       do i=1,10
	     if (abs(y_in(i)+999.25).lt.1.0e-5) then
	         w1_in(i+10)=0.0    !assign weight zero for absent channels
	     else 
	         w_sum=w_sum+w1_in(i+10)
               w_sum4= w_sum4+w1_in(i+10)               
	     endif    
       end do
       

!	 do i=1,20
!	     w_in(i)=w1_in(i)/max(w_sum,1.0)    !normalize weighting
!       enddo
       
       do i=1,10
	     w_in(i)=w1_in(i)/max(w_sum2,1.0)    !normalize weighting
           w_in(i+10)= w1_in(i+10)/max(w_sum4,1.0) 
	 enddo
!         write(1001,*) 'y_in', y_in  

!         write(1001,*) 'w_in', w_in ,w_sum  

 
        

c	
	kk2=0.0      
      kk4=0.0 
      do i=1,10
          if(w_in(i).ne.0.0) kk2=kk2+1
          if(w_in(i+10).ne.0.0) kk4=kk4+1
      enddo
      
      
      if (kk2.lt.Nx)  then
         rmo2=rm_in
         rto2=rt_in
         hdo2=dh_in
         epsro2=epsr_in
         erro2=1.0e8
         QC2=1.0e8 
         deallocate(x2)
           if (kk4.lt.Nx)  then
             rmo4=rm_in
             rto4=rt_in
             hdo4=dh_in
             epsro4=epsr_in
             erro4=1.0e8
             QC4=1.0e8  
             deallocate(x4)
           return
           endif
         
      endif
c     
       hdm=bsi(jj)
      
        if(tool_mode.eq.2) then
            
            call initial_diinv2M_mode2(is_OBM,id_tool,ibadtx,bsi(jj),   !rmi(jj) not used in the subroutine, bsi used for minimum initial borehole size
	1	ry_in,y_in,w_in,m,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig2,hd_ig2,rm_ig2,epsr_ig2,nig)
            
            call initial_diinv400k_mode2(is_OBM,id_tool,ibadtx,bsi(jj),   !rmi(jj) not used in the subroutine, bsi used for minimum initial borehole size
	1	ry_in,y_in,w_in,m,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig4,hd_ig4,rm_ig4,epsr_ig4,nig)
            
      else if (tool_mode.eq.1) then
          
          call initial_diinv2M_mode1(rm_in,is_OBM,id_tool,ibadtx,bsi(jj),   !rmi(jj) not used in the subroutine, bsi used for minimum initial borehole size
	1	ry_in,y_in,w_in,m,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig2,hd_ig2,rm_ig2,epsr_ig2,nig)
            
          call initial_diinv400k_mode1(rm_in,is_OBM,id_tool,ibadtx,bsi(jj),   !rmi(jj) not used in the subroutine, bsi used for minimum initial borehole size
	1	ry_in,y_in,w_in,m,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig4,hd_ig4,rm_ig4,epsr_ig4,nig)
          
         endif
          
c           call initial_gbhc2(id_tool,ibadtx,bsi(jj),   !rmi(jj) not used in the subroutine, bsi used for minimum initial borehole size
c	1	ry_in,y_in,w_in,m,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
c	2	rt_ig,hd_ig,rm_ig,nig)
c        end if	

c        epsr_ig=0.0  ! ?
      
C        if(inv_switch(4).eq.1) then ! if epsr is inverted, why I added it? or should I change to the following condition?

         if(tool_mode .eq. 2) then
         ! keep the first initial guess
             do i=1,3    ! use the first three serached local minimum 
             rt_ig2(i)=rt_ig2(i)
            rm_ig2(i)=rm_ig2(i)
            hd_ig2(i)=hd_ig2(i)
            epsr_ig2(i)= epsr_ig2(i)
        !   if(is_OBM.eq.1) rm_ig2(1)=rm_OBM   ! no OBM for mode 2
           
            rt_ig4(i)=rt_ig4(i)
            rm_ig4(i)=rm_ig4(i)
            hd_ig4(i)=hd_ig4(i)
            epsr_ig4(i)= epsr_ig4(i)
           
            enddo 
         !  if(is_OBM.eq.1) rm_ig4(1)=rm_OBM   ! no OBM for mode 2 
	   ! input model as the fourth initial
	
            rt_ig2(4)=rt_in
            rm_ig2(4)=rm_in
            hd_ig2(4)=dh_in
          !  epsr_ig2(4)= epsr_in
            epsr_ig2(4)=5.0 + 108.5/(rt_ig(4)**0.35)
            rt_ig4(4)=rt_in
            rm_ig4(4)=rm_in
            hd_ig4(4)=dh_in
          !  epsr_ig4(4)= epsr_in
            epsr_ig4(4)=5.0 + 280.0/(rt_ig(4)**0.46)
      ! Construct 64 seems more initial models, is that too many?
            if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then
	          hdmax=20.0  ! for v475 and arc475
            else if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then
                  hdmax=24.0   ! for v675 and arc675
             endif      
                      
            hdm=bsi(jj)
            delta_hd=(hdmax-hdm)/3.0
            delta_epsr=(10000.0-1.0)/3.0            
            delta_rt=(500.0-0.1)/3.0
            delta_rm=(10-0.1)/3.0
             ii=4   ! count of number of optimization
            do i=1,3          ! hd
               do j=1,3       ! epsr
                   do k=1,3   ! Rt
                        do l=1,3  ! Rm
                  ii=ii+1
!                
                  rt_ig2(ii)= 0.1+(k-1)*delta_rt
!                  if(rt_ig(ii).ge.100.0.or.rt_ig(ii).le.0.15) rt_ig(ii)=10.0  ! can modify to have better initial Rt?
!                  rm_ig(ii)=rm_ig(1)
                  rm_ig2(ii)=0.1+(l-1)*delta_rm      !  Rm loop
                 ! if (inv_switch(2).eq.1) rm_ig2(ii)=min(rm_ig2(ii),10.0)  !avoide too high initial Rm
                  hd_ig2(ii)=hdm+(i-1)*delta_hd
                 ! epsr_ig2(ii)=1+(j-1)*delta_epsr
                  ! here is the change on Nov 10, 2015
                  if (rt_ig2(ii) .lt. 10.0) then
                           epsr_ig2(ii)= decaderm(2*j)
                       else if (rt_ig2(ii) .lt. 100.0) then
                           epsr_ig2(ii)= 2*decaderm(2*j)
                       else   
                          epsr_ig2(ii)= 10.0* decaderm(2*j)
                       endif
                  rt_ig4(ii)= 0.1+(k-1)*delta_rt
!                  if(rt_ig(ii).ge.100.0.or.rt_ig(ii).le.0.15) rt_ig(ii)=10.0  ! can modify to have better initial Rt?
!                  rm_ig(ii)=rm_ig(1)
                  rm_ig4(ii)= 0.1+(l-1)*delta_rm 
                  hd_ig4(ii)=hdm+(i-1)*delta_hd
                 ! epsr_ig4(ii)=1+(j-1)*delta_epsr
                  ! here is the change on Nov 10, 2015
                  if (rt_ig4(ii) .lt. 10.0) then
                           epsr_ig4(ii)= decaderm(2*j)
                       else if (rt_ig4(ii) .lt. 100.0) then
                           epsr_ig4(ii)= 2*decaderm(2*j)
                       else   
                          epsr_ig4(ii)= 10.0* decaderm(2*j)
                       endif
                  enddo
                 enddo 
               enddo
            enddo
            nig=ii
               
         endif
         
         
           if(tool_mode .eq. 1) then
         ! keep the first three local minimum initial guess
            
                do i=1,3
            rt_ig2(i)=rt_ig2(i)
            rm_ig2(i)=rm_in
            hd_ig2(i)=dh_in
            epsr_ig2(i)= epsr_ig2(i)
           
            rt_ig4(i)=rt_ig4(i)
            rm_ig4(i)=rm_in
            hd_ig4(i)=dh_in
            epsr_ig4(i)= epsr_ig4(i)
                                               
                  if(is_OBM .eq.1) then
                 rm_ig2(i)=rm_OBM 
                 rm_ig4(i)=rm_OBM 
                  endif
            enddo   
	   ! input model as the fourth  initial
	
            rt_ig2(4)=rt_in
            rm_ig2(4)=rm_in
            hd_ig2(4)=dh_in
            ! epsr_ig2(4)= epsr_in
            epsr_ig2(4)= 5.0 + 108.5/(rt_ig(4)**0.35) 
            rt_ig4(4)=rt_in
            rm_ig4(4)=rm_in
            hd_ig4(4)=dh_in
            !epsr_ig4(4)= epsr_in
            epsr_ig4(4)= 5.0 + 280.0/(rt_ig(4)**0.46) 
            
            if(is_OBM .eq.1) then
                 rm_ig2(4)=rm_OBM 
                 rm_ig4(4)=rm_OBM 
             endif
         
           ! delta_hd=(hdmax-hdm)/5.0
!            delta_epsr=(100000.0-1.0)/5.0   ! get a better way of epsr distribution  
            
            ! epsr intial guesses
!       do j=1,10
!		eps(j)= 100.0*decaderm(j)
!       end do
!           delta_rt=(1500.0-0.1)/5.0
            
             ii=4   ! count of number of optimization
           
               do j=1,8       ! epsr
                   do k=1,7   ! Rt
                  ii=ii+1
!                for 2M
!                  rt_ig2(ii)= 0.1+(k-1)*delta_rt
       !           rt_ig2(ii)=tmp_rt*decadem(k)
                  rm_ig2(ii)=rm_in       !  Rm use input
                   hd_ig2(ii)=dh_in
                 ! epsr_ig2(ii)=1+(j-1)*delta_epsr
                  rt_ig2(ii)=  rt_in*decadem(k)   
                   ! add epsr restriction according to Rt
                       if (rt_ig2(ii) .lt. 10.0) then
                           epsr_ig2(ii)= 0.3*decaderm(j)
                       else if (rt_ig2(ii) .lt. 100.0) then
                           epsr_ig2(ii)= decaderm(j)
                       else   
                          epsr_ig2(ii)= 10.0* decaderm(j)
                       endif
		    !    epsr_ig2(ii)= 100.0*decaderm(j)
                  
                  
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
     
                  ! for 400k
                !  rt_ig4(ii)= 0.1+(k-1)*delta_rt
                   rt_ig4(ii)=  rt_in*decadem(k)     
                  rm_ig4(ii)=rm_in
                  hd_ig4(ii)=dh_in
   !               epsr_ig4(ii)=1+(j-1)*delta_epsr
                   ! add epsr restriction according to Rt
                       if (rt_ig4(ii) .lt. 10.0) then
                           epsr_ig4(ii)= 0.3*decaderm(j)
                       else if (rt_ig4(ii) .lt. 100.0) then
                           epsr_ig4(ii)= decaderm(j)
                       else   
                          epsr_ig4(ii)= 10.0* decaderm(j)
                       endif
                !  epsr_ig2(ii)= 100.0*decaderm(j)
                 
                  if(is_OBM .eq.1) then
                 rm_ig2(ii)=rm_OBM 
                 rm_ig4(ii)=rm_OBM 
                  endif
                  
                  enddo 
               enddo
                                 
             nig=ii
            
      endif
c !!!!!!!!!!!!!!!!!!             
   !    if(is_OBM.eq.1) rm_ig2=rm_OBM
       
        do ii=1,nig
c
      if(tool_mode .eq. 2) then      ! for mode 1, apply limit for hd 
          if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then  ! arc or venus 6
              if (hd_ig2(ii).ge.24.0) then
                  hd_ig2(ii)=24.00
              end if
          else if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then ! arc or venus 5
              if (hd_ig2(ii).ge.18.0) then
                  hd_ig2(ii)=18.00
              end if
          end if
          
          if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then
              if (hd_ig4(ii).ge.24.0) then
                  hd_ig4(ii)=24.00
              end if
          else if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then
              if (hd_ig4(ii).ge.18.0) then
                  hd_ig4(ii)=18.00
              end if
          end if
      endif
            
! Reassign the model parameters initial values (inverted ones only)
           if(inv_switch(1).eq.1)    bhc_para2(1)= hd_ig2(ii)
           if(inv_switch(2).eq.1)    bhc_para2(2)= rm_ig2(ii)
           if(inv_switch(3).eq.1)    bhc_para2(3)= rt_ig2(ii)
           if(inv_switch(4).eq.1)    bhc_para2(4)= epsr_ig2(ii)
           
           if(inv_switch(1).eq.1)    bhc_para4(1)= hd_ig4(ii)
           if(inv_switch(2).eq.1)    bhc_para4(2)= rm_ig4(ii)
           if(inv_switch(3).eq.1)    bhc_para4(3)= rt_ig4(ii)
           if(inv_switch(4).eq.1)    bhc_para4(4)= epsr_ig4(ii)
           
            i=0
           do j=1,4
             if(inv_switch(j).eq.1) then
                i=i+1
                x2(i)=bhc_para2(j)
                x4(i)=bhc_para4(j)
             endif
           enddo

           N = Nx
           
           
           
           ! change the previous 20 channels to 10 here for inversion
           ldfjac=10
           m=10
          
!         write(1001,*) 'initial  bhc_para', bhc_para , nig, ii  
           ! call SNLS1 for 2MHz inversion
           CALL SNLS1(FCN_bhcmodel2_2M,IOPT,M,N,X2,FVEC,FJAC,LDFJAC,FTOL,XTOL,
     1	 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
     2     INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
                      
            FNORM2 = ENORM(M,FVEC)
           
          !  update borehole parameters that were inverted
            i=0
        do j=1,4
           if(inv_switch(j).eq.1) then
              i=i+1
              bhc_para2(j)=x2(i)
           endif
        enddo
        
        
!          
!   recomputing fitting error for 2MHz inverted 
       hdi2=bhc_para2(1)
       rmi2=bhc_para2(2)
       rti2=bhc_para2(3)
       epsri2=bhc_para2(4)
       
       call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rmi2,rti2,hdi2,epsri2,6)

       call el2norm_2M(y_in,w_in,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,errpc)

! recording inversion results and fitting error
        hd_candi2(ii)=bhc_para2(1)
        rm_candi2(ii)=bhc_para2(2)
        rt_candi2(ii)=bhc_para2(3)
        epsr_candi2(ii)=bhc_para2(4)
        error_candi2(ii)=errpc
        
        
        ! initialize again to make sure SNLS1 is working
         IOPT = 1
	   FTOL = 1.0e-20
	   XTOL = 1.0e-20

	   GTOL = 0.E0
C C
	   MAXFEV = 400
	   EPSFCN = 0.0
	   MODE = 1
	   FACTOR = 1.E2
	   NPRINT = 0
           ! call SNLS1 for 400kHz inversion
           CALL SNLS1(FCN_bhcmodel2_400k,IOPT,M,N,X4,FVEC,FJAC,LDFJAC,FTOL,XTOL,
     1	GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
     2  INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)


	    FNORM4 = ENORM(M,FVEC)
          
           !  update borehole parameters that were inverted
         i=0
        do j=1,4
           if(inv_switch(j).eq.1) then
              i=i+1
              bhc_para4(j)=x4(i)
           endif
        enddo


!          
!   recomputing fitting error for 2MHz inverted 
       hdi4=bhc_para4(1)
       rmi4=bhc_para4(2)
       rti4=bhc_para4(3)
       epsri4=bhc_para4(4)
       
       call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rmi4,rti4,hdi4,epsri4,6)

       call el2norm_400k(y_in,w_in,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,errpc)

! recording inversion results and fitting error
        hd_candi4(ii)=bhc_para4(1)
        rm_candi4(ii)=bhc_para4(2)
        rt_candi4(ii)=bhc_para4(3)
        epsr_candi4(ii)=bhc_para4(4)
        error_candi4(ii)=errpc
!      
  
        
! 
        end do  !end of the ii loop for nig=18 initial guess

           
        ! Pick the one with smallest error from 18 initial guesses?
         error_min2=1.e20
         error_min4=1.e20
      	 do ii=1,nig     !nig=85?
	       if (error_candi2(ii).le.error_min2) then
	          error_min2=error_candi2(ii)
	          erro2(jj)=error_min2
		      hdo2(jj)=hd_candi2(ii)
		      rto2(jj)=rt_candi2(ii)
                rmo2(jj)=rm_candi2(ii)
                epsro2(jj)=epsr_candi2(ii)
             end if 
             
              if (error_candi4(ii).le.error_min4) then
                 error_min4=error_candi4(ii) 
                 erro4(jj)=error_min4
		       hdo4(jj)=hd_candi4(ii)
		       rto4(jj)=rt_candi4(ii)
                 rmo4(jj)=rm_candi4(ii)
                 epsro4(jj)=epsr_candi4(ii)
	        end if
           end do
           
c  set the i_init to 1 to prevent initialization again in the loop through depth 
c   get QC from sensitivity and fitting error
      rmo1=rmo2(jj)
      rto1=rto2(jj)
      hdo1=hdo2(jj)
      epsro1=epsro2(jj)
      erro1=erro2(jj)
      
!      call calculate_std_hessian2M(rmo1,rto1,hdo1,epsro1,std_hd2,std_rm2,std_rt2,std_epsr2)
      
      call calculateQC_2M(rmo1,rto1,hdo1,epsro1,erro1, QC1)   ! haven't looked into this QC subrountine yet
      QC2(jj)=QC1 
      
      rmo1=rmo4(jj)
      rto1=rto4(jj)
      hdo1=hdo4(jj)
      epsro1=epsro4(jj)
      erro1=erro4(jj)
      
!      call calculate_std_hessian400k(rmo1,rto1,hdo1,epsro1,std_hd4,std_rm4,std_rt4,std_epsr4)
      
      call calculateQC_400k(rmo1,rto1,hdo1,epsro1,erro1, QC1)  
      
      QC4(jj)=QC1 
!        QC(jj)=1.0   ! temp set, remove after calculateQC is done
        
        
       end do
	deallocate(x2)
      deallocate(x4)
	return
c
      END
c
c
      subroutine el2norm_2M(y_in,w_in,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,cpsm2_m,cadm2_m,cpsm4_m,cadm4_m,el2_o)
c
	parameter (ntr=5,max_m=20,max_ms=10)
	real*4 y_in(max_m),w_in(max_m)
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(max_m)
	real*4 cpsm2_m(ntr),cadm2_m(ntr),cpsm4_m(ntr),cadm4_m(ntr)
	real*4 el2_o,sumsqr,y_mod(max_m)
	real*4 sum_sig
c
c
        do i=1,ntr
          if(idx_ps2(i) .ne. 0) y_mod(idx_ps2(i))=cpsm2_m(i)
        end do
c
        do i=1,ntr
          if(idx_ad2(i) .ne. 0) y_mod(idx_ad2(i))=cadm2_m(i)
        end do

        do i=1,ntr
          if(idx_ps4(i) .ne. 0) y_mod(idx_ps4(i))=cpsm4_m(i)
        end do
c
        do i=1,ntr
          if(idx_ad4(i) .ne. 0) y_mod(idx_ad4(i))=cadm4_m(i)
        end do
c
      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)
        

	sumsqr=0.0
	sum_sig=0.0
c
	do i=1,max_ms
	  if(index_all(i).ne.0) then
	     sumsqr=sumsqr + w_in(i)*w_in(i)*(y_in(i)-y_mod(i))*(y_in(i)-y_mod(i))
	     sum_sig=sum_sig+w_in(i)*w_in(i)*(y_in(i))*(y_in(i))
	  endif
	end do
c
	el2_o=sqrt(sumsqr/sum_sig)
c
	return
      end
      
       subroutine el2norm_400k(y_in,w_in,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,cpsm2_m,cadm2_m,cpsm4_m,cadm4_m,el2_o)
c
	parameter (ntr=5,max_m=20,max_ms=10)
	real*4 y_in(max_m),w_in(max_m)
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(max_m)
	real*4 cpsm2_m(ntr),cadm2_m(ntr),cpsm4_m(ntr),cadm4_m(ntr)
	real*4 el2_o,sumsqr,y_mod(max_m)
	real*4 sum_sig
c
c
        do i=1,ntr
          if(idx_ps2(i) .ne. 0) y_mod(idx_ps2(i))=cpsm2_m(i)
        end do
c
        do i=1,ntr
          if(idx_ad2(i) .ne. 0) y_mod(idx_ad2(i))=cadm2_m(i)
        end do

        do i=1,ntr
          if(idx_ps4(i) .ne. 0) y_mod(idx_ps4(i))=cpsm4_m(i)
        end do
c
        do i=1,ntr
          if(idx_ad4(i) .ne. 0) y_mod(idx_ad4(i))=cadm4_m(i)
        end do
c
      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)
        

	sumsqr=0.0
	sum_sig=0.0
c
	do i=1,max_ms
	  if(index_all(i+10).ne.0) then
	     sumsqr=sumsqr + w_in(i+10)*w_in(i+10)*(y_in(i+10)-y_mod(i+10))*(y_in(i+10)-y_mod(i+10))
	     sum_sig=sum_sig+w_in(i+10)*w_in(i+10)*(y_in(i+10))*(y_in(i+10))
	  endif
	end do
c
	el2_o=sqrt(sumsqr/sum_sig)
c
	return
	end

      subroutine initial_diinv(id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,epsr_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m), rm(10), eps(10)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng)
c	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
c	real*4 el2(nrmax,nhdmax)
      real*4 rt(nrmax),hd(14),el2_min(ng)
	real*4 el2(nrmax,14,10,10)
	real*4 decadem(nrmax),decaderm(10)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 3.000,  2.044,  1.392,  0.949,
	1	0.646,  0.440,  0.300, 0.15 /
      
      data (decaderm(i),i=1,10) /1000.0, 500.00 , 300.000, 200.044,  100.392,  50.949,
	1	10.646, 1.440, 0.1, 0.01 /
	

	data rmax / 5000.0 /
	data rmin / 0.1 /
c
	save all
      
      epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12    continue

       ! Rt initial guesses 9
      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	do i=1,nr
		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(j)
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
      end do
c
!    Hd intial value guesses 14
 !     delta_hd=(20.0-hd_min)/(nhdmax-1)
      
  !      delta_hd=(20.0-hd_min)/14.0  ! for v475
      delta_hd=(24.0-hd_min)/14.0 
!      delta_hd=(42.0-real(id_tool))/(nhdmax-1)
!      tmp=max((hd_min-4.0),real(id_tool))
!      delta_hd=(42.0-tmp)/(nhdmax-1)
!	do i=1,nhdmax
        
        do i=1,14
	!	hd(i)=hd_min+(i-4)*delta_hd   ! for v475
           hd(i)=hd_min+(i-1)*delta_hd  ! v675 modify
      end do
c    
      ! Rm inital guesses  10
     
       do i=1,10
		rm(i)=0.015+ (i-1)*(10.0-0.015)/10.0
			if (rm(i) .gt. 10.0) rm(i)=10.0
			if (rm(i) .lt. 0.015) rm(i)=0.015
      end do
      
      ! epsr intial guesses
       do j=1,10
		eps(j)= 100.0*decaderm(j)
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
      end do


	 do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,10  !eps loop
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm(l),rt(j),hd(k),eps(m),ibadtx)
		 call el2norm(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,k,l,m))  
	            end do
	         end do
            end do
       end do

c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,10  !eps loop
                       if(el2(j,k,l,m) .le. el2_min(1)) then
                          el2_min(1)=el2(j,k,l,m)
                          rt_ig(1)=rt(j)
                          hd_ig(1)=hd(k)
                          rm_ig(1)=rm(l)
                          epsr_ig(1)=eps(m)
                          j1=j
                          k1=k
                          l1=l
                          m1=m
                      end if
                  end do
	         end do
            end do
       end do
c
      el2_min(2)=1.0e20
       do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,10  !eps loop
              if ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) 
	1          go to 2500
                       if(el2(j,k,l,m) .le. el2_min(2)) then
                         el2_min(2)=el2(j,k,l,m)
                          rt_ig(2)=rt(j)
                          hd_ig(2)=hd(k)
                          rm_ig(2)=rm(l)
                          epsr_ig(2)=eps(m)
                          j2=j
                          k2=k
                          l2=l
                          m2=m
                         end if
2500	      continue
             end do
	     end do
        end do
       end do
c
      el2_min(3)=1.0e20
      do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,10  !eps loop
            if ( ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) .or.
	1	((j.eq.j2).and.(k.eq.k2).and.(l.eq.l2).and.(m.eq.m2)))
	1        go to 3500
              if(el2(j,k,l,m) .le. el2_min(3)) then
                el2_min(3)=el2(j,k,l,m)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd(k)
                rm_ig(3)=rm(l)
                epsr_ig(3)=eps(m)
              end if
3500	      continue
             end do
	      end do
           end do
       end do

c

	return
      end   !end initial_diinv
      
       subroutine initial_diinv2M_mode2(is_OBM,id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,epsr_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m), rm(10), eps(16)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng)
c	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
c	real*4 el2(nrmax,nhdmax)
      real*4 rt(nrmax),hd(14),el2_min(ng)
	real*4 el2(nrmax,14,10,16)
	real*4 decadem(nrmax),decaderm(16)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 2.5,  2.0,  1.5,  1.0,
	1	0.75,  0.5,  0.25, 0.1 /
      
      data (decaderm(i),i=1,16) /1000.0, 750.0, 500.00 , 250.000, 100.0, 75.0,
	1	50.0, 25.0,10.0, 5.0, 2.50, 1.0, 0.5, 0.25, 0.1, 0.01 /
	

	data rmax / 1500.0 /
	data rmin / 0.1 /
c
	save all
      
 !     epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12    continue

       ! Rt initial guesses 9
      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	do i=1,nr
!		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(i)
              !rt(i)= 0.1+ (i-1)*(1500.0-0.1)/nr
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
              end do
c
!    Hd intial value guesses 14
 !     delta_hd=(20.0-hd_min)/(nhdmax-1)
          
              if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then
                  hd_min=5.25  ! for arc/v475
                  delta_hd=(24.0-hd_min)/14.0  ! for v475
             else if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then
                 hd_min=7.0   ! for arc/v675
                 delta_hd=(30.0-hd_min)/14.0 
             endif 
               do i=1,14
	!	hd(i)=hd_min+(i-4)*delta_hd   ! for v475
                hd(i)=hd_min+(i-1)*delta_hd  ! v675 modify
               end do       
         
      
!      delta_hd=(42.0-real(id_tool))/(nhdmax-1)
!      tmp=max((hd_min-4.0),real(id_tool))
!      delta_hd=(42.0-tmp)/(nhdmax-1)
!	do i=1,nhdmax
        
             
c    
      ! Rm inital guesses  10
                
       do i=1,10
 !          if (is_OBM.eq.1) then            ! mode 2 cannot have OBM
 !              rm(i)= 1000.0
 !              else
		rm(i)=0.015+ (i-1)*(10.0-0.015)/10.0
			if (rm(i) .gt. 10.0) rm(i)=10.0
			if (rm(i) .lt. 0.015) rm(i)=0.015
              
!          end if
      end do
      
      ! epsr intial guesses
       do j=1,16
		eps(j)= 100.0*decaderm(j)
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
      end do


	 do j=1,nr				!rt loop  9 
c        do k=1,nhdmax		!HD loop 14
           do k=1,14
               do l=1,10   ! rm loop  10
                   do m=1,16  !eps loop  10
                       
                         ! add epsr restriction according to Rt
                       if (rt(j) .lt. 10.0) then
                           eps(m)= 0.3*decaderm(m)
                       else if (rt(j) .lt. 100.0) then
                           eps(m)= decaderm(m)
                       else   
                          eps(m)= 10.0* decaderm(m)
                       endif
                       
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm(l),rt(j),hd(k),eps(m),ibadtx)
		 call el2norm_2M(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,k,l,m))  
	            end do
	         end do
            end do
       end do

c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
                       if(el2(j,k,l,m) .le. el2_min(1)) then
                          el2_min(1)=el2(j,k,l,m)
                          rt_ig(1)=rt(j)
                          hd_ig(1)=hd(k)
                          rm_ig(1)=rm(l)
                          epsr_ig(1)=eps(m)
                          j1=j
                          k1=k
                          l1=l
                          m1=m
                      end if
                  end do
	         end do
            end do
       end do
c
      el2_min(2)=1.0e20
       do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
              if ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) 
	1          go to 2500
                       if(el2(j,k,l,m) .le. el2_min(2)) then
                         el2_min(2)=el2(j,k,l,m)
                          rt_ig(2)=rt(j)
                          hd_ig(2)=hd(k)
                          rm_ig(2)=rm(l)
                          epsr_ig(2)=eps(m)
                          j2=j
                          k2=k
                          l2=l
                          m2=m
                         end if
2500	      continue
             end do
	     end do
        end do
       end do
c
      el2_min(3)=1.0e20
      do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
            if ( ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) .or.
	1	((j.eq.j2).and.(k.eq.k2).and.(l.eq.l2).and.(m.eq.m2)))
	1        go to 3500
              if(el2(j,k,l,m) .le. el2_min(3)) then
                el2_min(3)=el2(j,k,l,m)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd(k)
                rm_ig(3)=rm(l)
                epsr_ig(3)=eps(m)
              end if
3500	      continue
             end do
	      end do
           end do
          end do

         return
      end               ! initial_diinv2M_mode2
      

      
      subroutine initial_diinv400k_mode2(is_OBM,id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,epsr_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m), rm(10), eps(16)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng)
c	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
c	real*4 el2(nrmax,nhdmax)
      real*4 rt(nrmax),hd(14),el2_min(ng)
	real*4 el2(nrmax,14,10,16)
	real*4 decadem(nrmax),decaderm(16)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 2.5,  2.0,  1.5,  1.0,
	1	0.75,  0.5,  0.25, 0.1 /
      
      data (decaderm(i),i=1,16) /1000.0, 750.0, 500.00 , 250.000, 100.0, 75.0,
	1	50.0, 25.0,10.0, 5.0, 2.50, 1.0, 0.5, 0.25, 0.1, 0.01 /
	

	data rmax / 1500.0 /
	data rmin / 0.1 /
c
	save all
      
 !     epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12    continue

       ! Rt initial guesses 9
      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	do i=1,nr
!		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(i)
              ! rt(i)= 0.1+ (i-1)*(1500.0-0.1)/nr
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
      end do
c
!    Hd intial value guesses 14
 !     delta_hd=(20.0-hd_min)/(nhdmax-1)
      if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then
          hd_min=5.25  ! for v475
      else if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then
          hd_min=7.0 
      endif 
      
      if ((id_tool .eq. 4) .or. (id_tool .eq. 5)) then
        delta_hd=(24.0-hd_min)/14.0  ! for v475
      else if ((id_tool .eq. 7) .or. (id_tool .eq. 6)) then
          delta_hd=(30.0-hd_min)/14.0 
      endif 
      
!      delta_hd=(42.0-real(id_tool))/(nhdmax-1)
!      tmp=max((hd_min-4.0),real(id_tool))
!      delta_hd=(42.0-tmp)/(nhdmax-1)
!	do i=1,nhdmax
        
        do i=1,14
	!	hd(i)=hd_min+(i-4)*delta_hd   ! for v475
           hd(i)=hd_min+(i-1)*delta_hd  ! v675 modify
        end do
c    
      ! Rm inital guesses  10
       do i=1,10
!           if (is_OBM.eq.1) then 
!               rm(i)= 1000.0
!          else
		rm(i)=0.015+ (i-1)*(10.0-0.015)/10.0
			if (rm(i) .gt. 10.0) rm(i)=10.0
			if (rm(i) .lt. 0.015) rm(i)=0.015
!          end if 
      end do
      
      ! epsr intial guesses
       do j=1,16
		eps(j)= 100.0*decaderm(j)
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
      end do


	 do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
                         ! add epsr restriction according to Rt
                       if (rt(j) .lt. 10.0) then
                           eps(m)= 0.3*decaderm(m)
                       else if (rt(j) .lt. 100.0) then
                           eps(m)= decaderm(m)
                       else   
                          eps(m)= 10.0* decaderm(m)
                       endif
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm(l),rt(j),hd(k),eps(m),ibadtx)
		 call el2norm_400k(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,k,l,m))  
	            end do
	         end do
            end do
       end do

c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
                       if(el2(j,k,l,m) .le. el2_min(1)) then
                          el2_min(1)=el2(j,k,l,m)
                          rt_ig(1)=rt(j)
                          hd_ig(1)=hd(k)
                          rm_ig(1)=rm(l)
                          epsr_ig(1)=eps(m)
                          j1=j
                          k1=k
                          l1=l
                          m1=m
                      end if
                  end do
	         end do
            end do
       end do
c
      el2_min(2)=1.0e20
       do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
              if ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) 
	1          go to 2500
                       if(el2(j,k,l,m) .le. el2_min(2)) then
                         el2_min(2)=el2(j,k,l,m)
                          rt_ig(2)=rt(j)
                          hd_ig(2)=hd(k)
                          rm_ig(2)=rm(l)
                          epsr_ig(2)=eps(m)
                          j2=j
                          k2=k
                          l2=l
                          m2=m
                         end if
2500	      continue
             end do
	     end do
        end do
       end do
c
      el2_min(3)=1.0e20
      do j=1,nr				!rt loop
c        do k=1,nhdmax		!HD loop
           do k=1,14
               do l=1,10   ! rm loop
                   do m=1,16  !eps loop
            if ( ((j.eq.j1).and.(k.eq.k1).and.(l.eq.l1).and.(m.eq.m1)) .or.
	1	((j.eq.j2).and.(k.eq.k2).and.(l.eq.l2).and.(m.eq.m2)))
	1        go to 3500
              if(el2(j,k,l,m) .le. el2_min(3)) then
                el2_min(3)=el2(j,k,l,m)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd(k)
                rm_ig(3)=rm(l)
                epsr_ig(3)=eps(m)
              end if
3500	      continue
             end do
	      end do
           end do
       end do

c

	return
      end   !end initial_diinv400k_mode2
      
        subroutine initial_diinv2M_mode1(rm_fix,is_OBM,id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,epsr_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m), rm, eps(16)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min,rm_fix
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng)
c	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
c	real*4 el2(nrmax,nhdmax)
      real*4 rt(nrmax),hd(14),el2_min(ng)
	real*4 el2(nrmax,16)
	real*4 decadem(nrmax),decaderm(16)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 2.5,  2.0,  1.5,  1.0,
	1	0.75,  0.5,  0.25, 0.1 /
      
      data (decaderm(i),i=1,16) /1000.0, 750.0, 500.00 , 250.000, 100.0, 75.0,
	1	50.0, 25.0,10.0, 5.0, 2.50, 1.0, 0.5, 0.25, 0.1, 0.01 /
	

	data rmax / 1500.0 /
	data rmin / 0.1 /
c
	save all
      
 !     epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12    continue

       ! Rt initial guesses 9
      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	      do i=1,nr
!		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(i)    ! rt is not far from tmp_rt
             ! rt(i)= 0.1+ (i-1)*(1500.0-0.1)/nr
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
              end do
c
  
      ! Rm = rm_fix
       
           if (is_OBM.eq.1) then 
               rm= 1000.0
               else
		rm= rm_fix  ! it's actually rm_input   from input, for mode 1, rm and hd are fixed
          end if
      
      
      ! epsr intial guesses, can consider to increase the number to 15?
       do j=1,16
		eps(j)= 100.0*decaderm(j)
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
      end do


	 do j=1,nr				!rt loop  9 
              do m=1,16  !eps loop  10
                    ! add epsr restriction according to Rt
                       if (rt(j) .lt. 10.0) then
                           eps(m)= 0.3*decaderm(m)
                       else if (rt(j) .lt. 100.0) then
                           eps(m)= decaderm(m)
                       else   
                          eps(m)= 10.0* decaderm(m)
                       endif
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm,rt(j),hd_min,eps(m),ibadtx)
		 call el2norm_2M(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,m))  
	       end do
       end do

c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
             do m=1,16  !eps loop
                       if(el2(j,m) .le. el2_min(1)) then
                          el2_min(1)=el2(j,m)
                          rt_ig(1)=rt(j)
                          hd_ig(1)=hd_min
                          rm_ig(1)=rm
                          epsr_ig(1)=eps(m)
                          j1=j
                          m1=m
                      end if
           end do
       end do
c
      el2_min(2)=1.0e20
       do j=1,nr				!rt loop
                do m=1,16  !eps loop
              if ((j.eq.j1).and.(m.eq.m1)) 
	1          go to 2500
                       if(el2(j,m) .le. el2_min(2)) then
                         el2_min(2)=el2(j,m)
                          rt_ig(2)=rt(j)
                          hd_ig(2)=hd_min
                          rm_ig(2)=rm
                          epsr_ig(2)=eps(m)
                          j2=j
                          m2=m
                         end if
2500	      continue
            
        end do
       end do
c
      el2_min(3)=1.0e20
      do j=1,nr				!rt loop
             do m=1,16  !eps loop
            if ( ((j.eq.j1).and.(m.eq.m1)) .or.
	1	((j.eq.j2).and.(m.eq.m2)))
	1        go to 3500
              if(el2(j,m) .le. el2_min(3)) then
                el2_min(3)=el2(j,m)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd_min
                rm_ig(3)=rm
                epsr_ig(3)=eps(m)
              end if
3500	      continue
             
           end do
          end do

         return
        end               ! initial_diinv2M_mode1
      
      subroutine initial_diinv400k_mode1(rm_fix,is_OBM,id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,epsr_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m), rm, eps(16)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min,rm_fix
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng),epsr_ig(ng)
c	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
c	real*4 el2(nrmax,nhdmax)
      real*4 rt(nrmax),hd(14),el2_min(ng)
	real*4 el2(nrmax,16)
	real*4 decadem(nrmax),decaderm(16)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 2.5,  2.0,  1.5,  1.0,
	1	0.75,  0.5,  0.25, 0.1 /
      
      data (decaderm(i),i=1,16) /1000.0, 750.0, 500.00 , 250.000, 100.0, 75.0,
	1	50.0, 25.0,10.0, 5.0, 2.50, 1.0, 0.5, 0.25, 0.1, 0.01 /
	

	data rmax / 1500.0 /
	data rmin / 0.1 /
c
	save all
      
 !     epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12    continue

       ! Rt initial guesses 9
      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	do i=1,nr
!		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(i)
             ! rt(i)= 0.1+ (i-1)*(1500.0-0.1)/nr
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
      end do
c
!  
c    ! Rm = rm_fix
       
           if (is_OBM.eq.1) then 
               rm= 1000.0
               else
		     rm= rm_fix  ! it's actually rm_input   from input, for mode 1, rm and hd are fixed
               end if
               
      ! epsr intial guesses
       do j=1,16
		eps(j)= 100.0*decaderm(j)
!			if (eps(j) .gt. 100000.0) eps(j)=10.0
!			if (eps(j) .lt. 1) eps(j)=1
      end do


	 do j=1,nr				!rt loop
                do m=1,16  !eps loop
                      ! add epsr restriction according to Rt
                       if (rt(j) .lt. 10.0) then
                           eps(m)= 0.3*decaderm(m)
                       else if (rt(j) .lt. 100.0) then
                           eps(m)= decaderm(m)
                       else   
                          eps(m)= 10.0* decaderm(m)
                       endif
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm,rt(j),hd_min,eps(m),ibadtx)
		 call el2norm_400k(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,m))  
	            
            end do
       end do

c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
               do m=1,16  !eps loop
                       if(el2(j,m) .le. el2_min(1)) then
                          el2_min(1)=el2(j,m)
                          rt_ig(1)=rt(j)
                          hd_ig(1)=hd_min
                          rm_ig(1)=rm
                          epsr_ig(1)=eps(m)
                          j1=j
                          m1=m
                      end if
                  
            end do
       end do
c
      el2_min(2)=1.0e20
       do j=1,nr				!rt loop
             do m=1,16  !eps loop
              if ((j.eq.j1).and.(m.eq.m1)) 
	1          go to 2500
                       if(el2(j,m) .le. el2_min(2)) then
                         el2_min(2)=el2(j,m)
                          rt_ig(2)=rt(j)
                          hd_ig(2)=hd_min
                          rm_ig(2)=rm
                          epsr_ig(2)=eps(m)
                          j2=j
                          m2=m
                         end if
2500	      continue
             
        end do
       end do
c
      el2_min(3)=1.0e20
      do j=1,nr				!rt loop
            do m=1,16  !eps loop
            if ( ((j.eq.j1).and.(m.eq.m1)) .or.
	1	((j.eq.j2).and.(m.eq.m2)))
	1        go to 3500
              if(el2(j,m) .le. el2_min(3)) then
                el2_min(3)=el2(j,m)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd_min
                rm_ig(3)=rm
                epsr_ig(3)=eps(m)
              end if
3500	      continue
          end do
       end do

	return
      end   !end initial_diinv400k_mode1
      
c
      subroutine get_r_tf_table_di(ru,xu,x_nuncty,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	1	id_tool,ibadtx)
c
c this subroutine will fetch the raw transform table and do the mixing
c according to the ibadtx and put the mixed output into xu, used for converting raw signal (db/deg) into apparent resistivity
! xu(ntf,max_m) --- compensated channels (in db.deg) for each Rt
c
c	include 'arc_ttmbhctab_def.inc'
      include 'venus_ttmdivtab_def.inc'
c
	parameter (max_m=20) !max. number of measurements
	parameter (ntf=26) !5 mbhc spacings,25-point transform (same as nrt in the bhc table), while I have 26 point, rt=1500 I added
	real*4 ru(ntf),xu(ntf,max_m),x_nuncty(ntf,max_m)
	real*4 arc_nuncty(ntf,ntr+1),tmp(ntr),tmpm(ntr)
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),id_tool,ibadtx
c 
c ARC Normalized uncertainty Table
c
c 
c j =  1, TR spacing =     10 
c 
	data (arc_nuncty(i,1),i=1,ntf)  /  
	1   3.1799124 ,  2.2706003 ,  1.8904639 ,  1.5525811 ,  1.3090094,
	2   1.2112535 ,  1.1406040 ,  1.0875742 ,  1.0619618 ,  1.0373388,
	3   1.0189127 ,  1.0116578 ,  1.0066631 ,  1.0032030 ,  1.0016893,
	4   1.0003995 ,  0.9996229 ,  0.9994050 ,  0.9993166 ,  0.9993187,
	5   0.9993573 ,  0.9994377 ,  0.9995557 ,  0.9996340 ,  0.9999999, 
     6   1.0000000/ 
c 
c j =  2, TR spacing =     16 
c 
	data (arc_nuncty(i,2),i=1,ntf)  /  
	1   9.1878420 ,  5.1146661 ,  3.6801229 ,  2.5526448 ,  1.8308696,
	2   1.5624205 ,  1.3750438 ,  1.2369601 ,  1.1704479 ,  1.1059422,
	3   1.0565174 ,  1.0364189 ,  1.0221499 ,  1.0118803 ,  1.0071991,
	4   1.0030190 ,  1.0002772 ,  0.9993840 ,  0.9989020 ,  0.9987058,
	5   0.9987020 ,  0.9988029 ,  0.9990166 ,  0.9991765 ,  0.9999999, 
     6   1.0000000/ 
c 
c j =  3, TR spacing =     22 
c 
	data (arc_nuncty(i,3),i=1,ntf)  /  
	1  27.1863471 , 11.8944749 ,  7.4228914 ,  4.3597861 ,  2.6592791,
	2   2.0887151 ,  1.7123474 ,  1.4465203 ,  1.3216284 ,  1.2018978,
	3   1.1102908 ,  1.0726683 ,  1.0455733 ,  1.0256691 ,  1.0163838,
	4   1.0078721 ,  1.0020389 ,  1.0000150 ,  0.9988231 ,  0.9982113,
	5   0.9980730 ,  0.9981214 ,  0.9983970 ,  0.9986364 ,  0.9999999,
     6   1.0000000/ 
c 
c j =  4, TR spacing =     28 
c 
	data (arc_nuncty(i,4),i=1,ntf)  /  
	1  84.5199589 , 29.0974070 , 15.7294369 ,  7.7999587 ,  4.0263198,
	2   2.8999178 ,  2.2052261 ,  1.7406140 ,  1.5301236 ,  1.3327499,
	3   1.1839245 ,  1.1229038 ,  1.0786841 ,  1.0457914 ,  1.0301989,
	4   1.0156282 ,  1.0053184 ,  1.0015847 ,  0.9992690 ,  0.9979462,
	5   0.9975425 ,  0.9974284 ,  0.9977060 ,  0.9980143 ,  0.9999999,
     6   1.0000000/ 
c 
c i =  5, TR spacing =     34 
c 
	data (arc_nuncty(i,5),i=1,ntf)  / 
	1 267.5452000 , 73.2613691 , 34.3450288 , 14.3715696 ,  6.2660245,
	2   4.1312639 ,  2.9079405 ,  2.1388506 ,  1.8053206 ,  1.5015737,
	3   1.2779715 ,  1.1872480 ,  1.1215252 ,  1.0723499 ,  1.0488003,
	4   1.0264894 ,  1.0103219 ,  1.0042804 ,  1.0003979 ,  0.9980339,
	5   0.9972112 ,  0.9967971 ,  0.9969905 ,  0.9973438 ,  0.9999999, 
     6   1.0000000/ 
c 
c i =  6, TR spacing =     40 
c 
	data (arc_nuncty(i,6),i=1,ntf)  / 
	1 833.4802259 ,188.0179040 , 76.7684933 , 27.1230174 ,  9.9748024,
	2   6.0123833 ,  3.9111080 ,  2.6752695 ,  2.1648209 ,  1.7155213,
	3   1.3947910 ,  1.2668697 ,  1.1746986 ,  1.1056934 ,  1.0724651,
	4   1.0406864 ,  1.0172413 ,  1.0082668 ,  1.0023451 ,  0.9985772,
	5   0.9971612 ,  0.9962853 ,  0.9962852 ,  0.9966490 ,  0.9999999, 
     6   1.0000000/ 
c 
	do i=1,ntf
		ru(i)=rt(i)
		
		!2MHz phaseshift
		do j=1,ntr
			tmp(j)=ps_a(i,1,nrm,j,1)  ! **** iHD=1, irm=nrm, iecc=1)
!			write(170,*) j, tmp(j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do j=1,ntr
		 if(idx_ps2(j) .ne. 0) then
			xu(i,idx_ps2(j))=tmpm(j)
			if(id_tool .le. 5) then
				x_nuncty(i,idx_ps2(j))=arc_nuncty(i,j)
			else
				x_nuncty(i,idx_ps2(j))=arc_nuncty(i,j+1)
			end if
		 end if
		end do
		
		!2MHz attenuation
		do j=1,ntr
			tmp(j)=ad_a(i,1,nrm,j,1)    !******** 2MHz attenuation, need to add 400kHz
!			write(170,*) j, tmp(j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do j=1,ntr
		 if(idx_ad2(j) .ne. 0) then
			xu(i,idx_ad2(j))=tmpm(j)
			if(id_tool .le. 5) then
				x_nuncty(i,idx_ad2(j))=arc_nuncty(i,j)
			else
				x_nuncty(i,idx_ad2(j))=arc_nuncty(i,j+1)
			end if
		 end if
		end do
		

		!400kHz phaseshift
		do j=1,ntr
			tmp(j)=ps_b(i,1,nrm,j,1)  ! **** iHD=1, irm=nrm, iecc=1)
!			write(170,*) j, tmp(j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do j=1,ntr
		 if(idx_ps4(j) .ne. 0) then
			xu(i,idx_ps4(j))=tmpm(j)
		 end if
		end do
		
		!400kHz attenuation
		do j=1,ntr
			tmp(j)=ad_b(i,1,nrm,j,1)    !******** 2MHz attenuation, need to add 400kHz
!			write(170,*) j, tmp(j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do j=1,ntr
		 if(idx_ad4(j) .ne. 0) then
			xu(i,idx_ad4(j))=tmpm(j)
		 end if
		end do

		
		
	end do
	return
       end


      SUBROUTINE FCN_bhcmodel2(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
	include 'venus_ttmdivtab_def.inc'
c
	parameter (max_m=20)
!	parameter (rt_min=0.05,rt_max=5000.0)
	parameter (rt_min=0.1,rt_max=1500.0)
	parameter (rm_min=0.015,rm_max=10.0)
      real*4 w_in(max_m),y_in(max_m),y_m(max_m)
!	real*4 dpsm2(ntr),dadm2(ntr),dpsm4(ntr),dadm4(ntr)
!        real*4 psm2(ntr),adm2(ntr),rm0(18),error(18)
      real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr)
      real*4 hd_min, hd_max,epsr_min,epsr_max
      real*4 rt_x,rm_x,hd_x,epsr_x
      integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(max_m)
	integer n_unknow,m_data,k,j
	integer id_tool,No_x
	integer inv_switch_x(4)
	real*4 bhc_para(4) ! borehole model parameter (hd,rm,rt,ecc)
	real*4 hd_t,rt_t,rm_t,epsr_t
	integer is_OBM_x
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	1	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in
c
       INTEGER M,N,IFLAG
       REAL X(N),FVEC(M)
       INTEGER I
       REAL TMP1,TMP2,TMP3,TMP4



c for unphysical inputs, set high penality

        if ((id_tool .eq. 7)  .or. (id_tool .eq. 6)) then
                 hd_max=28.00
        else
                 hd_max=20.00
        end if
c       
c	  hd_min=hd(1)-0.25
        hd_min=hdm          ! changed because don't know how hd(1) is passed, from where, what value

        

        bhc_para(1)= hd_x
        bhc_para(2)= rm_x
        bhc_para(3)= rt_x
        bhc_para(4)= epsr_x

!       update borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              bhc_para(j)=x(i)
           endif
        enddo
  
!       write(110,*) 'bhc_para', bhc_para(1:3)

       epsr_min=1.0
       epsr_max=100000.0


! Enforce bounds
       
        if(inv_switch_x(1).eq.1) then   !hd
          bhc_para(1)=max(bhc_para(1),hd_min) 
          bhc_para(1)=min(bhc_para(1),hd_max)
        endif 

        if(inv_switch_x(2).eq.1) then  !rm
          bhc_para(2)=max(bhc_para(2),rm_min) 
          bhc_para(2)=min(bhc_para(2),rm_max)
        endif 

        if(inv_switch_x(3).eq.1) then  !rt
          bhc_para(3)=max(bhc_para(3),rt_min) 
          bhc_para(3)=min(bhc_para(3),rt_max)
        endif 
        

        if(inv_switch_x(4).eq.1) then  !epsr
          bhc_para(4)=max(bhc_para(4),epsr_min) 
          bhc_para(4)=min(bhc_para(4),epsr_max)
        endif 
        
!        if((bhc_para(2).GT.bhc_para(3)) .and. (inv_switch_x(2).eq.1) ) then
!	               bhc_para(2)=bhc_para(3)/1.2
!	  endif
        
        hd_t=bhc_para(1)
        rm_t=bhc_para(2)
        rt_t=bhc_para(3)
        epsr_t=bhc_para(4)

!        write(110,*) hd_t,rt_t,rm_t

	          
		call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm_t,rt_t,hd_t,epsr_t,ibadtx)
!		write(800,*) 'psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx'
!		write(800,*) psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx

!		do i=1,ntr
!			if(idx_ps2(i) .ne. 0) y_m(idx_ps2(i))=psm2(i)
!		end do
!		do i=1,ntr
!			if(idx_ad2(i) .ne. 0) y_m(idx_ad2(i))=adm2(i)
!		end do

		y_m(1:5)=psm2_vec(1:5)
		y_m(6:10)=adm2_vec(1:5)
		y_m(11:15)=psm4_vec(1:5)
		y_m(16:20)=adm4_vec(1:5)
		
		index_all(1:5)=idx_ps2(1:5)
		index_all(6:10)=idx_ad2(1:5)
		index_all(11:15)=idx_ps4(1:5)
		index_all(16:20)=idx_ad4(1:5)

!		do i=1,m
!			fvec(i)=sqrt(w_in(i))*(y_in(i)-y_m(i))
!		end do

        fvec=0.0
        ii=0
        
        
		do i=1,max_m
			if (index_all(i).ne.0) then
			   ii=ii+1 
			   fvec(ii)=sqrt(w_in(i))*(y_in(i)-y_m(i))
			end if
		end do
		
!!	end if
	

!       update x from borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              x(i)=bhc_para(j)
           endif
        enddo


!!	write(800,*) 'idx_ps2,idx_ad2'
!!	write(800,*) idx_ps2,idx_ad2
!      write(800,*) 'y_m'
!      write(800,*) y_m
!
!      write(800,*) 'fvec'
!      write(800,*) fvec


      RETURN
        END

      SUBROUTINE FCN_bhcmodel2_2M(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
	include 'venus_ttmdivtab_def.inc'
c
	parameter (max_m=20)
      parameter (max_ms=10)
!	parameter (rt_min=0.05,rt_max=5000.0)
	parameter (rt_min=0.1,rt_max=1500.0)
	parameter (rm_min=0.015,rm_max=10.0)
      real*4 w_in(max_m),y_in(max_m),y_m(max_m), w_in2(max_ms),y_in2(max_ms) 
!	real*4 dpsm2(ntr),dadm2(ntr),dpsm4(ntr),dadm4(ntr)
!        real*4 psm2(ntr),adm2(ntr),rm0(18),error(18)
      real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr)
      real*4 hd_min, hd_max,epsr_min,epsr_max
      real*4 rt_x,rm_x,hd_x,epsr_x
      integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(max_m)
	integer n_unknow,m_data,k,j
	integer id_tool,No_x
	integer inv_switch_x(4)
	real*4 bhc_para(4) ! borehole model parameter (hd,rm,rt,ecc)
	real*4 hd_t,rt_t,rm_t,epsr_t
	integer is_OBM_x
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	1	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in
c
       INTEGER M,N,IFLAG
       REAL X(N),FVEC(M)
       INTEGER I
       REAL TMP1,TMP2,TMP3,TMP4



c for unphysical inputs, set high penality

       if ((id_tool .eq. 7)  .or. (id_tool .eq. 6)) then
                 hd_max=28.00
        else
                 hd_max=20.00
        end if
c       
c	  hd_min=hd(1)-0.25
        hd_min=hdm          ! changed because don't know how hd(1) is passed, from where, what value

        

        bhc_para(1)= hd_x
        bhc_para(2)= rm_x
        bhc_para(3)= rt_x
        bhc_para(4)= epsr_x

!       update borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              bhc_para(j)=x(i)
           endif
        enddo
  
!       write(110,*) 'bhc_para', bhc_para(1:3)

       epsr_min=1.0
       epsr_max=100000.0


! Enforce bounds
       
        if(inv_switch_x(1).eq.1) then   !hd
          bhc_para(1)=max(bhc_para(1),hd_min) 
          bhc_para(1)=min(bhc_para(1),hd_max)
        endif 

        if(inv_switch_x(2).eq.1) then  !rm
          bhc_para(2)=max(bhc_para(2),rm_min) 
          bhc_para(2)=min(bhc_para(2),rm_max)
        endif 

        if(inv_switch_x(3).eq.1) then  !rt
          bhc_para(3)=max(bhc_para(3),rt_min) 
          bhc_para(3)=min(bhc_para(3),rt_max)
        endif 
        

        if(inv_switch_x(4).eq.1) then  !epsr
          bhc_para(4)=max(bhc_para(4),epsr_min) 
          bhc_para(4)=min(bhc_para(4),epsr_max)
        endif 
        
!        if((bhc_para(2).GT.bhc_para(3)) .and. (inv_switch_x(2).eq.1) ) then
!	               bhc_para(2)=bhc_para(3)/1.2
!	  endif
        
        hd_t=bhc_para(1)
        rm_t=bhc_para(2)
        rt_t=bhc_para(3)
        epsr_t=bhc_para(4)

!        write(110,*) hd_t,rt_t,rm_t

	          
		call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm_t,rt_t,hd_t,epsr_t,ibadtx)
!		write(800,*) 'psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx'
!		write(800,*) psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx

!		do i=1,ntr
!			if(idx_ps2(i) .ne. 0) y_m(idx_ps2(i))=psm2(i)
!		end do
!		do i=1,ntr
!			if(idx_ad2(i) .ne. 0) y_m(idx_ad2(i))=adm2(i)
!		end do

		y_m(1:5)=psm2_vec(1:5)
		y_m(6:10)=adm2_vec(1:5)
		y_m(11:15)=psm4_vec(1:5)
		y_m(16:20)=adm4_vec(1:5)
		
		index_all(1:5)=idx_ps2(1:5)
		index_all(6:10)=idx_ad2(1:5)
		index_all(11:15)=idx_ps4(1:5)
		index_all(16:20)=idx_ad4(1:5)

!		do i=1,m
!			fvec(i)=sqrt(w_in(i))*(y_in(i)-y_m(i))
!		end do

        fvec=0.0
        ii=0
        
        
		do i=1,max_ms
			if (index_all(i).ne.0) then
			   ii=ii+1 
			   fvec(ii)=sqrt(w_in(i))*(y_in(i)-y_m(i))
			end if
		end do
		
!!	end if
	

!       update x from borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              x(i)=bhc_para(j)
           endif
        enddo


!!	write(800,*) 'idx_ps2,idx_ad2'
!!	write(800,*) idx_ps2,idx_ad2
!      write(800,*) 'y_m'
!      write(800,*) y_m
!
!      write(800,*) 'fvec'
!      write(800,*) fvec


      RETURN
       END
 
      SUBROUTINE FCN_bhcmodel2_400k(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
	include 'venus_ttmdivtab_def.inc'
c
	parameter (max_m=20)
      parameter (max_ms=10)
!	parameter (rt_min=0.05,rt_max=5000.0)
	parameter (rt_min=0.1,rt_max=1500.0)
	parameter (rm_min=0.015,rm_max=10.0)
      real*4 w_in(max_m),y_in(max_m),y_m(max_m), w_in2(max_ms),y_in2(max_ms) 
!	real*4 dpsm2(ntr),dadm2(ntr),dpsm4(ntr),dadm4(ntr)
!        real*4 psm2(ntr),adm2(ntr),rm0(18),error(18)
      real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr)
      real*4 hd_min, hd_max,epsr_min,epsr_max
      real*4 rt_x,rm_x,hd_x,epsr_x
      integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(max_m)
	integer n_unknow,m_data,k,j
	integer id_tool,No_x
	integer inv_switch_x(4)
	real*4 bhc_para(4) ! borehole model parameter (hd,rm,rt,ecc)
	real*4 hd_t,rt_t,rm_t,epsr_t
	integer is_OBM_x
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	1	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in
c
       INTEGER M,N,IFLAG
       REAL X(N),FVEC(M)
       INTEGER I
       REAL TMP1,TMP2,TMP3,TMP4



c for unphysical inputs, set high penality
         if ((id_tool .eq. 7)  .or. (id_tool .eq. 6)) then
                 hd_max=28.00
        else
                 hd_max=20.00
        end if
c       
c	  hd_min=hd(1)-0.25
        hd_min=hdm          ! changed because don't know how hd(1) is passed, from where, what value


        bhc_para(1)= hd_x
        bhc_para(2)= rm_x
        bhc_para(3)= rt_x
        bhc_para(4)= epsr_x

!       update borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              bhc_para(j)=x(i)
           endif
        enddo
  
!       write(110,*) 'bhc_para', bhc_para(1:3)

       epsr_min=1.0
       epsr_max=100000.0


! Enforce bounds
       
        if(inv_switch_x(1).eq.1) then   !hd
          bhc_para(1)=max(bhc_para(1),hd_min) 
          bhc_para(1)=min(bhc_para(1),hd_max)
        endif 

        if(inv_switch_x(2).eq.1) then  !rm
          bhc_para(2)=max(bhc_para(2),rm_min) 
          bhc_para(2)=min(bhc_para(2),rm_max)
        endif 

        if(inv_switch_x(3).eq.1) then  !rt
          bhc_para(3)=max(bhc_para(3),rt_min) 
          bhc_para(3)=min(bhc_para(3),rt_max)
        endif 
        

        if(inv_switch_x(4).eq.1) then  !epsr
          bhc_para(4)=max(bhc_para(4),epsr_min) 
          bhc_para(4)=min(bhc_para(4),epsr_max)
        endif 
        
!        if((bhc_para(2).GT.bhc_para(3)) .and. (inv_switch_x(2).eq.1) ) then
!	               bhc_para(2)=bhc_para(3)/1.2
!	  endif
        
        hd_t=bhc_para(1)
        rm_t=bhc_para(2)
        rt_t=bhc_para(3)
        epsr_t=bhc_para(4)

!        write(110,*) hd_t,rt_t,rm_t

	          
		call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm_t,rt_t,hd_t,epsr_t,ibadtx)
!		write(800,*) 'psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx'
!		write(800,*) psm2_vec,adm2_vec,rm_t,rt_t,hd_t,ibadtx

!		do i=1,ntr
!			if(idx_ps2(i) .ne. 0) y_m(idx_ps2(i))=psm2(i)
!		end do
!		do i=1,ntr
!			if(idx_ad2(i) .ne. 0) y_m(idx_ad2(i))=adm2(i)
!		end do

		y_m(1:5)=psm2_vec(1:5)
		y_m(6:10)=adm2_vec(1:5)
		y_m(11:15)=psm4_vec(1:5)
		y_m(16:20)=adm4_vec(1:5)
		
		index_all(1:5)=idx_ps2(1:5)
		index_all(6:10)=idx_ad2(1:5)
		index_all(11:15)=idx_ps4(1:5)
		index_all(16:20)=idx_ad4(1:5)

!		do i=1,m
!			fvec(i)=sqrt(w_in(i))*(y_in(i)-y_m(i))
!		end do

        fvec=0.0
        ii=0
        
        ! get the 400kHz measured - modeled
		do i=1,max_ms
			if (index_all(i+10).ne.0) then
			   ii=ii+1 
			   fvec(ii)=sqrt(w_in(i+10))*(y_in(i+10)-y_m(i+10))
			end if
		end do
		
!!	end if
	

!       update x from borehole parameters that were inverted
        i=0
        do j=1,4
           if(inv_switch_x(j).eq.1) then
              i=i+1
              x(i)=bhc_para(j)
           endif
        enddo


!!	write(800,*) 'idx_ps2,idx_ad2'
!!	write(800,*) idx_ps2,idx_ad2
!      write(800,*) 'y_m'
!      write(800,*) y_m
!
!      write(800,*) 'fvec'
!      write(800,*) fvec


      RETURN
      END

	subroutine initial_gbhc2_OBM(id_tool,ibadtx,hd_min,
	1	rm_in,y_in,w_in,m_in,nps,nad,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
	2	rt_ig,hd_ig,rm_ig,nig)
c
	parameter (ntr=5,ng=30)
	parameter (max_m=20)
	parameter (nrmax=9,nhdmax=28)
	real*4 rm_in(max_m),y_in(max_m),w_in(max_m)
	real*4 psm2_vec(ntr),adm2_vec(ntr),psm4_vec(ntr),adm4_vec(ntr),hd_min
	integer idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),ibadtx,id_tool
	real*4 rt_ig(ng),hd_ig(ng),rm_ig(ng)
	real*4 rt(nrmax),hd(nhdmax),el2_min(ng)
	real*4 el2(nrmax,nhdmax)
	real*4 decadem(nrmax)
	real*4 epsri
	real*4 tmp_rt,tmp
	integer ind_rt(20)

	data (ind_rt(i),i=1,20) / 15,20,5,10,14,19,4,9,13,18,3,8,12,17,2,7,11,16,1,6/

	data rm_OBM /1000.0/
	
	data (decadem(i),i=1,nrmax) /5.00 , 3.000,  2.044,  1.392,  0.949,
	1	0.646,  0.440,  0.300, 0.15 /
	

	data rmax / 5000.0 /
	data rmin / 0.1 /
c
	save all
      
      epsri=10.0 ! this initialization only applies for centered case
      nig=3

c Computing base of initial rt

      do i=1,20
		   if (abs(rm_in(ind_rt(i))+999.25).gt.1.0e-5) then
	         tmp_rt=rm_in(ind_rt(i))
	         goto 12
	     endif
	enddo
12          continue

      ir_st=1
      ir_end=9
      nr=ir_end-ir_st+1
	do i=1,nr
		j=ir_st+(i-1)
			rt(i)=tmp_rt*decadem(j)
			if (rt(i) .gt. rmax) rt(i)=rmax
			if (rt(i) .lt. rmin) rt(i)=rmin
	end do
c

      delta_hd=(20.0-hd_min)/(nhdmax-1)
!      delta_hd=(42.0-real(id_tool))/(nhdmax-1)
!      tmp=max((hd_min-4.0),real(id_tool))
!      delta_hd=(42.0-tmp)/(nhdmax-1)


	do i=1,nhdmax
		hd(i)=hd_min+(i-4)*delta_hd
	end do
c


	 do j=1,nr				!rt loop
        do k=1,nhdmax		!HD loop
	     call venus_forward_model(psm2_vec,adm2_vec,psm4_vec,adm4_vec,rm_OBM,rt(j),hd(k),epsri,ibadtx)
		 call el2norm(y_in,w_in,m_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,psm2_vec,adm2_vec,psm4_vec,adm4_vec,el2(j,k))  
	   end do
	 end do


c
      el2_min(1)=1.0e20
	 do j=1,nr				!rt loop
         do k=1,nhdmax			!HD loop
              if(el2(j,k) .le. el2_min(1)) then
                el2_min(1)=el2(j,k)
                rt_ig(1)=rt(j)
                hd_ig(1)=hd(k)
                rm_ig(1)=rm_OBM
                j1=j
                k1=k
              end if
         end do
	 end do


c
      el2_min(2)=1.0e20
       
	 do j=1,nr				!rt loop
         do k=1,nhdmax		!HD loop
              if ((j.eq.j1).and.(k.eq.k1)) 
	1          go to 2500
              if(el2(j,k) .le. el2_min(2)) then
                el2_min(2)=el2(j,k)
                rt_ig(2)=rt(j)
                hd_ig(2)=hd(k)
                rm_ig(2)=rm_OBM
                j2=j
                k2=k
              end if
2500	      continue
         end do
	 end do
       
c
      el2_min(3)=1.0e20
       
	 do j=1,nr				!rt loop
         do k=1,nhdmax			!HD loop
            if ( ((j.eq.j1).and.(k.eq.k1)) .or.
	1	((j.eq.j2).and.(k.eq.k2)))
	1        go to 3500
              if(el2(j,k) .le. el2_min(3)) then
                el2_min(3)=el2(j,k)
                rt_ig(3)=rt(j)
                hd_ig(3)=hd(k)
                rm_ig(3)=rm_OBM
              end if
3500	      continue
          end do
	 end do
       
c
       j=0
	return
      end   !end initial_gbhc2_OBM
c

       subroutine calculateQC(rmo,rto,hdo,epsro, erro,QC)
cc Cindy 06/11 Caculate QC 
cc Cindy 07/05 Caculate combined sensitivity for Rt, Rm and HD
!  keli sept. 25, 2013, cmpute standard deviation according to jacobian and Hessian 
! copied from keli, modify to accomodate dielectric inversion change
	parameter (max_m=20) !max. number of measurements
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
      real*4 y_m1(max_m),y_m2(max_m)
	real*4 sensit,rmo,rto,hdo,epsro, erro, QC
	real*4 sensit_hd(max_m),sensit_rt(max_m),sensit_rm(max_m),sensit_epsr(max_m)
	real*4 hdo_1,hdo_2,rmo_1,rmo_2,rto_1,rto_2,epsro_1,epsro_2
      real*4 psm2_1(ntr),adm2_1(ntr),psm2_2(ntr),adm2_2(ntr)
      real*4 psm4_1(ntr),adm4_1(ntr),psm4_2(ntr),adm4_2(ntr)

      real*4 delta
	INTEGER M,N,IFLAG
      INTEGER I
c
	real*4  w_in(max_m),y_in(max_m)
	integer n_unknow,m_data,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer ibadtx
	real*4  rm_t,hdm
	integer id_tool,No_x
	integer inv_switch_x(4)
      real*4  rt_x,rm_x,hd_x,epsr_x
      integer is_OBM_x
      real*4, allocatable:: Jacobian(:,:), Hessian(:,:), Hessianinv(:,:),uncertainty(:)
      integer nn,mm
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      integer ii,index_all(max_m)
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
     &	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in


      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)

         
c
c calculate sensitivity to hd	
      hdo_1=hdo+0.1;
      hdo_2=hdo-0.1;
      delta = 0.2;
      if (hdo .le. 0) then
         sensit_dh = 1.e-3
 !        return
      elseif (hdo_2 < 0) then
           hdo_2 = hdo
           delta = 0.1
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo_1,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
	do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo_2,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
      do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
	end do


c	
      ii=0
	do i=1,max_m
         if (index_all(i).ne.0) then
           ii=ii+1
           sensit_hd(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c  
c calculate sensitivity to rt	
      rto_1=rto+0.05;
      rto_2=rto-0.05;
      delta = 0.1;
      if (rto .le. 0) then
         sensit_rt = 1.e-3
!         return
      elseif (rto_2 < 0) then
           rto_2 = rto
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto_1,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
	do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
	end do

c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto_2,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do

c	
      ii=0
	do i=1,max_m
	   if(index_all(i).ne.0) then
	      ii=ii+1
            sensit_rt(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c    
c calculate sensitivity to rm	
      rmo_1=rmo+0.05;
      rmo_2=rmo-0.05;
      delta = 0.1;
      if (rmo .le. 0) then
         sensit_rm = 1.e-3
!         return
      elseif (rmo_2 < 0) then
           rmo_2 = rmo
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo_1,rto,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
	do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo_2,rto,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
      do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
	end do
c	
      ii=0
	do i=1,max_m
	   if (index_all(i).ne.0) then
	     ii=ii+1
           sensit_rm(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do


c calculate sensitivity to epsr	
C      eccmax=(hdo-real(id_tool))/2
      epsro_1=epsro+0.1;
      epsro_2=epsro-0.1;
      delta = 0.2;
      if (epsro .le. 0) then
         sensit_epsr = 1.e-3
!         return
      elseif (epsro_2 < 0) then
           epsro_2 = epsro
           delta = 0.1
      end if  
c         
 
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo,epsro_1,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
	do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo,epsro_2,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
      do i=1,ntr
		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
	end do
c	
      ii=0
	do i=1,max_m
         if(index_all(i).ne.0) then
            ii=ii+1
            sensit_epsr(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do

!     building Jacobian matrix
      nn=0
      do i=1,max_m
         if (index_all(i).ne.0) nn=nn+1
      enddo   

      mm=inv_switch_x(1)+inv_switch_x(2)+inv_switch_x(3)+inv_switch_x(4)
      
      
      allocate(Jacobian(nn,mm),Hessian(mm,mm),Hessianinv(mm,mm)) !Jacobian in the order of dh, rm, rt and ecc
      allocate(uncertainty(mm))

      Jacobian=0.0
      do i=1,nn
         ii=0
         
         if (inv_switch_x(1).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_hd(i)
         end if
            
         if (inv_switch_x(2).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rm(i)
         end if

         if (inv_switch_x(3).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rt(i)
         end if

         if (inv_switch_x(4).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_epsr(i)
         end if
      enddo

! Computing Hessian
!      Hessian=0.0
!
!	do i=1,mm
!	  do  j=1,mm
!	     do k=1,nn
!	        Hessian(i,j)=Hessian(i,j)+Jacobian(k,i)*Jacobian(k,j)
!           enddo
!        enddo
!      enddo             

      Hessian=  MATMUL(TRANSPOSE(Jacobian),Jacobian) 

    	call findINV(Hessian,Hessianinv,mm,errorflag)

      uncertainty=0.d0
          
      do i=1,mm
         uncertainty(i)=sqrt( abs(Hessianinv(i,i)) )
      enddo

      uncertainty_cutoff(1)=5.0 !cutoff value: 2 in for hd 
      uncertainty_cutoff(2)=0.3*rmo !cutoff value: 30% for rm
      uncertainty_cutoff(3)=0.3*rto !cutoff value: 30% for rt
      uncertainty_cutoff(4)=0.5*epsro      !cutoff value: 50% for epsr
      
      misfit_cutff=0.2            !cutoff value,20% misfit 
       
      ii=0

      sensit=0.0      
      do i=1,4
        if(inv_switch_x(i).eq.1) then
           ii=ii+1
           sensit=sensit+min(uncertainty(ii)/uncertainty_cutoff(i),1.0)
        endif
      enddo

      sensit=sensit/mm
      
      QC=sensit*erro/misfit_cutff


      if(QC.gt.1.0) QC=1.0

!      write(113,'(100f15.4)') erro, erro/misfit_cutff, uncertainty, uncertainty_cutoff, QC

 
!      if(inv_switch_x(1).eq.0) sensit_hd=1.0
!      if(inv_switch_x(2).eq.0) sensit_rm=1.0
!      if(inv_switch_x(3).eq.0) sensit_rt=1.0
!      if(inv_switch_x(4).eq.0) sensit_ecc=1.0
!
!
!c    
!c Calculate the total sensitivity
!      sensit = 0;
!      do i=1,2*ntr
!         sensit = sensit+sensit_hd(i)*sensit_rt(i)*sensit_rm(i)*sensit_ecc(i)
!	end do

      deallocate(Jacobian,Hessian,Hessianinv) !Jacobian in the order of dh, rm, rt and ecc
      deallocate(uncertainty)

	
	
       return
      end	
	
      
       subroutine calculateQC_2M(rmo,rto,hdo,epsro,erro,QC)
cc Cindy 06/11 Caculate QC 
cc Cindy 07/05 Caculate combined sensitivity for Rt, Rm and HD
!  keli sept. 25, 2013, cmpute standard deviation according to jacobian and Hessian 
! copied from keli, modify to accomodate dielectric inversion change
	parameter (max_m=20) !max. number of measurements
      parameter (max_ms=10) !max. number of measurements for one frequency 2M 10 channels
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
      real*4 y_m1(max_ms),y_m2(max_ms)
	real*4 sensit,rmo,rto,hdo,epsro, erro, QC
	real*4 sensit_hd(max_ms),sensit_rt(max_ms),sensit_rm(max_ms),sensit_epsr(max_ms)
	real*4 hdo_1,hdo_2,rmo_1,rmo_2,rto_1,rto_2,epsro_1,epsro_2
      real*4 psm2_1(ntr),adm2_1(ntr),psm2_2(ntr),adm2_2(ntr)
      real*4 psm4_1(ntr),adm4_1(ntr),psm4_2(ntr),adm4_2(ntr)

      real*4 delta
	INTEGER M,N,IFLAG
      INTEGER I
c
	real*4  w_in(max_m),y_in(max_m)
	integer n_unknow,m_data,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer ibadtx
	real*4  rm_t,hdm
	integer id_tool,No_x
	integer inv_switch_x(4)
      real*4  rt_x,rm_x,hd_x,epsr_x
      integer is_OBM_x
      real*4, allocatable:: Jacobian(:,:), Hessian(:,:), Hessianinv(:,:),uncertainty(:)
      integer nn,mm
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      integer ii,index_all(max_m)
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
     &	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in


      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)

         
c
c calculate sensitivity to hd	
      hdo_1=hdo+0.1;
      hdo_2=hdo-0.1;
      delta = 0.2;
      if (hdo .le. 0) then
         sensit_dh = 1.e-3
 !        return
      elseif (hdo_2 < 0) then
           hdo_2 = hdo
           delta = 0.1
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo_1,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo_2,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do


c	
      ii=0
	do i=1,max_ms
         if (index_all(i).ne.0) then
           ii=ii+1
           sensit_hd(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c  
c calculate sensitivity to rt	
      rto_1=rto+0.05;
      rto_2=rto-0.05;
      delta = 0.1;
      if (rto .le. 0) then
         sensit_rt = 1.e-3
!         return
      elseif (rto_2 < 0) then
           rto_2 = rto
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto_1,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do

c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto_2,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
c	end do

c	
      ii=0
	do i=1,max_ms
	   if(index_all(i).ne.0) then
	      ii=ii+1
            sensit_rt(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c    
c calculate sensitivity to rm	
      rmo_1=rmo+0.05;
      rmo_2=rmo-0.05;
      delta = 0.1;
      if (rmo .le. 0) then
         sensit_rm = 1.e-3
!         return
      elseif (rmo_2 < 0) then
           rmo_2 = rmo
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo_1,rto,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo_2,rto,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
	   if (index_all(i).ne.0) then
	     ii=ii+1
           sensit_rm(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do


c calculate sensitivity to epsr	
C      eccmax=(hdo-real(id_tool))/2
      epsro_1=epsro+0.1;
      epsro_2=epsro-0.1;
      delta = 0.2;
      if (epsro .le. 0) then
         sensit_epsr = 1.e-3
!         return
      elseif (epsro_2 < 0) then
           epsro_2 = epsro
           delta = 0.1
      end if  
c         
 
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo,epsro_1,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo,epsro_2,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
         if(index_all(i).ne.0) then
            ii=ii+1
            sensit_epsr(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do

!     building Jacobian matrix
      nn=0
      do i=1,max_ms
         if (index_all(i).ne.0) nn=nn+1
      enddo   

      mm=inv_switch_x(1)+inv_switch_x(2)+inv_switch_x(3)+inv_switch_x(4)
      
      
      allocate(Jacobian(nn,mm),Hessian(mm,mm),Hessianinv(mm,mm)) !Jacobian in the order of dh, rm, rt and ecc
      allocate(uncertainty(mm))

      Jacobian=0.0
      do i=1,nn
         ii=0
         
         if (inv_switch_x(1).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_hd(i)
         end if
            
         if (inv_switch_x(2).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rm(i)
         end if

         if (inv_switch_x(3).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rt(i)
         end if

         if (inv_switch_x(4).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_epsr(i)
         end if
      enddo

! Computing Hessian
!      Hessian=0.0
!
!	do i=1,mm
!	  do  j=1,mm
!	     do k=1,nn
!	        Hessian(i,j)=Hessian(i,j)+Jacobian(k,i)*Jacobian(k,j)
!           enddo
!        enddo
!      enddo             

      Hessian=  MATMUL(TRANSPOSE(Jacobian),Jacobian) 

    	call findINV(Hessian,Hessianinv,mm,errorflag)

      uncertainty=0.d0
          
      do i=1,mm
         uncertainty(i)=sqrt( abs(Hessianinv(i,i)) )
      enddo

      uncertainty_cutoff(1)=5.0 !cutoff value: 2 in for hd 
      uncertainty_cutoff(2)=0.3*rmo !cutoff value: 30% for rm
      uncertainty_cutoff(3)=0.3*rto !cutoff value: 30% for rt
      uncertainty_cutoff(4)=0.5*epsro      !cutoff value: 50% for epsr
      
      misfit_cutff=0.2            !cutoff value,20% misfit 
       
      ii=0

      sensit=0.0      
      do i=1,4
        if(inv_switch_x(i).eq.1) then
           ii=ii+1
           sensit=sensit+min(uncertainty(ii)/uncertainty_cutoff(i),1.0)
        endif
      enddo

      sensit=sensit/mm
      
      QC=sensit*erro/misfit_cutff


      if(QC.gt.1.0) QC=1.0

!      write(113,'(100f15.4)') erro, erro/misfit_cutff, uncertainty, uncertainty_cutoff, QC

 
!      if(inv_switch_x(1).eq.0) sensit_hd=1.0
!      if(inv_switch_x(2).eq.0) sensit_rm=1.0
!      if(inv_switch_x(3).eq.0) sensit_rt=1.0
!      if(inv_switch_x(4).eq.0) sensit_ecc=1.0
!
!
!c    
!c Calculate the total sensitivity
!      sensit = 0;
!      do i=1,2*ntr
!         sensit = sensit+sensit_hd(i)*sensit_rt(i)*sensit_rm(i)*sensit_ecc(i)
!	end do

      deallocate(Jacobian,Hessian,Hessianinv) !Jacobian in the order of dh, rm, rt and ecc
      deallocate(uncertainty)

	
	
       return
      end
      
       subroutine calculateQC_400k(rmo,rto,hdo,epsro,erro,QC)
cc Cindy 06/11 Caculate QC 
cc Cindy 07/05 Caculate combined sensitivity for Rt, Rm and HD
!  keli sept. 25, 2013, cmpute standard deviation according to jacobian and Hessian 
! copied from keli, modify to accomodate dielectric inversion change
	parameter (max_m=20) !max. number of measurements
      parameter (max_ms=10) !max. number of measurements for one frequency 2M 10 channels
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
      real*4 y_m1(max_ms),y_m2(max_ms)
	real*4 sensit,rmo,rto,hdo,epsro, erro, QC
	real*4 sensit_hd(max_ms),sensit_rt(max_ms),sensit_rm(max_ms),sensit_epsr(max_ms)
	real*4 hdo_1,hdo_2,rmo_1,rmo_2,rto_1,rto_2,epsro_1,epsro_2
      real*4 psm2_1(ntr),adm2_1(ntr),psm2_2(ntr),adm2_2(ntr)
      real*4 psm4_1(ntr),adm4_1(ntr),psm4_2(ntr),adm4_2(ntr)

      real*4 delta
	INTEGER M,N,IFLAG
      INTEGER I
c
	real*4  w_in(max_m),y_in(max_m)
	integer n_unknow,m_data,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer ibadtx
	real*4  rm_t,hdm
	integer id_tool,No_x
	integer inv_switch_x(4)
      real*4  rt_x,rm_x,hd_x,epsr_x
      integer is_OBM_x
      real*4, allocatable:: Jacobian(:,:), Hessian(:,:), Hessianinv(:,:),uncertainty(:)
      integer nn,mm
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      integer ii,index_all(max_m)
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
     &	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in


      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)

         
c
c calculate sensitivity to hd	
      hdo_1=hdo+0.1;
      hdo_2=hdo-0.1;
      delta = 0.2;
      if (hdo .le. 0) then
         sensit_dh = 1.e-3
 !        return
      elseif (hdo_2 < 0) then
           hdo_2 = hdo
           delta = 0.1
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo_1,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo_2,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do


c	
      ii=0
	do i=1,max_ms
         if (index_all(i).ne.0) then
           ii=ii+1
           sensit_hd(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c  
c calculate sensitivity to rt	
      rto_1=rto+0.05;
      rto_2=rto-0.05;
      delta = 0.1;
      if (rto .le. 0) then
         sensit_rt = 1.e-3
!         return
      elseif (rto_2 < 0) then
           rto_2 = rto
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto_1,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do

c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto_2,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
c	end do

c	
      ii=0
	do i=1,max_ms
	   if(index_all(i).ne.0) then
	      ii=ii+1
            sensit_rt(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c    
c calculate sensitivity to rm	
      rmo_1=rmo+0.05;
      rmo_2=rmo-0.05;
      delta = 0.1;
      if (rmo .le. 0) then
         sensit_rm = 1.e-3
!         return
      elseif (rmo_2 < 0) then
           rmo_2 = rmo
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo_1,rto,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo_2,rto,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
	   if (index_all(i).ne.0) then
	     ii=ii+1
           sensit_rm(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do


c calculate sensitivity to epsr	
C      eccmax=(hdo-real(id_tool))/2
      epsro_1=epsro+0.1;
      epsro_2=epsro-0.1;
      delta = 0.2;
      if (epsro .le. 0) then
         sensit_epsr = 1.e-3
!         return
      elseif (epsro_2 < 0) then
           epsro_2 = epsro
           delta = 0.1
      end if  
c         
 
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo,epsro_1,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo,epsro_2,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
         if(index_all(i).ne.0) then
            ii=ii+1
            sensit_epsr(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do

!     building Jacobian matrix
      nn=0
      do i=1,max_ms
         if (index_all(i).ne.0) nn=nn+1
      enddo   

      mm=inv_switch_x(1)+inv_switch_x(2)+inv_switch_x(3)+inv_switch_x(4)
      
      
      allocate(Jacobian(nn,mm),Hessian(mm,mm),Hessianinv(mm,mm)) !Jacobian in the order of dh, rm, rt and ecc
      allocate(uncertainty(mm))

      Jacobian=0.0
      do i=1,nn
         ii=0
         
         if (inv_switch_x(1).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_hd(i)
         end if
            
         if (inv_switch_x(2).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rm(i)
         end if

         if (inv_switch_x(3).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rt(i)
         end if

         if (inv_switch_x(4).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_epsr(i)
         end if
      enddo

! Computing Hessian
!      Hessian=0.0
!
!	do i=1,mm
!	  do  j=1,mm
!	     do k=1,nn
!	        Hessian(i,j)=Hessian(i,j)+Jacobian(k,i)*Jacobian(k,j)
!           enddo
!        enddo
!      enddo             

      Hessian=  MATMUL(TRANSPOSE(Jacobian),Jacobian) 

    	call findINV(Hessian,Hessianinv,mm,errorflag)

      uncertainty=0.d0
          
      do i=1,mm
         uncertainty(i)=sqrt( abs(Hessianinv(i,i)) )
      enddo

      uncertainty_cutoff(1)=5.0 !cutoff value: 2 in for hd 
      uncertainty_cutoff(2)=0.3*rmo !cutoff value: 30% for rm
      uncertainty_cutoff(3)=0.3*rto !cutoff value: 30% for rt
      uncertainty_cutoff(4)=0.5*epsro      !cutoff value: 50% for epsr
      
      misfit_cutff=0.2            !cutoff value,20% misfit 
       
      ii=0

      sensit=0.0      
      do i=1,4
        if(inv_switch_x(i).eq.1) then
           ii=ii+1
           sensit=sensit+min(uncertainty(ii)/uncertainty_cutoff(i),1.0)
        endif
      enddo

      sensit=sensit/mm
      
      QC=sensit*erro/misfit_cutff


      if(QC.gt.1.0) QC=1.0

!      write(113,'(100f15.4)') erro, erro/misfit_cutff, uncertainty, uncertainty_cutoff, QC

 
!      if(inv_switch_x(1).eq.0) sensit_hd=1.0
!      if(inv_switch_x(2).eq.0) sensit_rm=1.0
!      if(inv_switch_x(3).eq.0) sensit_rt=1.0
!      if(inv_switch_x(4).eq.0) sensit_ecc=1.0
!
!
!c    
!c Calculate the total sensitivity
!      sensit = 0;
!      do i=1,2*ntr
!         sensit = sensit+sensit_hd(i)*sensit_rt(i)*sensit_rm(i)*sensit_ecc(i)
!	end do

      deallocate(Jacobian,Hessian,Hessianinv) !Jacobian in the order of dh, rm, rt and ecc
      deallocate(uncertainty)

	
	
       return
      end	
	 
      subroutine calculate_std_hessian2M(rmo,rto,hdo,epsro,std_hd,std_rm,std_rt,std_epsr)
cc Cindy 06/11 Caculate QC 
cc Cindy 07/05 Caculate combined sensitivity for Rt, Rm and HD
!  keli sept. 25, 2013, cmpute standard deviation according to jacobian and Hessian 
! copied from keli, modify to accomodate dielectric inversion change
	parameter (max_m=20) !max. number of measurements
      parameter (max_ms=10) !max. number of measurements for one frequency 2M 10 channels
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
      real*4 y_m1(max_ms),y_m2(max_ms)
	real*4 sensit,rmo,rto,hdo,epsro, std_hd,std_rm,std_rt,std_epsr
      real*4 std(4)
	real*4 sensit_hd(max_ms),sensit_rt(max_ms),sensit_rm(max_ms),sensit_epsr(max_ms)
	real*4 hdo_1,hdo_2,rmo_1,rmo_2,rto_1,rto_2,epsro_1,epsro_2
      real*4 psm2_1(ntr),adm2_1(ntr),psm2_2(ntr),adm2_2(ntr)
      real*4 psm4_1(ntr),adm4_1(ntr),psm4_2(ntr),adm4_2(ntr)

      real*4 delta
	INTEGER M,N,IFLAG
      INTEGER I
c
	real*4  w_in(max_m),y_in(max_m)
	integer n_unknow,m_data,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer ibadtx
	real*4  rm_t,hdm
	integer id_tool,No_x
	integer inv_switch_x(4)
      real*4  rt_x,rm_x,hd_x,epsr_x
      integer is_OBM_x
      real*4, allocatable:: Jacobian(:,:), Hessian(:,:), Hessianinv(:,:),uncertainty(:)
      integer nn,mm
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      integer ii,index_all(max_m)
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
     &	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in


      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)

         
c
c calculate sensitivity to hd	
      hdo_1=hdo+0.1;
      hdo_2=hdo-0.1;
      delta = 0.2;
      if (hdo .le. 0) then
         sensit_dh = 1.e-3
 !        return
      elseif (hdo_2 < 0) then
           hdo_2 = hdo
           delta = 0.1
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo_1,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo_2,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do


c	
      ii=0
	do i=1,max_ms
         if (index_all(i).ne.0) then
           ii=ii+1
           sensit_hd(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c  
c calculate sensitivity to rt	
      rto_1=rto+0.05;
      rto_2=rto-0.05;
      delta = 0.1;
      if (rto .le. 0) then
         sensit_rt = 1.e-3
!         return
      elseif (rto_2 < 0) then
           rto_2 = rto
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto_1,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do

c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto_2,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
c	end do

c	
      ii=0
	do i=1,max_ms
	   if(index_all(i).ne.0) then
	      ii=ii+1
            sensit_rt(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c    
c calculate sensitivity to rm	
      rmo_1=rmo+0.5;
      rmo_2=rmo-0.5;
      delta = 1.0;
      if (rmo .le. 0) then
         sensit_rm = 1.e-3
!         return
      elseif (rmo_2 < 0) then
           rmo_2 = rmo
           delta = 0.5
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo_1,rto,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo_2,rto,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
	   if (index_all(i).ne.0) then
	     ii=ii+1
           sensit_rm(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do


c calculate sensitivity to epsr	
C      eccmax=(hdo-real(id_tool))/2
      epsro_1=epsro+50.0;
      epsro_2=epsro-50.0;
      delta = 100.0;
      if (epsro .le. 0) then
         sensit_epsr = 1.e-3
!         return
      elseif (epsro_2 < 0) then
           epsro_2 = epsro
           delta = 50.0
      end if  
c         
 
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo,epsro_1,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm2_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm2_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo,epsro_2,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
         if(index_all(i).ne.0) then
            ii=ii+1
            sensit_epsr(ii)=sqrt(w_in(i))*(y_m1(i)-y_m2(i))/delta
         end if
	end do

!     building Jacobian matrix
      nn=0
      do i=1,max_ms
         if (index_all(i).ne.0) nn=nn+1
      enddo   

      mm=inv_switch_x(1)+inv_switch_x(2)+inv_switch_x(3)+inv_switch_x(4)
      
      
      allocate(Jacobian(nn,mm),Hessian(mm,mm),Hessianinv(mm,mm)) !Jacobian in the order of dh, rm, rt and ecc
      allocate(uncertainty(mm))

      Jacobian=0.0
      do i=1,nn
         ii=0
         
         if (inv_switch_x(1).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_hd(i)
         end if
            
         if (inv_switch_x(2).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rm(i)
         end if

         if (inv_switch_x(3).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rt(i)
         end if

         if (inv_switch_x(4).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_epsr(i)
         end if
      enddo

! Computing Hessian
!      Hessian=0.0
!
!	do i=1,mm
!	  do  j=1,mm
!	     do k=1,nn
!	        Hessian(i,j)=Hessian(i,j)+Jacobian(k,i)*Jacobian(k,j)
!           enddo
!        enddo
!      enddo             

      Hessian=  MATMUL(TRANSPOSE(Jacobian),Jacobian) 

    	call findINV(Hessian,Hessianinv,mm,errorflag)

      uncertainty=0.d0
          
      do i=1,mm
         uncertainty(i)=sqrt( abs(Hessianinv(i,i)) )
      enddo

!      uncertainty_cutoff(1)=5.0 !cutoff value: 2 in for hd 
!      uncertainty_cutoff(2)=0.3*rmo !cutoff value: 30% for rm
!      uncertainty_cutoff(3)=0.3*rto !cutoff value: 30% for rt
!      uncertainty_cutoff(4)=0.5*epsro      !cutoff value: 50% for epsr
      
!      misfit_cutff=0.2            !cutoff value,20% misfit 
       std(1:4)=0.0
      ii=0

 !     sensit=0.0      
      do i=1,4
        if(inv_switch_x(i).eq.1) then
           ii=ii+1
           std(i)= uncertainty(ii)
!           sensit=sensit+min(uncertainty(ii)/uncertainty_cutoff(i),1.0)
        endif
      enddo
       std_hd=std(1)
       std_rm=std(2)
       std_rt=std(3)
       std_epsr=std(4)

!      sensit=sensit/mm
      
 !     QC=sensit*erro/misfit_cutff


!      if(QC.gt.1.0) QC=1.0

!      write(113,'(100f15.4)') erro, erro/misfit_cutff, uncertainty, uncertainty_cutoff, QC

 
!      if(inv_switch_x(1).eq.0) sensit_hd=1.0
!      if(inv_switch_x(2).eq.0) sensit_rm=1.0
!      if(inv_switch_x(3).eq.0) sensit_rt=1.0
!      if(inv_switch_x(4).eq.0) sensit_ecc=1.0
!
!
!c    
!c Calculate the total sensitivity
!      sensit = 0;
!      do i=1,2*ntr
!         sensit = sensit+sensit_hd(i)*sensit_rt(i)*sensit_rm(i)*sensit_ecc(i)
!	end do

      deallocate(Jacobian,Hessian,Hessianinv) !Jacobian in the order of dh, rm, rt and ecc
      deallocate(uncertainty)

	
	
       return
      end
      
       subroutine calculate_std_hessian400k(rmo,rto,hdo,epsro,std_hd,std_rm,std_rt,std_epsr)
cc Cindy 06/11 Caculate QC 
cc Cindy 07/05 Caculate combined sensitivity for Rt, Rm and HD
!  keli sept. 25, 2013, cmpute standard deviation according to jacobian and Hessian 
! copied from keli, modify to accomodate dielectric inversion change
	parameter (max_m=20) !max. number of measurements
      parameter (max_ms=10) !max. number of measurements for one frequency 2M 10 channels
	parameter (max_n=4)  !max. number of unknown parameters
	parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
      real*4 y_m1(max_ms),y_m2(max_ms)
	real*4 sensit,rmo,rto,hdo,epsro, std_hd,std_rm,std_rt,std_epsr
      real*4 std(4)
	real*4 sensit_hd(max_ms),sensit_rt(max_ms),sensit_rm(max_ms),sensit_epsr(max_ms)
	real*4 hdo_1,hdo_2,rmo_1,rmo_2,rto_1,rto_2,epsro_1,epsro_2
      real*4 psm2_1(ntr),adm2_1(ntr),psm2_2(ntr),adm2_2(ntr)
      real*4 psm4_1(ntr),adm4_1(ntr),psm4_2(ntr),adm4_2(ntr)

      real*4 delta
	INTEGER M,N,IFLAG
      INTEGER I
c
	real*4  w_in(max_m),y_in(max_m)
	integer n_unknow,m_data,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr)
	integer ibadtx
	real*4  rm_t,hdm
	integer id_tool,No_x
	integer inv_switch_x(4)
      real*4  rt_x,rm_x,hd_x,epsr_x
      integer is_OBM_x
      real*4, allocatable:: Jacobian(:,:), Hessian(:,:), Hessianinv(:,:),uncertainty(:)
      integer nn,mm
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      integer ii,index_all(max_m)
c
	common /data_in_bhc/ y_in,n_unknow,m_data,idx_ps2,idx_ad2,idx_ps4,idx_ad4,
     &	ibadtx,rm_t,hdm,id_tool,No_x,rt_x,rm_x,hd_x,epsr_x,inv_switch_x,is_OBM_x
	common /weight_in/ w_in


      index_all(1:ntr)= idx_ps2(1:ntr)
      index_all(ntr+1:2*ntr)= idx_ad2(1:ntr)
      index_all(2*ntr+1:3*ntr)= idx_ps4(1:ntr)
      index_all(3*ntr+1:4*ntr)= idx_ad4(1:ntr)

         
c
c calculate sensitivity to hd	
      hdo_1=hdo+0.1;
      hdo_2=hdo-0.1;
      delta = 0.2;
      if (hdo .le. 0) then
         sensit_dh = 1.e-3
 !        return
      elseif (hdo_2 < 0) then
           hdo_2 = hdo
           delta = 0.1
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo_1,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo_2,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do


c	
      ii=0
	do i=1,max_ms
         if (index_all(i).ne.0) then
           ii=ii+1
           sensit_hd(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c  
c calculate sensitivity to rt	
      rto_1=rto+0.05;
      rto_2=rto-0.05;
      delta = 0.1;
      if (rto .le. 0) then
         sensit_rt = 1.e-3
!         return
      elseif (rto_2 < 0) then
           rto_2 = rto
           delta = 0.05
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto_1,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do

c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto_2,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm2_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm2_2(i)
c	end do

c	
      ii=0
	do i=1,max_ms
	   if(index_all(i).ne.0) then
	      ii=ii+1
            sensit_rt(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do
c    
c calculate sensitivity to rm	
      rmo_1=rmo+0.5;
      rmo_2=rmo-0.5;
      delta = 1.0;
      if (rmo .le. 0) then
         sensit_rm = 1.e-3
!         return
      elseif (rmo_2 < 0) then
           rmo_2 = rmo
           delta = 0.5
      end if  
c      
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo_1,rto,hdo,epsro,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo_2,rto,hdo,epsro,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
	   if (index_all(i).ne.0) then
	     ii=ii+1
           sensit_rm(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do


c calculate sensitivity to epsr	
C      eccmax=(hdo-real(id_tool))/2
      epsro_1=epsro+50.0;
      epsro_2=epsro-50.0;
      delta = 100.0;
      if (epsro .le. 0) then
         sensit_epsr = 1.e-3
!         return
      elseif (epsro_2 < 0) then
           epsro_2 = epsro
           delta = 50.0
      end if  
c         
 
      call venus_forward_model(psm2_1,adm2_1,psm4_1,adm4_1,rmo,rto,hdo,epsro_1,ibadtx)
	do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m1(idx_ps2(i))=psm4_1(i)
      end do
      do i=1,ntr
		  if(idx_ad2(i) .ne. 0) y_m1(idx_ad2(i))=adm4_1(i)
	end do
c	do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m1(idx_ps4(i))=psm4_1(i)
c      end do
c      do i=1,ntr
c		  if(idx_ad4(i) .ne. 0) y_m1(idx_ad4(i))=adm4_1(i)
c	end do
c
      call venus_forward_model(psm2_2,adm2_2,psm4_2,adm4_2,rmo,rto,hdo,epsro_2,ibadtx)
      do i=1,ntr
		  if(idx_ps2(i) .ne. 0) y_m2(idx_ps2(i))=psm4_2(i)
	end do
	do i=1,ntr
	    if(idx_ad2(i) .ne. 0) y_m2(idx_ad2(i))=adm4_2(i)
	end do
c      do i=1,ntr
c		  if(idx_ps4(i) .ne. 0) y_m2(idx_ps4(i))=psm4_2(i)
c	end do
c	do i=1,ntr
c	    if(idx_ad4(i) .ne. 0) y_m2(idx_ad4(i))=adm4_2(i)
c	end do
c	
      ii=0
	do i=1,max_ms
         if(index_all(i).ne.0) then
            ii=ii+1
            sensit_epsr(ii)=sqrt(w_in(i+10))*(y_m1(i)-y_m2(i))/delta
         end if
	end do

!     building Jacobian matrix
      nn=0
      do i=1,max_ms
         if (index_all(i).ne.0) nn=nn+1
      enddo   

      mm=inv_switch_x(1)+inv_switch_x(2)+inv_switch_x(3)+inv_switch_x(4)
      
      
      allocate(Jacobian(nn,mm),Hessian(mm,mm),Hessianinv(mm,mm)) !Jacobian in the order of dh, rm, rt and ecc
      allocate(uncertainty(mm))

      Jacobian=0.0
      do i=1,nn
         ii=0
         
         if (inv_switch_x(1).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_hd(i)
         end if
            
         if (inv_switch_x(2).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rm(i)
         end if

         if (inv_switch_x(3).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_rt(i)
         end if

         if (inv_switch_x(4).eq.1) then
            ii=ii+1
            Jacobian(i,ii)=sensit_epsr(i)
         end if
      enddo

! Computing Hessian
!      Hessian=0.0
!
!	do i=1,mm
!	  do  j=1,mm
!	     do k=1,nn
!	        Hessian(i,j)=Hessian(i,j)+Jacobian(k,i)*Jacobian(k,j)
!           enddo
!        enddo
!      enddo             

      Hessian=  MATMUL(TRANSPOSE(Jacobian),Jacobian) 

    	call findINV(Hessian,Hessianinv,mm,errorflag)

      uncertainty=0.d0
          
      do i=1,mm
         uncertainty(i)=sqrt( abs(Hessianinv(i,i)) )
      enddo

!      uncertainty_cutoff(1)=5.0 !cutoff value: 2 in for hd 
!      uncertainty_cutoff(2)=0.3*rmo !cutoff value: 30% for rm
!      uncertainty_cutoff(3)=0.3*rto !cutoff value: 30% for rt
!      uncertainty_cutoff(4)=0.5*epsro      !cutoff value: 50% for epsr
      
!      misfit_cutff=0.2            !cutoff value,20% misfit 
       
      ii=0
       std(1:4)=0.0
!      sensit=0.0      
      do i=1,4
        if(inv_switch_x(i).eq.1) then
           ii=ii+1
           std(i)= uncertainty(ii)
!           sensit=sensit+min(uncertainty(ii)/uncertainty_cutoff(i),1.0)
        endif
      enddo
       std_hd=std(1)
       std_rm=std(2)
       std_rt=std(3)
       std_epsr=std(4)

!      sensit=sensit/mm
      
!      QC=sensit*erro/misfit_cutff


!      if(QC.gt.1.0) QC=1.0

!      write(113,'(100f15.4)') erro, erro/misfit_cutff, uncertainty, uncertainty_cutoff, QC

 
!      if(inv_switch_x(1).eq.0) sensit_hd=1.0
!      if(inv_switch_x(2).eq.0) sensit_rm=1.0
!      if(inv_switch_x(3).eq.0) sensit_rt=1.0
!      if(inv_switch_x(4).eq.0) sensit_ecc=1.0
!
!
!c    
!c Calculate the total sensitivity
!      sensit = 0;
!      do i=1,2*ntr
!         sensit = sensit+sensit_hd(i)*sensit_rt(i)*sensit_rm(i)*sensit_ecc(i)
!	end do

      deallocate(Jacobian,Hessian,Hessianinv) !Jacobian in the order of dh, rm, rt and ecc
      deallocate(uncertainty)

	
	
       return
      end	
     
      
      subroutine calculateQC_V2(dho,rmo,rto,epsro,erro,std_hd,std_rm,std_rt,std_epsr,QC1)
!  keli sept. 25, 2013, cmpute QC based onstandard deviation and data mismatch   
!     parameter (max_m=20) !max. number of measurements
!     parameter (max_n=4)  !max. number of unknown parameters
!     parameter (ntr=5,ntf=25) !5 mbhc spacings,25-point transform
!      real*4 y_m1(max_m),y_m2(max_m)

      !input
      real*4 dho,rmo,rto,epsro,erro,std_hd,std_rm,std_rt,std_epsr
      
      !output
      real*4 QC1
      
      !local
      real*4 sensit
      INTEGER i,j,k,ii
c
!     integer inv_switch_x(4)
      real*4 uncertainty_cutoff(4), misfit_cutff
      integer errorflag
      real*4 std(4)
      real*4 weight(4), weightall
      
      std(1)=std_hd
      std(2)=std_rm
      std(3)=std_rt
      std(4)=std_epsr       


      uncertainty_cutoff(1)=20.0 !cutoff value: 2 in for hd 
      uncertainty_cutoff(2)= rmo !cutoff value: 30% for rm
      uncertainty_cutoff(3)=0.4*rto !cutoff value: 30% for rt
      uncertainty_cutoff(4)= 0.4*epsro      !cutoff value: 25 for epsr, but it could be bigger?
      
      misfit_cutff=0.5            !cutoff value,20% misfit 
      
      QC1=min(erro/misfit_cutff,1.0)

      do i=1,4
           ratio=std(i)/uncertainty_cutoff(i)
          QC1=max(min(ratio,1.0),QC1)
      enddo
      
     
      return
      end   

        SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
        IMPLICIT NONE
        !Declarations
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
        REAL*4, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
        REAL*4, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
        
        LOGICAL :: FLAG = .TRUE.
        INTEGER :: i, j, k, l
        REAL*4 :: m
        REAL*4, DIMENSION(n,2*n) :: augmatrix !augmented matrix
        
        !Augment input matrix with an identity matrix
        DO i = 1, n
                DO j = 1, 2*n
                        IF (j <= n ) THEN
                                augmatrix(i,j) = matrix(i,j)
                        ELSE IF ((i+n) == j) THEN
                                augmatrix(i,j) = 1
                        Else
                                augmatrix(i,j) = 0
                        ENDIF
                END DO
        END DO
        
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
                IF (augmatrix(k,k) == 0) THEN
                    FLAG = .FALSE.
                    DO i = k+1, n
                         IF (augmatrix(i,k) /= 0) THEN
                         DO j = 1,2*n
                          augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                         END DO
                         FLAG = .TRUE.
                         EXIT
                         ENDIF
                           IF (FLAG .EQV. .FALSE.) THEN
!                                  PRINT*, "Matrix is non - invertible"
                                      inverse = 0
                                        errorflag = -1
                                        return
                            ENDIF
                        END DO
                ENDIF
                DO j = k+1, n                        
                        m = augmatrix(j,k)/augmatrix(k,k)
                      DO i = k, 2*n
                      augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
                      END DO
                END DO
        END DO
        
        !Test for invertibility
        DO i = 1, n
                IF (augmatrix(i,i) == 0) THEN
!                        PRINT*, "Matrix is non - invertible"
                        inverse = 0
                        errorflag = -1
                        return
                ENDIF
        END DO
        
        !Make diagonal elements as 1
        DO i = 1 , n
                m = augmatrix(i,i)
                DO j = i , (2 * n)                                
                           augmatrix(i,j) = (augmatrix(i,j) / m)
                END DO
        END DO
        
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
                DO i =1, k
                m = augmatrix(i,k+1)
                 DO j = k, (2*n)
                  augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
                 END DO
                END DO
        END DO                                
        
        !store answer
        DO i =1, n
                DO j = 1, n
                        inverse(i,j) = augmatrix(i,j+n)
                END DO
        END DO
        errorflag = 0
	END SUBROUTINE FINDinv 	 