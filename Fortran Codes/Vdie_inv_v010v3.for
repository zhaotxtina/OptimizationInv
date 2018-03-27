

!      subroutine Vdie_inv_v001(id_tool,bit_size,mode,   ! tool ID, bit_size, inversion mode
!     &     nch,index_channel,weight_channel,          ! number of input channels, channels index and weighting 
!     &     total_number,depth,meas,nch_desire,i_auto, ! No. measurement points, depth and measurement channels
!     &     dh_input,rm_input,rt_input,epsr_input,      ! model parameters in case they are not inverted
!     &     Data_FilePath,path_size,                   ! location of data files
!     &     is_OBM,misfit_epsr,                         ! OBM indicator, misfit threshhold for Epsr inversion?
!     &     dh_bhc,rm_bhc,rt_bhc,epsr_bhc,err_bhc,QC)   ! inversion results

      subroutine Vdie_inv_v003(id_tool,bit_size,mode,   ! tool ID, bit_size, inversion mode
     &     nch,index_channel,weight_channel,          ! number of input channels, channels index and weighting 
     &     total_number,depth,meas,nch_desire,i_auto, ! No. measurement points, depth and measurement channels
     &     dh_input,rm_input,rt_input,epsr_input,      ! model parameters in case they are not inverted
     &     Data_FilePath,path_size,                   ! location of data files
     &     is_OBM,misfit_epsr,                         ! OBM indicator, misfit threshhold for Epsr inversion?
     &     dh_bhc2,rm_bhc2,rt_bhc2,epsr_bhc2,err_bhc2,QC2,   ! inversion results
     &     dh_bhc4,rm_bhc4,rt_bhc4,epsr_bhc4,err_bhc4,QC4)   ! inversion results
!     &     std_hd2,std_rm2, std_rt2,std_epsr2,   ! uncertainty results
!     &     std_hd4,std_rm4, std_rt4,std_epsr4)   ! uncertainty results
      
      
cDEC$ ATTRIBUTES DLLEXPORT :: Vdie_inv_v003
! parameters
      integer,parameter:: nd=100000
	integer,parameter:: ntr=5,nframe=1,nevent=5,ntf2=29,ntf4=31 !5 mbhc spacings
	real*4, parameter:: abst = -999.25 	!absent value
	integer, parameter::max_m=20 

cccccccccccccc Input 
	integer id_tool,i400k
      real*4 rm,bit_size 
	integer mode ! inversion mode
      integer nch
      integer index_channel(nch)
      real*4  weight_channel(nch)
      integer total_number
	real*4 depth(total_number)
	real*4 meas(total_number,nch) ! channel orders defined by index_channel

	real*4  arc_dhsw_ver
      integer path_size
	CHARACTER Data_FilePath*path_size
!	CHARACTER*200 Data_FilePath
	integer is_OBM ! OBM indicator
	real*4  misfit_epsr !Epsr inversion threshhold??
	real*4 dh_input(total_number),rm_input(total_number),rt_input(total_number),epsr_input(total_number)
	character*300 filename
	
      integer nch_desire  !selecting 10 channels for inversion, for each frequency 
      integer i_auto


!     output 
      
!	real*4 dh_bhc(total_number),rm_bhc(total_number),rt_bhc(total_number),epsr_bhc(total_number),
!     &       err_bhc(total_number),QC(total_number)

      real*4 dh_bhc2(total_number),rm_bhc2(total_number),rt_bhc2(total_number),epsr_bhc2(total_number),
     &       err_bhc2(total_number),QC2(total_number)
      real*4 std_hd2(total_number),std_rm2(total_number), std_rt2(total_number),std_epsr2(total_number)
      real*4 dh_bhc4(total_number),rm_bhc4(total_number),rt_bhc4(total_number),epsr_bhc4(total_number),
     &       err_bhc4(total_number),QC4(total_number)
      real*4 std_hd4(total_number),std_rm4(total_number), std_rt4(total_number),std_epsr4(total_number)


!     intermediate variables
      ! this portion is probably for inverse transform from Rt to AD/PS measurement
	real*4 rps2(ntr,nd),rad2(ntr,nd),rps4(ntr,nd),rad4(ntr,nd)  ! ntr=5, nd=100000
!	real*4 ps2_am(ntr),ad2_am(ntr),ps4_am(ntr),ad4_am(ntr)
!     real*4 rps2_tmp(ntr,nd),rad2_tmp(ntr,nd)
	real*4 p2(2*ntr,nd),a2(2*ntr,nd),p4(2*ntr,nd),a4(2*ntr,nd) !needed for dhsw_ver.6.3 only, ??? don't know yet
      integer icontig
      real*4 bs(nframe),t0,rm0
      real*4 rto(nframe),hdo(nframe),erro(nframe),rmo(nframe),epsro(nframe),QCo(nframe) ! nframe=1
	real*4 ps2(ntr,nd),ad2(ntr,nd),ps4(ntr,nd),ad4(ntr,nd)
	real*4 psm2(ntr,nd),adm2(ntr,nd),psm4(ntr,nd),adm4(ntr,nd)
	real*4 psm2p(ntr),adm2p(ntr)
	real*4 ru2(ntf2),ru4(ntf4),p_airt2(2*ntr),p_airt4(2*ntr)    ! ntf2=29(2M), ntf4=31(400k),  for resistivity transform table points
	real*4 rmlog(1),dipi(nd),pu2(ntf2,ntr*2),pu4(ntf4,ntr*2)
	real*4 rpsb(ntr,nd),radb(ntr,nd),fraction
	integer Nx, inv_switch(4), inv_switch_tmp(4) ! switch for whether inverting dh, rm, rt and epsr, 1 for invert, 0 for not inverting.

!	real*4 measi_raw(max_m,nframe)
	
      !  I haven't figured out the function of those parameters, seem from post processing of ARCwizard, the flags from 4 answer products
	integer iflag(nd,nevent),itemp(nd,nevent)    ! nevent=5
	real*4 frh(nd),frv(nd),frxo(nd),frt(nd),fdi(nd)         ! don't know those variables will be used or not, are they from ARCWizard?
	real*4 frt_bhc(nd),fhd_bhc(nd),ferr_ani(nd),ferr_inv(nd)
	real*4 ferr_die(nd),ferr_bhc(nd),tmp(nd),tmplog(nd)
	real*4 rxo_low(nd),rw_low(nd)
	real*8 test(855),temp(20)
	integer i,n,j,iz,nchar,istart,iend,ich
	integer i_init

      real*4 psm2i(ntr,nframe),adm2i(ntr,nframe),psm4i(ntr,nframe),adm4i(ntr,nframe)
      real*4 psm2_vec(5),adm2_vec(5),psm4_vec(5),adm4_vec(5)
	real*4 rps2_vec(5),rad2_vec(5),rps4_vec(5),rad4_vec(5)
	real*4, parameter:: adps2_min=0.10000001, adps4_min=0.050000001
	real*4, parameter:: adps2_max=1999.000000, adps4_max=1999.000000
      real*4 Rapp_max(20),Rapp_min(20)

!temperarily loading in raw measurement signal, for testing purpose only , will be removed 
      real*4 resp(nd,24)  ! raw signal measurements in db/deg, same convention as meas 
      real*4 hdi,rmi,rti,epsri
      real*4 tmpv(24)
!      integer,parameter:: i_test_transform=0
      
! local variables for removing clipped channels
      real*4,allocatable:: Rapp_vec(:)      
      integer nch_short,nch_out,ii

      integer,allocatable:: ind_short(:)
      real*4, allocatable:: Rapp_short(:), weight_short(:)

! locan variables for channel selection
      integer, allocatable:: ind_select_local(:)
      
      !revised index and weighting vector used for inversion (removing clipped channels, and select the best subset)
      integer nch_select
      real*4, allocatable:: weight_forinv(:),Rapp_forinv(:)
      integer, allocatable::ind_forinv(:) 
      
            
      allocate (Rapp_vec(nch))
      allocate (ind_select_local(nch_desire))   !selected nch_out (=nch_desire) channels ind_select(1:nch_out), ind_select(nch_out+1:nch) are not used

      
        ! defining min and max vector, following eCAl standard order
      Rapp_min(1:10)=adps2_min   ! defined as 0.1
      Rapp_min(11:20)=adps4_min  ! 0.05
      Rapp_max(1:10)=adps2_max   ! defined as 1999.0
      Rapp_max(11:20)=adps4_max  !  1999.0


	bs(1)=bit_size
	inv_switch(1)=0  ! inverting for dh
	inv_switch(2)=0  ! inverting for rm
	inv_switch(3)=0  ! inverting for rt
	inv_switch(4)=0  ! inverting for epsr


      if (mode .eq. 1) then      !rt and epsr  only
          inv_switch(3)=1
          inv_switch(4)=1
      elseif (mode .eq. 2) then  !invert for dh and rm,rt,epsr
          inv_switch(1)=1
          inv_switch(2)=1
          inv_switch(3)=1
          inv_switch(4)=1

      endif
      
      if(is_OBM.eq.1) then

          if(mode.eq.2) then
             write(*,*) 'Rm cannot be inverted for OBM'
             return
           endif

          elseif (mode .eq. 1) then      !rt and epsr  only
              inv_switch(3)=1
              inv_switch(4)=1

      endif

      Nx = sum(inv_switch)

	rps2 = abst 
	rad2 = abst
	rps4 = abst 
	rad4 = abst
      
      ! Assign (available) rps2,rad2,rps4,rad4 values from meas
	! Default order of measuremetn channels(channel index): ps2,ad2,ps4,ad4 
      
      do i=1,total_number
         tmpv=abst
         do j=1,nch
            tmpv(index_channel(j))=meas(i,j)
         enddo
         
         do j=1,5
            rad4(j,i)=tmpv(15+j)
            rps4(j,i)=tmpv(10+j)
            rad2(j,i)=tmpv(5+j)
            rps2(j,i)=tmpv(j)
         enddo
      
         ! Assign the clipped channels to be absent  ****differnt cutoff for 2MHz and 400kHz
         do j=1,5
            if(rad4(j,i).le.adps4_min) rad4(j,i)=abst
            if(rad4(j,i).ge.adps4_max) rad4(j,i)=abst
            if(rad2(j,i).le.adps2_min) rad2(j,i)=abst
            if(rad2(j,i).ge.adps2_max) rad2(j,i)=abst
            
            if(rps4(j,i).le.adps4_min) rps4(j,i)=abst
            if(rps4(j,i).ge.adps4_max) rps4(j,i)=abst
            if(rps2(j,i).le.adps2_min) rps2(j,i)=abst
            if(rps2(j,i).ge.adps2_max) rps2(j,i)=abst
            
         enddo
      enddo  
      
c      call arc_load_ttmbhctbl(id_tool,Data_FilePath,path_size)    !load table for (Rt,Rm,HD) --> PS and AD
      call venus_load_ttmdivtbl(id_tool,Data_FilePath,path_size)    !load table for (Rt,Rm,HD,epsr) --> PS and AD
      
      !! debug success Jan 21
      
        i400k = 0
        icontig = 6
c        arc_dhsw_ver = 6.0
        rm0 = 0.05
        t0 = 127
	  arc_dhsw_ver = 7.3 

!        if(arc_dhsw_ver.lt.6.3)then
!          call get_arc_ps_ad(ps2,ad2,ps4,ad4,p2,a2,p4,a4,   !Compute PS and AD from P and A 
!     1	id_tool,arc_dhsw_ver,i400k,nd)
!        end if
      

        i_init = 0
       
        ! probably testing only, for the following one point, inverse transform
       rps2_vec(1:5)=rps2(1:5,1)
	 rad2_vec(1:5)=rad2(1:5,1)
	 rps4_vec(1:5)=rps4(1:5,1)
	 rad4_vec(1:5)=rad4(1:5,1)

      !      converting Rapp back to db/deg so that they can be used in the inversion
         call inv_ra_1n(psm2_vec,adm2_vec,psm4_vec,adm4_vec,
	1	id_tool,i_init,abst,rps2_vec,rad2_vec,rps4_vec,rad4_vec,icontig)


       psm2(1:5,1)=psm2_vec(1:5)
       adm2(1:5,1)=adm2_vec(1:5)
       psm4(1:5,1)=psm4_vec(1:5)
       adm4(1:5,1)=adm4_vec(1:5)
       
        ! first convert the whole input measurement data (Rapp) ==> AD/PS

       do i=1,total_number
      
       rps2_vec(1:5)=rps2(1:5,i)
	 rad2_vec(1:5)=rad2(1:5,i)
	 rps4_vec(1:5)=rps4(1:5,i)
	 rad4_vec(1:5)=rad4(1:5,i)

!      converting compensated non-borehle corrected Rapp (rps2_vec,rad2_vec,rps4_vec,rad4_vec)  back to db/deg (psm2_vec,adm2_vec,psm4_vec,adm4_vec)so that they can be used in the inversion

	   call inv_ra_1n(psm2_vec,adm2_vec,psm4_vec,adm4_vec,              
	1	id_tool,i-1,abst,rps2_vec,rad2_vec,rps4_vec,rad4_vec, icontig)         !    rps2(1,i),rad2(1,i),rps4(1,i),rad4(1,i)
       
 
      ! psm2i,adm2i,psm4i,adm4i  are the PS/AD measurements now, from Rapp to measurement
       psm2i(1:5,1)=psm2_vec(1:5)
       adm2i(1:5,1)=adm2_vec(1:5)
       psm4i(1:5,1)=psm4_vec(1:5)
       adm4i(1:5,1)=adm4_vec(1:5)
      
c   
      
       nch_select=nch
       
       
!       if (i_auto.eq.1) then   ! select channels, remove clipped channels,  use i_auto=0 to avoid this part
           
!****************************************
! Removing clipped channels and absent channels
!****************************************
!         Rapp_vec(1:nch)=meas(i,1:nch)
        
!         nch_short=0
!         do j=1,nch
!           if(Rapp_vec(j).lt.Rapp_max(index_channel(j)) .and. Rapp_vec(j).gt. Rapp_min(index_channel(j)) )  nch_short=nch_short+1
!         enddo

!         allocate(Rapp_short(nch_short), weight_short(nch_short))
!         allocate(ind_short(nch_short))

!         Rapp_short=abst 
!         ind_short=0
!         weight_short=0.0

!         ii=0        
!         do j=1,nch

!           if(Rapp_vec(j).lt.Rapp_max(index_channel(j)) .and. Rapp_vec(j).gt. Rapp_min(index_channel(j)) ) then
!             ii=ii+1
!             Rapp_short(ii)=Rapp_vec(j)
!             ind_short(ii) = index_channel(j)
!             weight_short(ii)=weight_channel(j)
!           endif
!         enddo

        
! finding the best nch_desire channels to be used in the inversion
! *********************************************************
! Selecting best k channels from list of cond_short(1:nch_short)
! local index for cond_short: 1,2,3,..., nch_short 
! *********************************************************
        
!       ! based on channel list Rapp_short and ind_short, finding the combination of nch_select channels that gives the biggest variance (in terms of conductivity)
!         if (nch_short.le.nch_desire) then
!              nch_select=nch_short
!              do j=1,nch_select
!                ind_select_local(j)=j
 !             enddo
       
!         else
       
!              call finding_best_chan_list(nch_short,Rapp_short,nch_desire,ind_select_local) 
!              !ind_select_local: index of the selected channels assuming the full channel index 1,2,3,...,nch
!              nch_select=nch_desire

!         endif
        
!         ! preparing index and weighting for the inversion 
!         allocate (weight_forinv(nch_select),Rapp_forinv(nch_select))
!         allocate (ind_forinv(nch_select))
  
!         do j=1,nch_select
!            ind_forinv(j)=ind_short(ind_select_local(j))
!            weight_forinv(j)=weight_short(ind_select_local(j))
!            Rapp_forinv(j)=Rapp_short(ind_select_local(j))    
!         enddo
    
     
!       endif  
       
 
! Inversion 
 
        inv_switch_tmp=inv_switch
!        inv_switch_tmp(4)=0
      
!       hdo=dh_bhc2(i)
!       rto=rt_bhc2(i)
!       rmo=rm_bhc2(i)
!       epsro=epsr_bhc2(i)
       
 !      psm2i(1:5,i)=psm2_vec(1:5)
 !      adm2i(1:5,i)=adm2_vec(1:5)
 !      psm4i(1:5,i)=psm4_vec(1:5)
 !      adm4i(1:5,i)=adm4_vec(1:5)
      
!!          if(err_bhc(i).gt.misfit_epsr.and.inv_switch(4).eq.1) then
          
!!              hdi=dh_bhc(i)
!!              rmi=rm_bhc(i)
!!              rti=rt_bhc(i)
!!              ecci=ecc_bhc(i)
              hdi=dh_input(i)
              rmi=rm_input(i)
              rti=rt_input(i)
              epsri=epsr_input(i)

          
 !             inv_switch_tmp=inv_switch

 !             if (i_auto.eq.1) then   ! with selected channels

!                 call Ven_diinv_sub(dh_bhc(i),rm_bhc(i),rt_bhc(i),epsr_bhc(i),err_bhc(i),QC(i),
!	1             id_tool,icontig,i-1,nframe,bs,nch_select,ind_forinv, weight_forinv,psm2i,adm2i,psm4i,adm4i,
!	2             is_OBM,inv_switch,hdi,rmi,rti,epsri)
                 
!            call Ven_diinv_sub(dh_bhc2(i),rm_bhc2(i),rt_bhc2(i),epsr_bhc2(i),std_hd2(i),std_rm2(i),
!     1      std_rt2(i),std_epsr2(i),err_bhc2(i),dh_bhc4(i),rm_bhc4(i),rt_bhc4(i),epsr_bhc4(i),
!	2      std_hd4(i),std_rm4(i),std_rt4(i),std_epsr4(i),err_bhc4(i),id_tool,icontig,i-1,
!	3      nframe,bs,nch_select,nd_forinv, weight_forinv,psm2i,adm2i,psm4i,adm4i,
!     4       is_OBM,inv_switch,hdi,rmi,rti,epsri,mode)
	
!             else                    ! no selection of channels
              
!                 call Ven_diinv_sub(dh_bhc(i),rm_bhc(i),rt_bhc(i),epsr_bhc(i),err_bhc(i),QC(i),
!	1             id_tool,icontig,i-1,nframe,bs,nch,index_channel, weight_channel,psm2i,adm2i,psm4i,adm4i,
!	2             is_OBM,inv_switch,hdi,rmi,rti,epsri)
              
          call Ven_diinv_sub(dh_bhc2(i),rm_bhc2(i),rt_bhc2(i),epsr_bhc2(i),
     1      err_bhc2(i),QC2(i),dh_bhc4(i),rm_bhc4(i),rt_bhc4(i),epsr_bhc4(i),
	2      err_bhc4(i),QC4(i),id_tool,icontig,i-1,
	3      nframe,bs,nch,index_channel, weight_channel,psm2i,adm2i,psm4i,adm4i,
     4       is_OBM,inv_switch,hdi,rmi,rti,epsri,mode)        
                 
                 
!             endif	
             
             
!              call calculateQC_V2(dh_bhc2(i),rm_bhc2(i),rt_bhc2(i),epsr_bhc2(i),
!     &          err_bhc2(i),std_hd2(i),std_rm2(i),std_rt2(i),std_epsr2(i),QC2(i))

!            call calculateQC_V2(dh_bhc4(i),rm_bhc4(i),rt_bhc4(i),epsr_bhc4(i),
!     &          err_bhc4(i),std_hd4(i),std_rm4(i),std_rt4(i),std_epsr4(i),QC4(i))

!          end if    

!       endif

!!!       write(6,*)'after bhc,rt,dh,err='
!!!       write(6,*)rt_bhc(i),dh_bhc(i),rm_bhc(i),err_bhc(i)
!c  filename='wizard2_arc475.txt'  11   real inversion output

!       write(11,1001) depth(i),rt_bhc(i),dh_bhc(i),rm_bhc(i),epsr_bhc(i),err_bhc(i),
!     1                QC(i)
       
!c  filename='wizard1_arc475.txt'  12
!       write(12,1002) depth(i),rps2(1,i),rps2(2,i),rps2(3,i),rps2(4,i),
!     1      rps2(5,i),rad2(1,i),rad2(2,i),rad2(3,i),rad2(4,i),rad2(5,i),
!     1	              rt_input(i),rm_input(i)


      if(i_auto.eq.1) then
       deallocate(Rapp_short, weight_short)
       deallocate(ind_short)

       deallocate (weight_forinv,Rapp_forinv)
       deallocate (ind_forinv)
      endif
       
       
       enddo
           
      deallocate (Rapp_vec)
      deallocate (ind_select_local)   
      
!      close(11)
!      close(12)
!1000  format(a106)
!1001  format(7(1pe16.7))  
!1002  format(13f10.4)  
      return
      
      end
      
      
      
      subroutine finding_best_chan_list(nch,Rapp_vec,nch_desire,ind_select)
        implicit none
        ! Input
        integer nch
        real*4 Rapp_vec(nch) !input nch measurement channels in apparent resistivity
        integer nch_desire
        
        !output
        integer ind_select(nch_desire) !index of the selected channels, assuming index for Rapp_vec is 1,2,3,...,nch 
        
        !local
        real*4,allocatable::  cond_vec(:),ch_set_k(:)
        integer, allocatable:: ind_set_k(:)
        integer i,j,k,m,n,ii
        real*4 variance_max,var_tmp
        
        integer m2, h
        logical*1 mtc
               
        allocate(cond_vec(nch))
        allocate(ind_set_k(nch_desire))
        allocate(ch_set_k(nch_desire))

        do i=1,nch
           Cond_vec(i)=1/Rapp_vec(i)
        enddo

        !circulating for all possible nch_select combinations
        ind_select=0  
        variance_max=0.0

        n=nch
        k=nch_desire
        
        do i=1,k
          ind_set_k(i)=i   !the 1st k-combination out of 1:nch_short
          ch_set_k(i)=Cond_vec(ind_set_k(i))  
        enddo
        
                
        call comp_variance(k,ch_set_k,var_tmp)
        
        if (var_tmp.gt.variance_max) then
            variance_max=var_tmp
            ind_select(1:k)=ind_set_k(1:k)
        endif

        m2=ind_set_k(k)
        h=0
        mtc=1
                
        do while (ind_set_k(1).ne.(n-k+1))  !mtc: more to come, meaning haven't exhausted yet
            
            call nexksb(n,k,m2,h,ind_set_k,mtc)
            
            do i=1,k
              ch_set_k(i)=Cond_vec(ind_set_k(i))
            enddo
            
            
            call comp_variance(k,ch_set_k,var_tmp)

!            write(112,'(12I8, 20f20.8)') ind_set_k,m2,h,var_tmp,variance_max 
            
            if (var_tmp .gt. variance_max) then
                variance_max=var_tmp
                do i=1,k
                   ind_select(i)=ind_set_k(i)
                enddo
            endif

            
        enddo


        deallocate(cond_vec)
        deallocate(ind_set_k)
        deallocate(ch_set_k)
         
        return
      end

      subroutine comp_variance(n,x,var)
        implicit none
        !input
        integer n
        real*4 x(n)
        !output
        real*4 var
        !local
        integer i,j,k
        real*4 mean, tmp
        
        mean=0.0
        do i=1,n
          mean=mean+x(i)
        enddo
        
        mean=mean/n
        
        tmp=0.0
        do i=1,n
          tmp=tmp+(x(i)-mean)*(x(i)-mean)
        enddo
        
        var=sqrt(tmp/max(n-1,1))
        return
      end
      
      
       subroutine nexksb(n,k,m2,h,a,mtc)
        implicit none
        integer n,k !finding a(1:k) from 1:n
        integer m2,h,a(k)
        logical*1 mtc  
        
        integer i,j
        
        
        if (m2 .lt. n-h)  h=0
        
        h=h+1
        m2=a(k+1-h);
        do j=1,h
            a(k+j-h)=m2+j;
        enddo
        
        mtc=(a(1).ne.(n-k+1));
        
        return
        end