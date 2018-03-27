!      subroutine weights(w1_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,nps,nad,
!     1   wps_factor)

	subroutine weights(w1_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,nps,nad,
     &      nch,index_channel,weight_channel,wps_factor)


!calculates the weight matrix applied to the transmitters in the inversion codes

      parameter(ntr=5,ssmall=1.e-8)
      integer i,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),index_all(20),nps,nad
      real*4 t1_wm,t2_wm,t3_wm,t4_wm,t5_wm,wps_factor,w_sum
      real*4 w1_in(4*ntr),tx_wm(ntr),ch(ntr),tx_wm_norm(ntr)
      integer nch
      integer index_channel(nch)
      real*4  weight_channel(nch)


!filling the array w_in with weights normalized to 1 
c   for split dual freq weight assignment
c      do i=1,ntr
c         w1_in(i)=0.3
c         w1_in(i+ntr)=1.0
c        ! w1_in(i+ntr)=wps_factor    ! set the AT weights  small
c         w1_in(i+2*ntr)=0.7
c        ! w1_in(i+3*ntr)=wps_factor   ! set the AT weights  small
c         w1_in(i+3*ntr)= 1.0
c      enddo
 
c  the following was used before Nov 10, 2015, then I changed the weights to make AT bigger
c       do i=1,ntr
c         w1_in(i)= wps_factor
c         w1_in(i+ntr)=1.0
c        ! w1_in(i+ntr)=wps_factor    ! set the AT weights  small
c         w1_in(i+2*ntr)=wps_factor
c        ! w1_in(i+3*ntr)=wps_factor   ! set the AT weights  small
c         w1_in(i+3*ntr)= 1.0
c      enddo

c     The following small section was copied from Rappin_ns_0929dll, so weights are the same as non-split freq
        do i=1,ntr
         w1_in(i)= 1.0
         w1_in(i+ntr)=wps_factor
        ! w1_in(i+ntr)=wps_factor    ! set the AT weights  small
         w1_in(i+2*ntr)= 1.0
        ! w1_in(i+3*ntr)=wps_factor   ! set the AT weights  small
         w1_in(i+3*ntr)= wps_factor
       enddo

       idx_ps2=0
       idx_ad2=0
       idx_ps4=0
       idx_ad4=0
      
       do j=1,nch
          if (index_channel(j).le.5) idx_ps2(index_channel(j))=index_channel(j)
          if (index_channel(j).gt.5.and.index_channel(j).le.10) idx_ad2(index_channel(j)-5) = index_channel(j)
          if (index_channel(j).gt.10.and.index_channel(j).le.15) idx_ps4(index_channel(j)-10)=index_channel(j)
          if (index_channel(j).gt.15) idx_ad4(index_channel(j)-15) = index_channel(j)
       enddo

       index_all(1:5)=idx_ps2(1:5)
       index_all(6:10)=idx_ad2(1:5)
       index_all(11:15)=idx_ps4(1:5)
       index_all(16:20)=idx_ad4(1:5)


       nps=0
       do j=1,5
         if(idx_ps2(j).ne.0) nps=nps+1
         if(idx_ps4(j).ne.0) nps=nps+1
       enddo
         
       nad=0
       do j=1,5
         if(idx_ad2(j).ne.0) nad=nad+1
         if(idx_ad4(j).ne.0) nad=nad+1
       enddo
      
       
       do j=1,nch
          w1_in(index_channel(j))= w1_in(index_channel(j))*weight_channel(j)
       enddo
 
!       Set weighting to be zero for the missing channels
       do j=1,20
         if(index_all(j).eq.0) w1_in(j)=0.0
       enddo

      return
      end
      



!      subroutine weights_old(w1_in,idx_ps2,idx_ad2,idx_ps4,idx_ad4,nps,nad,
!     1   t1_wm,t2_wm,t3_wm,t4_wm,t5_wm,wps_factor)
!
!
!!calculates the weight matrix applied to the transmitters in the inversion codes
!
!
!!	t1_wm	- weight multiplier for transmitter 1 (t1_wm>0)
!!	t2_wm	- weight multiplier for transmitter 2 (t2_wm>0)
!!	t3_wm	- weight multiplier for transmitter 3 (t3_wm>0)
!!	t4_wm	- weight multiplier for transmitter 4 (t4_wm>0)
!!	t5_wm	- weight multiplier for transmitter 5 (t5_wm>0)
!
!      parameter(ntr=5,ssmall=1.e-8)
!      integer i,idx_ps2(ntr),idx_ad2(ntr),idx_ps4(ntr),idx_ad4(ntr),
!     1  nps,nad
!      real*4 t1_wm,t2_wm,t3_wm,t4_wm,t5_wm,wps_factor,w_sum
!      real*4 w1_in(4*ntr),tx_wm(ntr),ch(ntr),tx_wm_norm(ntr)
!
!      tx_wm(1)=t1_wm
!      tx_wm(2)=t2_wm
!      tx_wm(3)=t3_wm
!      tx_wm(4)=t4_wm
!      tx_wm(5)=t5_wm
!
!      do i = 1,ntr
!      if (tx_wm(i).lt.ssmall)tx_wm(i)=ssmall 
!      end do
!!do the mixing by channel
!
!c     ch(1)=0.75*t1_wm+0.5*t2_wm-0.25*t3_wm
!c     ch(2)=0.25*t1_wm+0.5*t2_wm+0.25*t3_wm
!c     ch(3)=0.25*t2_wm+0.5*t3_wm+0.25*t4_wm
!c     ch(4)=0.25*t3_wm+0.5*t4_wm+0.25*t5_wm
!c     ch(5)=0.5*t4_wm+0.75*t5_wm-0.25*t3_wm
!c     ch(1)=0.75*t1_wm+0.5*t2_wm-0.25*t3_wm
!c     ch(2)=0.25*t1_wm+0.5*t2_wm+0.25*t3_wm
!c     ch(3)=0.25*t2_wm+0.5*t3_wm+0.25*t4_wm
!c     ch(4)=0.25*t3_wm+0.5*t4_wm+0.25*t5_wm
!c     ch(5)=0.5*t4_wm+0.75*t5_wm-0.25*t3_wm
!c Qiming May 20th, 2002 
!cNo mixing by channel. The weights are now strictly interpreted as
!c the weight for the BHC curves. Bad Transmitter has its own
!c switch icontig so no need to use weight to controll the transmitter weight 
!      do i =1,ntr
!      ch(i)=tx_wm(i)
!      end do 
!cc normalize such that the total weight for the ntr tran is ntr
!      w_sum = 0.0 
!      do i = 1,ntr
!      w_sum = w_sum + ch(i)
!      end do
!      do i=1,ntr
!      tx_wm_norm(i)=ntr*ch(i)/w_sum
!      end do
!
!!filling the array w_in with weights normalized to 1 
!!such as the sum of the elements of w_in<=20
!      do i=1,ntr
!      w1_in(i)=tx_wm_norm(i)*wps_factor
!      w1_in(i+ntr)=tx_wm_norm(i)
!!ccc make the weight of the 400 khz measurement almost zero
!!      w1_in(i+2*ntr)=ssmall
!!      w1_in(i+3*ntr)=ssmall
!ccc make the weight of the 400 khz measurement almost zero
!      w1_in(i+2*ntr)=tx_wm_norm(i)*wps_factor
!      w1_in(i+3*ntr)=tx_wm_norm(i)
!
!
!
!      idx_ps2(i)=i
!      idx_ad2(i)=i+ntr
!      idx_ps4(i)=i+2*ntr
!      idx_ad4(i)=i+3*ntr   !index orders: ps2, ad2, ps4, ad4
!      end do
!
!      nps=ntr
!      nad=ntr
!      return
!      end
