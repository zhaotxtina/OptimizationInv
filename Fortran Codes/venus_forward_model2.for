c subroutine arc_bhc_model(psmo,admo,psmbo,admbo,rmi,rti,hdi,ibadtx)
c
c To obtain 2MHz psm and adm data for ARC675 in a 2-layer TTM BHC model with
c parameters, rmi, rti, (in ohm-m) and hdi (in inch) 
c
c Inputs:	rmi- resistivity in mud layer (ohm-m), 0.015<rxoi<10
c		rti - resistivity in virgin layer (ohm-m), 0.1<rti<1000  
c		hdi - the hole diameter (inch) BS<hdi<18(A5),30(A6,8,9)
c		ibadtx - integer flag for contingency mixing (=1,2,3,4,5,6) 6 for use all 
c Outputs:psmo- 1x5 real array for mbhc PS (deg) at the 5 ARC675 spacing. 
c		admo- 1x5 real array for mbhc AD (dB) at the 5 ARC675 spacing.
c         psmbo- 1x5 real array for mbhc PS (deg) 400 kHz at the 5 ARC675 spacing
c         admbo- 1x5 real array for mbhc AD (dB) 400 kHz at the ARC675 spacing.
c Author:	Peter Wu/SPC/Nov 20,2001
c
	subroutine venus_forward_model(psmo,admo,psmbo,admbo,rmi,rti,hdi,epsri,ibadtx)
c
	parameter (r_min = 0.05, rm_min=0.01, rm_max=5000.0)
	include 'venus_ttmdivtab_def.inc'
	real*4 psmo(ntr),admo(ntr),psmbo(ntr),admbo(ntr),rmi,rti,hdi,epsri
	real*4 psmu(ntr),admu(ntr),psmub(ntr),admub(ntr)
	integer idx(2,4),ibadtx
	real*4 ps1(ntr,2,2,2,2),ps2(ntr,2,2,2),ps3(ntr,2,2),ps4(ntr,2) !index orders: ntr,hd,rt,rm,ecc
	real*4 ad1(ntr,2,2,2,2),ad2(ntr,2,2,2),ad3(ntr,2,2),ad4(ntr,2)
	real*4 psb1(ntr,2,2,2,2),psb2(ntr,2,2,2),psb3(ntr,2,2),psb4(ntr,2) !index orders: ntr,hd,rt,rm,ecc
	real*4 adb1(ntr,2,2,2,2),adb2(ntr,2,2,2),adb3(ntr,2,2),adb4(ntr,2)
	
	
	real*8 y1l,y2l,yl,xl,x1l,x2l,y1,y2,y,x1,x2,x
c	real*4, allocatable:: ecc_tmp(:)
	
	data inite / 0 /
	
	if(rti.lt.r_min) rti=r_min
	if(rmi.lt.rm_min) rmi=rm_min

c
	call find_index_v(iout,iflag,hdi,hd,nhd_t)
	if(iflag .eq. 0 ) then
	    idx(1,1)=iout
	    idx(2,1)=iout+1
	else
	    idx(1,1)=min(iout,iout+iflag)
	    idx(2,1)=max(iout,iout+iflag)
	end if
c
	call find_index_v(iout,iflag,rti,rt,nrt_t)
	if(iflag .eq. 0 ) then
	    idx(1,2)=iout
	    idx(2,2)=iout+1
	else
	    idx(1,2)=min(iout,iout+iflag)
	    idx(2,2)=max(iout,iout+iflag)
	end if
c
	call find_index_v(iout,iflag,rmi,rm,nrm_t)
	if(iflag .eq. 0 ) then
	    idx(1,3)=iout
	    idx(2,3)=iout+1
	else
	    idx(1,3)=min(iout,iout+iflag)
	    idx(2,3)=max(iout,iout+iflag)
	end if


c      ihd=min(idx(1,1),idx(2,1))
c      necc_tmp=necc_hd(ihd)
c      allocate(ecc_tmp(necc_tmp))
c      ecc_tmp(1:necc_tmp)=ecc(ihd,1:necc_tmp)

c	call find_index_v(iout,iflag,ecci,ecc_tmp,necc_tmp)
      call find_index_v(iout,iflag,epsri,epsr,nepsr_t)
	if(iflag .eq. 0 ) then
	    idx(1,4)=iout
	    idx(2,4)=iout+1
	else
	    idx(1,4)=min(iout,iout+iflag)
	    idx(2,4)=max(iout,iout+iflag)
	end if

c
	do i=1,2 !hd
		do j=1,2  !rt
			do k=1,2 !rm
				do l=1,ntr !spacing
				   do m=1,2 !epsr
					ps1(l,i,j,k,m)=ps_a(idx(j,2),idx(i,1),idx(k,3),l,idx(m,4))    ! ps_a(nrt,nhd,nrm,ntr,nepsr), ps1 order(tr,hd,rt,rm,epsr)
					ad1(l,i,j,k,m)=ad_a(idx(j,2),idx(i,1),idx(k,3),l,idx(m,4))  ! 2M is for ad_a, ps_a
					psb1(l,i,j,k,m)=ps_b(idx(j,2),idx(i,1),idx(k,3),l,idx(m,4))    ! ps_b(nrt,nhd,nrm,ntr,nepsr), psb1 order(tr,hd,rt,rm,epsr)
					adb1(l,i,j,k,m)=ad_b(idx(j,2),idx(i,1),idx(k,3),l,idx(m,4))   ! 400k is for ad_b, ps_b
                   end do
				end do
			end do
		end do
	end do
	
c interpolate the HD grid
	x1=hd(idx(1,1))
	x2=hd(idx(2,1))
	x=hdi
	if(hdi .lt. (hd(1)-1.0) ) x=hdi-1.0
	do j=1,2
		do k=1,2
			do l=1,ntr
			    do m=1,2
      		    	!2MHz
      		    	y1=ps1(l,1,j,k,m)
				    y2=ps1(l,2,j,k,m)
				    if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then   ! PS >0 positive, use log interpolation
					   y1l=dlog10(y1)
					   y2l=dlog10(y2)
					   yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
					   ps2(l,j,k,m)=10.0d0**(yl)
				    else                                  ! if any PS <0 , linear interpolation
					  y=(y1-y2)/(x2-x1)*(x2-x)+y2
					  ps2(l,j,k,m)=y
				    end if

                    y1=abs(ad1(l,1,j,k,m))
                    y2=abs(ad1(l,2,j,k,m))
                    y1l=dlog10(y1)
                    y2l=dlog10(y2)
                    yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
                    ad2(l,j,k,m)=10.0d0**(yl)              !  log interpolation, AD
                    
                    ! 400kHz
      		    	y1=psb1(l,1,j,k,m)
				    y2=psb1(l,2,j,k,m)
				    if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
					   y1l=dlog10(y1)
					   y2l=dlog10(y2)
					   yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
					   psb2(l,j,k,m)=10.0d0**(yl)
				    else
					  y=(y1-y2)/(x2-x1)*(x2-x)+y2
					  psb2(l,j,k,m)=y
				    end if

                    y1=abs(adb1(l,1,j,k,m))
                    y2=abs(adb1(l,2,j,k,m))
                    y1l=dlog10(y1)
                    y2l=dlog10(y2)
                    yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
                    adb2(l,j,k,m)=10.0d0**(yl)
	          end do
			end do
		end do
	end do
c interpolate the rt grid
	x1=(rt(idx(1,2)))
	x2=(rt(idx(2,2)))
	x1=dlog10(x1)
	x2=dlog10(x2)
	if(rti .gt. r_min ) then
		x=(rti)
		x=dlog10(x)
	else
		x=(r_min)
		x=dlog10(x)
	end if
	do k=1,2
		do l=1,ntr
		   do m=1,2
		     !2MHz
			 y1=ps2(l,1,k,m)
			 y2=ps2(l,2,k,m)
			 if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
		 		y1l=dlog10(y1)
		 		y2l=dlog10(y2)
				yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
				ps3(l,k,m)=10.0d0**(yl)
			 else
				y=(y1-y2)/(x2-x1)*(x2-x)+y2
				ps3(l,k,m)=y
			 end if
c
			 y1=abs(ad2(l,1,k,m))
			 y2=abs(ad2(l,2,k,m))
			 y1l=dlog10(y1)
			 y2l=dlog10(y2)
			 yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			 ad3(l,k,m)=10.0d0**(yl)
		     !400kHz
			 y1=psb2(l,1,k,m)
			 y2=psb2(l,2,k,m)
			 if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
		 		y1l=dlog10(y1)
		 		y2l=dlog10(y2)
				yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
				psb3(l,k,m)=10.0d0**(yl)
			 else
				y=(y1-y2)/(x2-x1)*(x2-x)+y2
				psb3(l,k,m)=y
			 end if
c
			 y1=abs(adb2(l,1,k,m))
			 y2=abs(adb2(l,2,k,m))
			 y1l=dlog10(y1)
			 y2l=dlog10(y2)
			 yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			 adb3(l,k,m)=10.0d0**(yl)
		   
		   
		   end do
		end do
	end do
	
!     interpolate the rm grid
      if(rmi.le.rm_max) then
	    x1=(rm(idx(1,3)))
	    x2=(rm(idx(2,3)))
	    x1=dlog10(x1)
	    x2=dlog10(x2)
	    if(rmi .gt. rm_min ) then
		    x=(rmi)
		    x=dlog10(x)
	    else
		    x=(rm_min)
		    x=dlog10(x)
	    end if
	    do l=1,ntr
	       do m=1,2
	         !2MHz
		       y1=ps3(l,1,m)
		       y2=ps3(l,2,m)
		       if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
			     y1l=dlog10(y1)
			     y2l=dlog10(y2)
			     yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			     ps4(l,m)=10.0d0**(yl)
		       else
			     y=(y1-y2)/(x2-x1)*(x2-x)+y2
			     ps4(l,m)=y
		       end if

		       y1=abs(ad3(l,1,m))
		       y2=abs(ad3(l,2,m))
		       y1l=dlog10(y1)
		       y2l=dlog10(y2)
		       yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
		       ad4(l,m)=10.0d0**(yl)

               !400kHz

		       y1=psb3(l,1,m)
		       y2=psb3(l,2,m)
		       if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
			     y1l=dlog10(y1)
			     y2l=dlog10(y2)
			     yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			     psb4(l,m)=10.0d0**(yl)
		       else
			     y=(y1-y2)/(x2-x1)*(x2-x)+y2
			     psb4(l,m)=y
		       end if

		       y1=abs(adb3(l,1,m))
		       y2=abs(adb3(l,2,m))
		       y1l=dlog10(y1)
		       y2l=dlog10(y2)
		       yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
		       adb4(l,m)=10.0d0**(yl)

	       end do
	    end do
    	
    	else

	    do l=1,ntr
	       do m=1,2
    	
        	    ad4(l,m)=ad3(l,2,m)
	          ps4(l,m)=ps3(l,2,m)
        	    adb4(l,m)=adb3(l,2,m)
	          psb4(l,m)=psb3(l,2,m)

	        
	        enddo
	    enddo
   	
    	endif
c

c interpolate the epsr grid
	x1=(epsr(idx(1,4)))
	x2=(epsr(idx(2,4)))
      x=epsri
	
	do l=1,ntr
	  !2MHz
		y1=ps4(l,1)
		y2=ps4(l,2)
		if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
			y1l=dlog10(y1)
			y2l=dlog10(y2)
			yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			psmu(l)=10.0d0**(yl)
		else
			y=(y1-y2)/(x2-x1)*(x2-x)+y2
			psmu(l)=y
		end if
c
		y1=ad4(l,1)
		y2=ad4(l,2)
		y1l=dlog10(y1)
		y2l=dlog10(y2)
		yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
		admu(l)=10.0d0**(yl)

	  !400kHz

		y1=psb4(l,1)
		y2=psb4(l,2)
		if ((y1 .gt. 0.0d0) .and. (y2 .gt. 0.0d0) ) then
			y1l=dlog10(y1)
			y2l=dlog10(y2)
			yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
			psmub(l)=10.0d0**(yl)
		else
			y=(y1-y2)/(x2-x1)*(x2-x)+y2
			psmub(l)=y
		end if
c
		y1=adb4(l,1)
		y2=adb4(l,2)
		y1l=dlog10(y1)
		y2l=dlog10(y2)
		yl=(y1l-y2l)/(x2-x1)*(x2-x)+y2l
		admub(l)=10.0d0**(yl)
		
		
	end do

c      deallocate(ecc_tmp)

	call mix_bhc_contingency(psmo,psmu,ibadtx)
	call mix_bhc_contingency(admo,admu,ibadtx)
	call mix_bhc_contingency(psmbo,psmub,ibadtx)
	call mix_bhc_contingency(admbo,admub,ibadtx)

c
	return
	end






c subroutine arc_bhc_model(psmo,admo,rmi,rti,hdi,ibadtx)
c
c To obtain 2MHz psm and adm data for ARC675 in a 2-layer TTM BHC model with
c parameters, rmi, rti, (in ohm-m) and hdi (in inch) 
c
c Inputs:	rmi- resistivity in mud layer (ohm-m), 0.015<rxoi<10
c		rti - resistivity in virgin layer (ohm-m), 0.1<rti<1000  
c		hdi - the hole diameter (inch) BS<hdi<18(A5),30(A6,8,9)
c		ibadtx - integer flag for contingency mixing (=1,2,3,4,5,6) 6 for use all 
c Outputs:psmo- 1x5 real array for mbhc PS (deg) at the 5 ARC675 spacing. 
c		admo- 1x5 real array for mbhc AD (dB) at the 5 ARC675 spacing.
c Author:	Gongli Wang/Keli Sun
c
	subroutine arc_bhc_model_1(psmo,admo,psmbo,admbo,rmi,rti,hdi,ecci,ibadtx)
c
	parameter (r_min = 0.05, rm_min=0.01,rm_max=5000.0)
	include 'arc_ttmbhctab_def.inc'
	real*4 psmo(ntr),admo(ntr),psmbo(ntr),admbo(ntr),rmi,rti,hdi,ecci
	real*4 psmu(ntr),admu(ntr),psmub(ntr),admub(ntr)
	integer idx(2,4),ibadtx
	real*4 ps1(ntr,2,2,2,2),ps2(ntr,2,2,2),ps3(ntr,2,2),ps4(ntr,2) !index orders: ntr,hd,rt,rm,ecc
	real*4 ad1(ntr,2,2,2,2),ad2(ntr,2,2,2),ad3(ntr,2,2),ad4(ntr,2)
	real*4 psb1(ntr,2,2,2,2),psb2(ntr,2,2,2),psb3(ntr,2,2),psb4(ntr,2) !index orders: ntr,hd,rt,rm,ecc
	real*4 adb1(ntr,2,2,2,2),adb2(ntr,2,2,2),adb3(ntr,2,2),adb4(ntr,2)


	real*8 y1l,y2l,yl,xl,x1l,x2l,y1,y2,y,x1,x2,x
	real*4, allocatable:: ecc_tmp(:)
	
	logical :: Uselinear
      INTEGER :: ISIGM, ISIGT, IRBH, IDECC, ID 
      INTEGER :: PSIGM, PSIGT, PRBH, PDECC
      INTEGER :: PSIGM1, PSIGM2, KSIGM, LSIGM
      INTEGER :: PSIGT1, PSIGT2, KSIGT, LSIGT
      INTEGER :: PRBH1, PRBH2, KRBH, LRBH
      INTEGER :: PDECC1, PDECC2, KDECC, LDECC, IDELSIG

      INTEGER, PARAMETER :: MLagInterp = 6
      REAL*8, DIMENSION(MLagInterp) :: SIGMGRID, SIGTGRID, RBHGRID, DECCGRID
      REAL*8, DIMENSION(MLagInterp) :: SIGASIGM, SIGASIGT, SIGARBH, SIGADECC
	REAL*8 :: DECC, RBH, SIGT, SIGM, SIGAInterp

	
	
	
	data inite / 0 /
	
	
	if(rti.lt.r_min) rti=r_min
	if(rmi.lt.rm_min) rmi=rm_min
c
	call find_index_v(iout,iflag,hdi,hd,nhd_t)
	if(iflag .eq. 0 ) then
	    idx(1,1)=iout
	    idx(2,1)=iout+1
	else
	    idx(1,1)=min(iout,iout+iflag)
	    idx(2,1)=max(iout,iout+iflag)
	end if
c
	call find_index_v(iout,iflag,rti,rt,nrt_t)
	if(iflag .eq. 0 ) then
	    idx(1,2)=iout
	    idx(2,2)=iout+1
	else
	    idx(1,2)=min(iout,iout+iflag)
	    idx(2,2)=max(iout,iout+iflag)
	end if
c
	call find_index_v(iout,iflag,rmi,rm,nrm_t)
	if(iflag .eq. 0 ) then
	    idx(1,3)=iout
	    idx(2,3)=iout+1
	else
	    idx(1,3)=min(iout,iout+iflag)
	    idx(2,3)=max(iout,iout+iflag)
	end if

      UseLinear=0
      PSIGM=min(idx(1,3),idx(2,3)) !index of node left of Rm
      NSIGM=nrm_t-1   !left out the last rm (which is for OBM)
      CALL FindUppLwrPnts(PSIGM, NSIGM, PSIGM1, PSIGM2, KSIGM, LSIGM, UseLinear)

      if(rmi.gt.rm_max) then
        PSIGM=nrm_t
        PSIGM1=nrm_t
        PSIGM2=nrm_t
        KSIGM=1
      endif

      PSIGT=min(idx(1,2),idx(2,2))
      NSIGT=nrt_t
      CALL FindUppLwrPnts(PSIGT, NSIGT, PSIGT1, PSIGT2, KSIGT, LSIGT, UseLinear)

      PRBH=min(idx(1,1),idx(2,1))
      NRBH=nhd_t
      CALL FindUppLwrPnts(PRBH, NRBH, PRBH1, PRBH2, KRBH, LRBH, UseLinear)
    !  
      SIGMGRID(1 : KSIGM) = rm(PSIGM1 : PSIGM2)
      SIGTGRID(1 : KSIGT) = rt(PSIGT1 : PSIGT2)
      RBHGRID(1 : KRBH) = hd(PRBH1 : PRBH2)
      DECC=ecci
      RBH=hdi
      SIGT=rti
      SIGM=rmi


      
      Do itr=1,ntr
         do i_adps=1,4 !1 for ad 2MHz, 2 for ps 2MHz, 3 for ad 400kHz, 4 for ps 400kHz
          !
          ! Lagrange interpolation
            DO ISIGM = PSIGM1, PSIGM2
               DO ISIGT = PSIGT1, PSIGT2
                  DO IRBH = PRBH1, PRBH2

                        necc_tmp=necc_hd(IRBH)
                        allocate(ecc_tmp(necc_tmp))
                        ecc_tmp(1:necc_tmp)=ecc(IRBH,1:necc_tmp)

                        call find_index_v(iout,iflag,ecci,ecc_tmp,necc_tmp)
                        if(iflag .eq. 0 ) then
                            idx(1,4)=iout
                            idx(2,4)=iout+1
                        else
                            idx(1,4)=min(iout,iout+iflag)
                            idx(2,4)=max(iout,iout+iflag)
                        end if
                        
                        PDECC=min(idx(1,4),idx(2,4))
                        NDECC=necc_tmp
                        
                        CALL FindUppLwrPnts(PDECC, NDECC, PDECC1, PDECC2, KDECC, LDECC, UseLinear)
                        DECCGRID(1 : KDECC) = ecc_tmp(PDECC1 : PDECC2)
                        DO IDECC = PDECC1, PDECC2

                           if (i_adps.eq.1) SIGADECC(IDECC - PDECC1 + 1)=ad_a(ISIGT,IRBH,ISIGM,itr,IDECC)
                           if (i_adps.eq.2) SIGADECC(IDECC - PDECC1 + 1)=ps_a(ISIGT,IRBH,ISIGM,itr,IDECC)
                           if (i_adps.eq.3) SIGADECC(IDECC - PDECC1 + 1)=ad_b(ISIGT,IRBH,ISIGM,itr,IDECC)
                           if (i_adps.eq.4) SIGADECC(IDECC - PDECC1 + 1)=ps_b(ISIGT,IRBH,ISIGM,itr,IDECC)

                           
      !                     SIGADECC(IDECC - PDECC1 + 1) = SIGATBL(ITOOL, IFREQ, IRBH, IDECC, ISIGM, ISIGT)
                        END DO
                        
                        CALL LagInterpol(DECCGRID(1 : KDECC), SIGADECC(1 : KDECC), DECC, KDECC, SIGARBH(IRBH - PRBH1 + 1))
                        
                        deallocate(ecc_tmp)
                  END DO
                  
                  CALL LagInterpol(RBHGRID(1 : KRBH), SIGARBH(1 : KRBH), RBH, KRBH, SIGASIGT(ISIGT - PSIGT1 + 1))          
               END DO
                  
               CALL LagInterpol(DLOG(SIGTGRID(1 : KSIGT)), SIGASIGT(1 : KSIGT), DLOG(SIGT), KSIGT,SIGASIGM(ISIGM - PSIGM1+1))
            END DO
            
            if (rmi.le.rm(nrm_t-1)) then  ! WBM
                 CALL LagInterpol(DLOG(SIGMGRID(1 : KSIGM)), SIGASIGM(1 : KSIGM), DLOG(SIGM), KSIGM, SIGAInterp)
            else
                 SIGAInterp=SIGASIGM(KSIGM)  !OBM
            endif
          
          
            if(i_adps.eq.1)  admu(itr)= SIGAInterp
            if(i_adps.eq.2)  psmu(itr)= SIGAInterp
            if(i_adps.eq.3)  admub(itr)= SIGAInterp
            if(i_adps.eq.4)  psmub(itr)= SIGAInterp

           
          enddo !i_adps
      enddo !itr


!      deallocate(ecc_tmp)

	call mix_bhc_contingency(psmo,psmu,ibadtx)
	call mix_bhc_contingency(admo,admu,ibadtx)
	call mix_bhc_contingency(psmbo,psmub,ibadtx)
	call mix_bhc_contingency(admbo,admub,ibadtx)

c
	return
	end








       SUBROUTINE FindUppLwrPnts(PSIGM, NSIGM, PSIGM1, PSIGM2, KSIGM, LSIGM, UseLinear)
!------------------------------------------------------------------------------           
! Subroutine to find grid points around the interpolaiton point
!    
! Author:   Gong Li Wang
! Time:     8/26/2013 
!------------------------------------------------------------------------------           
       IMPLICIT NONE
       LOGICAL, INTENT(IN) :: UseLinear ! linear interpolation (i.e. 2nd order Langrange interpolation
       INTEGER, INTENT(IN) :: PSIGM, NSIGM !PSIGM: the index of the grid position left (lower than) of the actual point.  NSIGM: length of the grid
       INTEGER, INTENT(OUT) :: PSIGM1, PSIGM2, KSIGM, LSIGM !nodes (PSIGM1,...,PSIGM2) are used for interpolation, KSIGM=PSIGM2-PSIGM1+1 is the length of nodes (or the order of the Langrage interpolation)
!
       IF(UseLinear) THEN
       PSIGM1 = PSIGM
       PSIGM2 = PSIGM + 1
       KSIGM = 2
       LSIGM = 1 ! no use
      ELSE   
         PSIGM1 = PSIGM - 1
         IF(PSIGM == 1) PSIGM1 = PSIGM  
         PSIGM2 = PSIGM + 2
         IF(PSIGM == NSIGM - 1) PSIGM2 = NSIGM  
         
         PSIGM1=max(PSIGM1,1)
         PSIGM2=min(PSIGM2,NSIGM)
!
!------------
! 4 points for points close to edge
!        PSIGM1 = PSIGM - 1
!        PSIGM2 = PSIGM + 2
!        IF(PSIGM == 1) then
!           PSIGM1 = PSIGM  
!           PSIGM2 = PSIGM + 3
!        end if   
!        IF(PSIGM == NSIGM - 1) then
!           PSIGM1 = NSIGM - 3  
!           PSIGM2 = NSIGM  
!        end if   
!------------
!
         KSIGM = PSIGM2 - PSIGM1 + 1
         LSIGM = PSIGM - PSIGM1 + 1 ! used by CH interpolation
      END IF   
!       
      END SUBROUTINE FindUppLwrPnts



      SUBROUTINE LagInterpol(X, SIGAGRID, XI, N, SIGA)
      !------------------------------------------------------------------------------           
      ! Subroutine to find function value at XI with Lagrange interpolation
      !    
      ! Author:   Gong Li Wang
      ! Time:     8/26/2013 
      !------------------------------------------------------------------------------           
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N !number of teh node points, i.e. the order of Langrange interpolation
        REAL*8, DIMENSION(N), INTENT(IN) :: X, SIGAGRID ! X(N), SIGAGRID(N), node position and signal
        REAL*8, INTENT(IN) :: XI !position to be interpolated
        REAL*8, INTENT(OUT) :: SIGA !function value at XI
        REAL*8, DIMENSION(N) :: L
        INTEGER :: J, M
      !
        IF(N > 4) STOP 'Order of Lagrange Interpolation > 4 ...'
        SIGA = 0.0D0
        L = 1.0D0
        DO J = 1, N
           DO M = 1, N
           IF(M /= J) L(J) = L(J) * (XI - X(M)) / (X(J) - X(M))
           END DO
           SIGA = SIGA + L(J) * SIGAGRID(J)
        END DO   
      !
      END SUBROUTINE LagInterpol

