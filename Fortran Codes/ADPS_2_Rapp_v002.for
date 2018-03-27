c	subroutine ADPS_2_Rapp(rps2,rad2,rps4,rad4,
c	1	mode_t,i_init,abst,psm2,adm2,psm4,adm4,ibadtx)
c
c Purpose:To do transform from  PS,AD to Rapp resistivity   
c		  PS, AD 2MHz and 400KHz for ARC5,6,8,9 
c		one depth frame at a time.
c Inputs: mode_t	- integer switch for tool :
c				mode_t=3 for 3.00 in. ARC3 
c                 mode_t=4 for 4.75 in. Venus5 
c				mode_t=5 for 4.75 in. ARC5 
c				mode_t=6 for 6.75 in. ARC6 
c                 mode_t=7 for 6.75 in. Venus6 
c				mode_t=8 for 8.25 in. ARC8 
c				mode_t=9 for 9.00 in. ARC9 
c               mode_t=16 for EcoScope
c		i_init	- integer switch, =0 for the first time call to this 
c			  subroutine to initialize the transform table    
c			  switch >0 for subsequent call for a different depth frame.
c		abst - 1x1 real constant for absent value
coutputs: rps2- 1xntrm array(ntrm) for the MBHC 2MHz apparent Rps (ohm-m)
c		rad2- 1xntrm array(ntrm) for the MBHC 2MHz apparent Rad (ohm-m)
c		rps4- 1xntrm array(ntrm) for the MBHC 400KHz apparent Rps (ohm-m)
c		rad4- 1xntrm array(ntrm) for the MBHC 400KHz apparent Rad (ohm-m)
c		ibadtx- integer flag =1,2,3,4,5,6 for NG tx #, 6=all good 
c
c Inputs: psm2- 1xntrm array for input ntrm channels of ps (ntrm) at 2MHz
c		adm2- 1xntrm array for input ntrm channels of ad (ntrm) at 2MHz
c		psm4- 1xntrm array for input ntrm channels of ps (ntrm) at 400KHz
c		adm4- 1xntrm array for input ntrm channels of ad (ntrm) at 400KHz
c		abst value will be outputed for the corresponding inputs which
c		have abst value 
c Author:	Peter T. Wu/SPC		Feb 3. 2000
c		revision on 4/25/2000 for adding ARC3 transform table
c		revision on 11/01/00 to load the new A3 table which has different 
c		array dimension and 412.04KHz instead of the standard array and freq.
c		revision on 10/31/01 to add contingency mixing
c       revision on 01/03/11 to add EcoScope
c       for test, tzhao  Feb 10, 2015
      subroutine ADPS_2_Rapp(rps2,rad2,rps4,rad4,
	1	mode_t,i_init,abst,psm2,adm2,psm4,adm4,ibadtx)
	 ! transfer raw signal (db/deg) to apparent resistivity 
c
	parameter (ntr=5,ntrm=5,ntf2=29,ntf4=31)
c
	integer mode_t,i1(1),i2(1),i,i_init,j,il1(1),il2(1)
      real*4  rt2(ntf2), rt4(ntf4), xu2(ntf2, 10), xu4(ntf4,10)
	real*4 rps2(ntrm),rad2(ntrm),rps4(ntrm),rad4(ntrm)
	real*4 ps2(ntr),ad2(ntr),psa2(ntr),ada2(ntr)
	real*4 ps4(ntr),ad4(ntr),psa4(ntr),ada4(ntr)
	real*4 psm2(ntrm),adm2(ntrm),psm4(ntrm),adm4(ntrm)
	real*4 ru2(ntf2),psat2(ntr),adat2(ntr)
	real*4 psu2(ntf2,ntr),adu2(ntf2,ntr)
	real*4 pu2(ntf2,ntr*2),p_airt2(2*ntr)
	real*4 psum2(ntf2,ntrm),adum2(ntf2,ntrm)
	real*4 ru4(ntf4),psat4(ntr),adat4(ntr)
	real*4 psu4(ntf4,ntr),adu4(ntf4,ntr)
	real*4 pu4(ntf4,ntr*2),p_airt4(2*ntr)
	real*4 psum4(ntf4,ntrm),adum4(ntf4,ntrm)
	real*4 ps2_offset(ntr),ps4_offset(ntr)
	real*4 abst
      
      data (rt2(i),i=1,29)/ 0.10,0.15,0.20,0.30,0.50,0.70,1.0,1.5,2.0,
     &  3.0,5.0,7.0,10.0,15.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,
     &  300.0,500.0,700.0,1000.0,1500.0, 2000.0, 3000.0, 5000.0/
      
      data (rt4(i),i=1,31)/ 0.05, 0.07, 0.10,0.15,0.20,0.30,0.50,0.70,
     &  1.0,1.5,2.0,3.0,5.0,7.0,10.0,15.0,20.0,30.0,50.0,70.0,100.0,
     & 150.0,200.0,300.0,500.0,700.0,1000.0,1500.0,2000.0,3000.0,5000.0/
c
	include 'a3a_n_tf_table.inc'	
	include 'a3b_n_tf_table.inc'	
	include 'a475a_tf_table_n.inc'	
	include 'a475b_tf_table_n.inc'	
      include 'v475a_tf_table_n.inc'	
	include 'v475b_tf_table_n.inc'
c      include 'mxw475a_tf_table_n.inc'	
c	include 'mxw475b_tf_table_n.inc'	
	include 'a65a_tf_table_n.inc'	
	include 'a65b_tf_table_n.inc'
      include 'v65a_tf_table_n.inc'	
	include 'v65b_tf_table_n.inc'
	include 'a8a_tf_table.inc'	
	include 'a8b_tf_table.inc'	
	include 'a9a_tf_table.inc'	
	include 'a9b_tf_table.inc'
	include 'ea_tf_table.inc'	
	include 'eb_tf_table.inc'	
	save
c
c load the right table according to mode_t 
c
	if ( i_init .eq. 0 ) then
	  if ( mode_t .eq. 3) then
c load 2MHz transform table for ARC3
	    do i=1,ntf2
	      ru2(i)=ru3a(i)
	      do j=1,ntr
			psu2(i,j)=psu3a(i,j)
			adu2(i,j)=adu3a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400 KHz transform table for ARC3
	    do i=1,ntf4
	      ru4(i)=ru3b(i)
	      do j=1,ntr
			psu4(i,j)=psu3b(i,j)
			adu4(i,j)=adu3b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
          
	  else if ( mode_t .eq. 4) then
c load 2MHz transform table for ARC5
	    do i=1,ntf2
	      ru2(i)=ru475a(i)
	      do j=1,ntr
			psu2(i,j)=psu475a(i,j)
			adu2(i,j)=adu475a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400 KHz transform table for ARC5
	    do i=1,ntf4
	      ru4(i)=ru475b(i)
	      do j=1,ntr
			psu4(i,j)=psu475b(i,j)
			adu4(i,j)=adu475b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
          
          else if ( mode_t .eq. 5) then
c load 2MHz transform table for venus5
	    do i=1,ntf2
	      ru2(i)=ruv475a(i)
	      do j=1,ntr
			psu2(i,j)=psuv475a(i,j)
			adu2(i,j)=aduv475a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400 KHz transform table for venus5
	    do i=1,ntf4
	      ru4(i)=ruv475b(i)
	      do j=1,ntr
			psu4(i,j)=psuv475b(i,j)
			adu4(i,j)=aduv475b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)

       else if ( mode_t .eq. 6) then
c load 2MHz transform table for venus6
	    do i=1,ntf2
	      ru2(i)=ruv65a(i)
	      do j=1,ntr
			psu2(i,j)=psuv65a(i,j)
			adu2(i,j)=aduv65a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400KHz transform table for venus6
	    do i=1,ntf4
	      ru4(i)=ruv65b(i)
	      do j=1,ntr
			psu4(i,j)=psuv65b(i,j)
			adu4(i,j)=aduv65b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
          
       else if ( mode_t .eq. 7) then
c load 2MHz transform table for ARC6
	    do i=1,ntf2
	      ru2(i)=ru65a(i)
	      do j=1,ntr
			psu2(i,j)=psu65a(i,j)
			adu2(i,j)=adu65a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400KHz transform table for ARC6
	    do i=1,ntf4
	      ru4(i)=ru65b(i)
	      do j=1,ntr
			psu4(i,j)=psu65b(i,j)
			adu4(i,j)=adu65b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
          
	  else if ( mode_t .eq. 8) then
c load 2MHz transform table for ARC8
	    do i=1,ntf2
	      ru2(i)=ru8a(i)
	      do j=1,ntr
			psu2(i,j)=psu8a(i,j)
			adu2(i,j)=adu8a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400KHz transform table for ARC8
	    do i=1,ntf4
	      ru4(i)=ru8b(i)
	      do j=1,ntr
			psu4(i,j)=psu8b(i,j)
			adu4(i,j)=adu8b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
	  else if ( mode_t .eq. 9) then
c load 2MHz transform table for ARC9
	    do i=1,ntf2
	      ru2(i)=ru9a(i)
	      do j=1,ntr
			psu2(i,j)=psu9a(i,j)
			adu2(i,j)=adu9a(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400KHz transform table for ARC9
	    do i=1,ntf4
	      ru4(i)=ru9b(i)
	      do j=1,ntr
			psu4(i,j)=psu9b(i,j)
			adu4(i,j)=adu9b(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
		else if ( mode_t .eq. 16) then
c load 2MHz transform table for EcoScope
	    do i=1,ntf2
	      ru2(i)=rua(i)
	      do j=1,ntr
			psu2(i,j)=psua(i,j)
			adu2(i,j)=adua(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum2,psu2,ntf2,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum2,adu2,ntf2,ntr,2,ibadtx)
c load 400KHz transform table for EcoScope
	    do i=1,ntf4
	      ru4(i)=rub(i)
	      do j=1,ntr
			psu4(i,j)=psub(i,j)
			adu4(i,j)=adub(i,j)
	      end do
	    end do
		call mix_bhc_2d_contingency(psum4,psu4,ntf4,ntr,2,ibadtx)
		call mix_bhc_2d_contingency(adum4,adu4,ntf4,ntr,2,ibadtx)
	  else
	    stop 'Wrong tool mode, mode_t=3,4,5,6,7,8,9,16 only'
	  end if
      end if
      
c      !  assign them to xu2 and xu4, let's see if I need to use the MBHC table
c      do i=1,ntf2
c       do j=1,5
c          xu2(i,j)= psu2(i,j)
c          xu2(i,j+5) = adu2(i,j)
c       end do
c      end do
      
c      do i=1,ntf4
c       do j=1,5
c          xu4(i,j)= psu4(i,j)
c          xu4(i,j+5) = adu4(i,j)
c       end do
c      end do
       
      
c
c   transform from raw signal PS/AD to apparent resistivity Rapp
c
      do i=1,5
          if (abs(psm2(i)-abst).lt.1.0e-5) then	!absent in absent out
		rps2(i)=abst
	  else	
         call linter_log_v(rps2(i),il1,il2,psum2(:,i),ru2,psm2(i),0,ntf2,1)
        end if
          
        if (abs(adm2(i)-abst).lt.1.0e-5) then	!absent in absent out
		rad2(i)=abst
	   else	
        call linter_log_v(rad2(i),il1,il2,adum2(:,i),ru2,adm2(i),0,ntf2,1)
         end if
        
          if (abs(psm4(i)-abst).lt.1.0e-5) then	!absent in absent out
		rps4(i)=abst
	  else	
        call linter_log_v(rps4(i),il1,il2,psum4(:,i),ru4,psm4(i),0,ntf4,1)
        end if
        
         if (abs(adm4(i)-abst).lt.1.0e-5) then	!absent in absent out
		rad4(i)=abst
	  else	
        call linter_log_v(rad4(i),il1,il2,adum4(:,i),ru4,adm4(i),0,ntf4,1)
        end if
        
      end do
      
      
c	do i=1,ntrm
c!	  if (rps2(i) .eq. abst) then	!absent in absent out
c	  if (abs(rps2(i)-abst).lt.1.0e-5) then	!absent in absent out
c		psm2(i)=abst
c	  else							!normal inverse transform
c	    call linter_log_v(psm2(i),i1,i2,ru2,psum2(1,i),rps2(i),
c	1	0,ntf2,1)
c	  end if
c
c	  if (abs(rad2(i)-abst).lt.1.0e-5) then	!absent in absent out
c		adm2(i)=abst
c	  else							!normal inverse transform
c	    call linter_log_v(adm2(i),i1,i2,ru2,adum2(1,i),rad2(i),
c	1	0,ntf2,1)
c	  end if
c
c	  if (abs(rps4(i)-abst).lt.1.0e-5) then	!absent in absent out
c		psm4(i)=abst
c	  else							!normal inverse transform
c	    call linter_log_v(psm4(i),i1,i2,ru4,psum4(1,i),rps4(i),
c	1	0,ntf4,1)
c	  end if
c
c	  if (abs(rad4(i)-abst).lt.1.0e-5) then	!absent in absent out
c		adm4(i)=abst
c	  else							!normal inverse transform
c	    call linter_log_v(adm4(i),i1,i2,ru4,adum4(1,i),rad4(i),
c	1	0,ntf4,1)
c	  end if
c
c	end do
c
	return
	end
