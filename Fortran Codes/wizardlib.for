c
c subroutine arc5_mduanis_th2(modet,ps,ad,kh,beta)
c Purpose:	to compute the output phase shift, attenuation and phase ave. 
c		due to a point magnetic dipole at ARC spacing and frequency.
c		The formation is uniform TI anisotropic medium without borehole.
c		No constitutive relation between sigma and epslom is used
c		and the computation is based on Kh and Beta
c Inputs:	modet - integer switch for tool, =0,1,2,3,4 for arc3,5,6,8,9
c		kh	-complex horizontal wave number (MKS unit)
c		beta	-complex anisotropy factor (MKS unit)
c Outputs:	ps	-phase shift (deg)
c		ad	-attenuation (dB)
c Note:	all arguments are double precision
c Author:	Peter T. Wu/Anadrill, Sept. 1996
c-----------------------------------------------------------------
	subroutine arc5_mduanis_th2(modet,ps,ad,kh,beta)
c
	parameter(ntr_r=5,nr=6)
	integer modet
	complex*16 kh,beta,ikl,ai,v(nr)
	real*8 ps(ntr_r),ad(ntr_r)
	real*8 phase(nr),p_uw(nr),amp_db(nr)
	real*8 tr,xold,flip,tr_min
        data pi/ 3.1415927/
c
	ai=dcmplx(0.0d0,1.0d0)
	tr_min=13.0d0
	if(modet .le. 1 ) tr_min=7.0d0
	do i=1,6
	  iflag=i-1
	  tr=(tr_min+float(i-1)*6.0d0)*0.0254d0
	  ikl=ai*kh*tr
	  v(i)=-ai*((-2.0d0+ikl)*exp(ikl)+ikl*exp(ikl*beta))
	1	/(tr*tr*tr)
	  amp_db(i)=20.d0*log10(abs(v(i)))
	  phase(i)=(180.0/pi)*datan2(dimag(v(i)),dreal(v(i)))
	  call unwrap1d(p_uw(i),xold,flip,phase(i),
	1	300.0d0,-20.0d0,360.0d0,iflag)
c	  write(6,10) i,v(i),amp_db(i),phase(i),p_uw(i)
c10	  format(1x,'i,v,amp,p:',i2,5f12.7)
	end do
c
	do i=1,5
	  ps(i)=(p_uw(i+1)-p_uw(i))
	  ad(i)=amp_db(i)-amp_db(i+1)
c	  pa(i)=(p_uw(i+1)+p_uw(i))*0.5d0-90.0d0
	end do
	return
	end
c
c subroutine unwrap1d(xout,xold,flip,xin,th_u,th_d,step,iflag)
c
c Purpose: Double precision version of unwrap phase one point at a time. 
c For a long sequences of phase data, this routine need to be called consequetively.
c Inputs:	xin	- input phase angle (deg)
c		th_u	- upward step threshold (deg)
c		th_d	- downward step threshold (deg)
c		step	- wraping size (deg)
c		iflag	- =0 for starting
c			  =1 for continuing
c		xold	- old phase angle (It is input for iflag=1)
c		flip	- the cumulated amount of flip in phase angle 
c			  for the sequence (input for iflag=1)
c Output:	xout	- unwrapped phase (deg)
c		xold	- old phase angle 
c		flip	- the cumulated amount of flip in phase angle 
c			  for the sequence 
c Author:	Peter Wu/Anadrill/July 26, 1996
c-------------------------------------------------------------
c
	subroutine unwrap1d(xout,xold,flip,xin,th_u,th_d,step,iflag)
c
	real*8 xout,xold,flip,xin,th_u,th_d,step
c
	if (iflag .eq. 0) then
	  xout=xin
	  xold=xout
	  flip = 0.0d0
	else
	  xnew=xin+flip
	  if( (xnew-xold) .le. th_d ) then
	    xout=xnew + step 
	    flip =flip + step
	    xold=xout
	  else if ( (xnew-xold) .ge. th_u ) then
	    xout=xnew - step 
	    flip = flip - step
	    xold=xout
	  else
	    xout=xnew
	    xold=xout
	  end if
	end if
	return
	end
c
c
c
	subroutine mpdpanis2(modet,kh_in,beta_in,psm_o,adm_o,ibadtx)
c
      parameter (ntr_r=5)
	parameter (ntr=5)
	integer modet,ibadtx
	complex*16 kh_in,beta_in
	real*4 psm_o(ntr),adm_o(ntr)
	real*4 ps_s(ntr_r),ad_s(ntr_r)
	real*8 ps(ntr_r),ad(ntr_r)
c
	call arc5_mduanis_th2(modet,ps,ad,kh_in,beta_in)
c
	do i=1,ntr_r
	  ps_s(i)=ps(i)
	  ad_s(i)=ad(i)
	end do
	call mix_bhc_contingency(psm_o,ps_s,ibadtx)
	call mix_bhc_contingency(adm_o,ad_s,ibadtx)
c
	return
	end
c
c
	subroutine five_to_two(rh_in,rv_in,ang_in,epsh2_in,epsv2_in,
	1	epsh4_in,epsv4_in,modef,kh2,beta2,kh4,beta4)
c
	complex*16 kh2,beta2,csigh2,csigv2
	complex*16 kh4,beta4,csigh4,csigv4
	real*8 sdip,cdip,sdip2,cdip2
	real*4 rh_in,rv_in,ang_in,epsh2_in,epsv2_in,epsh4_in,epsv4_in
	real*4 epsh2,epsv2,epsh4,epsv4
	integer modef
c
      real*8 c0,pi,eps0,f2,f4,fmhz2,fmhz4
	real*8 omega2,omegamu2,omegaeps2,omega4,omegamu4,omegaeps4
      complex*16 un,ai,ciomegamu2,ciomegaeps2,ciomegamu4,ciomegaeps4
c
	common /constant/  c0,pi,eps0,un,ai,f2,fmhz2,f4,fmhz4,
	1	omega2,omegamu2,omegaeps2,omega4,omegamu4,omegaeps4,
     2	ciomegamu2,ciomegaeps2,ciomegamu4,ciomegaeps4
c
	sdip=dsin(ang_in*pi/180.0)
	sdip2=sdip*sdip
	cdip=dcos(ang_in*pi/180.0)
	cdip2=cdip*cdip
	if(modef .eq. 1 ) then
		csigh2=dcmplx(1.d0/rh_in,-omegaeps2*epsh2_in)
		csigv2=dcmplx(1.d0/rv_in,-omegaeps2*epsv2_in)
		kh2=cdsqrt(ai*omegamu2*csigh2)
		beta2=cdsqrt(cdip2+csigv2/csigh2*sdip2)
	else if( modef .eq. 2) then
		csigh4=dcmplx(1.d0/rh_in,-omegaeps4*epsh4_in)
		csigv4=dcmplx(1.d0/rv_in,-omegaeps4*epsv4_in)
		kh4=cdsqrt(ai*omegamu4*csigh4)
		beta4=cdsqrt(cdip2+csigv4/csigh4*sdip2)
	else if (modef .eq. 3) then
		csigh2=dcmplx(1.d0/rh_in,-omegaeps2*epsh2_in)
		csigv2=dcmplx(1.d0/rv_in,-omegaeps2*epsv2_in)
		kh2=cdsqrt(ai*omegamu2*csigh2)
		beta2=cdsqrt(cdip2+csigv2/csigh2*sdip2)
		csigh4=dcmplx(1.d0/rh_in,-omegaeps4*epsh4_in)
		csigv4=dcmplx(1.d0/rv_in,-omegaeps4*epsv4_in)
		kh4=cdsqrt(ai*omegamu4*csigh4)
		beta4=cdsqrt(cdip2+csigv4/csigh4*sdip2)
	else if(modef .eq. 4 ) then
		epsh2=108.5*rh_in**(-0.35)+5.0
		epsv2=108.5*rv_in**(-0.35)+5.0
		csigh2=dcmplx(1.d0/rh_in,-omegaeps2*epsh2)
		csigv2=dcmplx(1.d0/rv_in,-omegaeps2*epsv2)
		kh2=cdsqrt(ai*omegamu2*csigh2)
		beta2=cdsqrt(cdip2+csigv2/csigh2*sdip2)
	else if( modef .eq. 5) then
		epsh4=280.0*rh_in**(-0.46)+5.0
		epsv4=280.0*rv_in**(-0.46)+5.0
		csigh4=dcmplx(1.d0/rh_in,-omegaeps4*epsh4)
		csigv4=dcmplx(1.d0/rv_in,-omegaeps4*epsv4)
		kh4=cdsqrt(ai*omegamu4*csigh4)
		beta4=cdsqrt(cdip2+csigv4/csigh4*sdip2)
	else if (modef .eq. 6) then
		epsh2=108.5*rh_in**(-0.35)+5.0
		epsv2=108.5*rv_in**(-0.35)+5.0
		epsh4=280.0*rh_in**(-0.46)+5.0
		epsv4=280.0*rv_in**(-0.46)+5.0
		csigh2=dcmplx(1.d0/rh_in,-omegaeps2*epsh2)
		csigv2=dcmplx(1.d0/rv_in,-omegaeps2*epsv2)
		kh2=cdsqrt(ai*omegamu2*csigh2)
		beta2=cdsqrt(cdip2+csigv2/csigh2*sdip2)
		csigh4=dcmplx(1.d0/rh_in,-omegaeps4*epsh4)
		csigv4=dcmplx(1.d0/rv_in,-omegaeps4*epsv4)
		kh4=cdsqrt(ai*omegamu4*csigh4)
		beta4=cdsqrt(cdip2+csigv4/csigh4*sdip2)
	end if
c
	return
	end
c
c
	subroutine set_constant
c
      real*8 c0,pi,eps0,f2,f4,fmhz2,fmhz4
	real*8 omega2,omegamu2,omegaeps2,omega4,omegamu4,omegaeps4
      complex*16 un,ai,ciomegamu2,ciomegaeps2,ciomegamu4,ciomegaeps4
c
	common /constant/  c0,pi,eps0,un,ai,f2,fmhz2,f4,fmhz4,
	1	omega2,omegamu2,omegaeps2,omega4,omegamu4,omegaeps4,
     2	ciomegamu2,ciomegaeps2,ciomegamu4,ciomegaeps4
c
      c0=2.99792458d+08
      un=dcmplx(1.0,0.0)
      ai=dcmplx(0.0,1.0)
      PI=4.0d0*DATAN(1.0d0)
      eps0=1.0/(c0*c0*4.0d-07*pi)
      fmhz2=2.0d0
      f2=fmhz2*1.0d6
      omega2=2.0d0*pi*f2
	omegaeps2=omega2*eps0
      omegamu2=omega2*4.0d0*pi*1.0d-7
	ciomegamu2=ai*omegamu2
	ciomegaeps2=ai*omega2*eps0
      fmhz4=0.4d0
      f4=fmhz4*1.0d6
      omega4=2.0d0*pi*f4
	omegaeps4=omega4*eps0
      omegamu4=omega4*4.0d0*pi*1.0d-7
	ciomegamu4=ai*omegamu4
	ciomegaeps4=ai*omega4*eps0
	return
	end
c
c
c	subroutine mix_bhc_2d_contingency(xm,x,n1,n2,id,ibadtx)
c Purpose:	do mix bhc on input 2D matrix x dimensioned by x(n1,n2)
c		id =1 or 2 for the dimension that represents the TR spacing
c		the id_th dimension is 5 for 5 TR spacings of ARC.
c		The output will be xm(n1,n2) 
c Author: P. Wu/SPC Oct. 25, 2001
c
	subroutine mix_bhc_2d_contingency(xm,x,n1,n2,id,ibadtx)
c
	parameter (ntr = 5)
	real*4 x(n1,n2),xm(n1,n2),tmp(ntr),tmpm(ntr)
	integer ibadtx,id
c
	if(id .eq. 1) then
	  do i=1,n2
		call mix_bhc_contingency(xm(1,i),x(1,i),ibadtx)
	  end do
	else if (id .eq. 2) then
	  do i=1,n1
		do j=1,ntr
			tmp(j)=x(i,j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do j=1,ntr
			xm(i,j)=tmpm(j)
		end do
	  end do
	end if
c
	return
	end
c subroutine mix_bhc_3d_contingency(xm,x,n1,n2,n3,id,ibadtx)
c Purpose:	do mix bhc on input 3D matrix x dimensioned by x(n1,n2,n3)
c		id =1,2, or 3 for the dimension that represents the TR spacing
c		the id_th dimension is 5 for 5 TR spacings of ARC.
c		The output will be xm(n1,n2,n3)
c Author: P. Wu/SPC Oct. 25, 2001
c
	subroutine mix_bhc_3d_contingency(xm,x,n1,n2,n3,id,ibadtx)
c
	parameter (ntr = 5)
	integer ibadtx,id
	real*4 x(n1,n2,n3),xm(n1,n2,n3),tmp(ntr),tmpm(ntr)
c
	if(id .eq. 1) then
	do j=1,n3
	  do i=1,n2
		call mix_bhc_contingency(xm(1,i,j),x(1,i,j),ibadtx)
	  end do
	end do
c
	else if(id .eq. 2 ) then
	do j=1,n3
	  do i=1,n1
		do k=1,ntr
		  tmp(k)=x(i,k,j)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do k=1,ntr
		  xm(i,k,j)=tmpm(k)
		end do
	  end do
	end do
c
	else if(id .eq. 3 ) then
	do j=1,n2
	  do i=1,n1
		do k=1,ntr
		  tmp(k)=x(i,j,k)
		end do
		call mix_bhc_contingency(tmpm,tmp,ibadtx)
		do k=1,ntr
		  xm(i,j,k)=tmpm(k)
		end do
	  end do
	end do
	end if
	return
	end
c
c subroutine mix_bhc_contingency(xm,x,nbad)
c Purpose:	do mix bhc on input array x dimensioned by 5.with
c		contingency mixing algorithm where the excluded bad tx is denoted 
c		by nbad, nbed=1,2,...,5.  For normal mixing, nbed=6
c		The output will be xm(5), here we drop the two interpolated
c		spacings where there is no physical sensors
c Author: P. Wu/Anadrill/June 13, 2001
c
	subroutine mix_bhc_contingency(xm,x,nbad)
c
	real*4 x(5),xm(5)
	real*4 b(5,5,6)
	data  (b(i,1,6),i=1,5) / 0.75, 0.5 , -0.25,   0.0,  0.0  /
	data  (b(i,2,6),i=1,5) / 0.25, 0.5 ,  0.25,   0.0,  0.0  /
	data  (b(i,3,6),i=1,5) / 0.0 , 0.25,   0.5,  0.25,  0.0  /
	data  (b(i,4,6),i=1,5) / 0.0 , 0.0 ,  0.25,   0.5,  0.25 /
	data  (b(i,5,6),i=1,5) / 0.0 , 0.0 , -0.25,   0.5,  0.75 /
c
	data  (b(i,1,1),i=1,5) / 0.00, 0.75,  1.00, -0.25, -0.50 /
	data  (b(i,2,1),i=1,5) / 0.00, 0.5 ,  0.75,   0.0, -0.25 /
	data  (b(i,3,1),i=1,5) / 0.0 , 0.25,   0.5,  0.25,  0.0  /
	data  (b(i,4,1),i=1,5) / 0.0 , 0.0 ,  0.25,   0.5,  0.25 /
	data  (b(i,5,1),i=1,5) / 0.0 , 0.0 , -0.25,   0.5,  0.75 /
c
	data  (b(i,1,2),i=1,5) /0.875, 0.00, -0.25,   0.5, -0.125/
	data  (b(i,2,2),i=1,5) / 0.50, 0.00,  0.25,   0.5, -0.25 /
	data  (b(i,3,2),i=1,5) /0.125, 0.00,   0.5,   0.5, -0.125/
	data  (b(i,4,2),i=1,5) / 0.0 , 0.0 ,  0.25,   0.5,  0.25 /
	data  (b(i,5,2),i=1,5) / 0.0 , 0.0 , -0.25,   0.5,  0.75 /
c
	data  (b(i,1,3),i=1,5) /0.625, 0.5 ,  0.00,   0.0, -0.125/
	data  (b(i,2,3),i=1,5) /0.375, 0.5 ,  0.00,   0.0,  0.125/
	data  (b(i,3,3),i=1,5) / 0.25, 0.25,  0.00,  0.25,  0.25 /
	data  (b(i,4,3),i=1,5) /0.125, 0.0 ,  0.00,   0.5,  0.375/
	data  (b(i,5,3),i=1,5) /-0.125, 0.0,  0.00,   0.5,  0.625/
c
	data  (b(i,1,4),i=1,5) / 0.75, 0.5 , -0.25,   0.0,  0.0  /
	data  (b(i,2,4),i=1,5) / 0.25, 0.5 ,  0.25,   0.0,  0.0  /
	data  (b(i,3,4),i=1,5) /-0.125,0.5 ,  0.5 ,   0.0,  0.125/
	data  (b(i,4,4),i=1,5) /-0.25, 0.5 ,  0.25,   0.0,  0.50 /
	data  (b(i,5,4),i=1,5) /-0.25, 0.5 , -0.25,   0.0,  1.00 /
c
	data  (b(i,1,5),i=1,5) / 0.75, 0.5 , -0.25,   0.0,  0.0  /
	data  (b(i,2,5),i=1,5) / 0.25, 0.5 ,  0.25,   0.0,  0.0  /
	data  (b(i,3,5),i=1,5) / 0.0 , 0.25,   0.5,  0.25,  0.0  /
	data  (b(i,4,5),i=1,5) /-0.125,0.0 , 0.625,   0.5,  0.0  /
	data  (b(i,5,5),i=1,5) /-0.375,0.0 , 0.875,   0.5,  0.0  /
c
cccc Qiming 03/31/02 when invalid selction,default to normal selection
        if(nbad.lt.1.or.nbad.gt.6) nbad= 6
ccccccccc 
	do i=1,5
	  xm(i)=0.0
	  do j=1,5
		xm(i)=xm(i)+x(j)*b(j,i,nbad)
c               write(6,*)'i,j,x(j),xm(i)=',i,j,x(j),xm(i)
	  end do
	end do
	return
	end
c
c subroutine linter_v(yi,i1,i2,x,y,xi,iextra,n,ni)
c To obtain point-vise linearly interpolated/extrapolated value
c yi from a 1D table x vs y when a value xi is given, x needs
c to be monotonically increasing or decresaing.
c Inputs:	x,y - vectors of x- and y-axis values forming a 1D table
c		xi  - a vector of x-axis values for which the corresponding y-axis 
c			value are to be interpolated
c		iextra - switch for extrapolation, =1 for yes, =0 for no
c		n - dimensin of the x or y array
c		ni - dimension of the xi array
c Outputs:	yi  - vector for linearly interpolated/extrapolated y-axis value 
c			corresponding to x1
c		i1,i2 - vectors for indices of x such that 
c			x(i1(j))<= xi(j) <=x(i2(j)) for interp., j=1,...,ni
c			i1(j)=1,i2(j)=2, if xi(j)<=x(k) for all k (for extrapolation)
c			i1(j)=n-1,i2(j)=n, if xi(j)>=x(k) for all k(for extrapolation)
c Author:	Peter Wu/Anadrill/Feb. 24, 1995
	subroutine linter_v(yi,i1,i2,x,y,xi,iextra,n,ni)
	parameter (big=1.0e10)
	real*4 yi(1),x(1),y(1),xi(1)
	integer i1(1),i2(1)
	real*8 y1,y2,dxi,x1,x2,d1,d2
	do i=1,ni
	  call find_index_v(ix,iflag,xi(i),x,n)
	  if (iflag .eq. -1) then
		i1(i)=ix-1
		i2(i)=ix
	  else
                i1(i)=ix
                i2(i)=ix+1
	  endif
	  y1=y(i1(i))
          y2=y(i2(i))
          x1=x(i1(i))
          x2=x(i2(i))
          dxi=xi(i)
	  if(iextra .eq. 0) then
	    if(iflag .ne. 0) then
		yi(i)=y(ix)
	    else
	      if( x2 .ne. x1) then
		d2=y1+(y2-y1)*(dxi-x1)/(x2-x1)
		yi(i)=d2
	      else
		tmp=(y2-y1)*(dxi-x1)
		yi(i)=sign(big,tmp)
	      end if
	    endif
	  else
	    if( x2 .ne. x1) then
	      d2=y1+(y2-y1)*(dxi-x1)/(x2-x1)
	      yi(i)=d2
	    else
	      tmp=(y2-y1)*(dxi-x1)
	      yi(i)=sign(big,tmp)
	    end if
	  endif
	enddo
	return
	end
c
c subroutine linter_log_v(yi,i1,i2,x,y,xi,iextra,n,ni)
c To obtain point-vise linearly interpolated/extrapolated value
c yi in log domain from a 1D table x vs y when a value xi is given, x needs
c to be monotonically increasing or decresaing.
c Inputs:	x,y - vectors of x- and y-axis values forming a 1D table
c		xi  - a vector of x-axis values for which the corresponding y-axis 
c			value are to be interpolated
c		iextra - switch for extrapolation, =1 for yes, =0 for no
c		n - dimensin of the x or y array
c		ni - dimension of the xi array
c Outputs:	yi  - vector for linearly interpolated/extrapolated y-axis value 
c			corresponding to x1
c		i1,i2 - vectors for indices of x such that 
c			x(i1(j))<= xi(j) <=x(i2(j)) for interp., j=1,...,ni
c			i1(j)=1,i2(j)=2, if xi(j)<=x(k) for all k (for extrapolation)
c			i1(j)=n-1,i2(j)=n, if xi(j)>=x(k) for all k(for extrapolation)
c Author:	Peter Wu/Anadrill/Feb. 24, 1995
	subroutine linter_log_v(yi,i1,i2,x,y,xi,iextra,n,ni)
	parameter( exp_max=10.0, big=1.0e10)
	real*4 yi(ni),x(n),y(n),xi(ni)
	integer i2(ni),i1(ni)
	real*8 y1,y2,dxi,x1,x2,y1l,y2l,dxil,x1l,x2l,d1,d2
	do i=1,ni
	  call find_index_v(ix,iflag,xi(i),x,n)
	  if (iflag .eq. -1) then
		i1(i)=ix-1
		i2(i)=ix
	  else
                i1(i)=ix
                i2(i)=ix+1
	  endif
	  y1=y(i1(i))
	  y2=y(i2(i))
	  x1=x(i1(i))
	  x2=x(i2(i))
	  dxi=xi(i)
	  if(iextra .eq. 0) then
	    if(iflag .ne. 0) then
		yi(i)=y(ix)
	    else
	      if((x1.gt.0).and.(x2.gt.0).and.(dxi.gt.0).and.
	1	(y1.gt.0).and.(y2.gt.0) ) then
		y1l=dlog10(y1)
		y2l=dlog10(y2)
		xil=dlog10(dxi)
		x1l=dlog10(x1)
		x2l=dlog10(x2)
		if (x2l .ne. x1l) then
		  d1=y1l+(y2l-y1l)*(xil-x1l)/(x2l-x1l)
		  if(d1.gt.exp_max) then
			d1=exp_max
		  elseif(d1.lt.-exp_max) then
			d1=-exp_max
		  endif
		  d2=10**d1
		else
		  tmp=(y2l-y1l)*(xil-x1l)
		  d2=sign(big,tmp)
		end if
	      else
		if (x2 .ne. x1) then
		  d2=y1+(y2-y1)*(dxi-x1)/(x2-x1)
		else
		  tmp=(y2-y1)*(dxi-x1)
		  d2=sign(big,tmp)
		end if
	      end if	
	      yi(i)=d2
	    endif
	  else
	      if((x1.gt.0).and.(x2.gt.0).and.(dxi.gt.0).and.
	1	(y1.gt.0).and.(y2.gt.0) ) then
		y1l=dlog10(y1)
		y2l=dlog10(y2)
		xil=dlog10(dxi)
		x1l=dlog10(x1)
		x2l=dlog10(x2)
		if (x2l .ne. x1l) then
		  d1=y1l+(y2l-y1l)*(xil-x1l)/(x2l-x1l)
		  if(d1.gt.exp_max) then
			d1=exp_max
		  elseif(d1.lt.-exp_max) then
			d1=-exp_max
		  endif
		  d2=10**d1
		else
		  tmp=(y2l-y1l)*(xil-x1l)
                  d2=sign(big,tmp)
                end if
	    else
		if (x2 .ne. x1) then
		  d2=y1+(y2-y1)*(dxi-x1)/(x2-x1)
		else
		  tmp=(y2-y1)*(dxi-x1)
                  d2=sign(big,tmp)
                end if
	    end if	
	    yi(i)=d2
	  endif
	enddo
	return
	end
c subroutine find_index_v(iout,iflag,xi,x,n) find the index iout such that 
c xi's value is in between x(iout) and x(iout+1).  The input n-element vector 
c x needs to be monotonic increasing or decreasing.
c If the xi is out side the range of x, iflag=+-1, otherwise 0
c For extrapolation, using the two points at indices : iout, iout+iflag
c i.e. i1=min(iout,iout+iflag), i2=max(iout,iout+iflag)
c Author:       Peter Wu/Anadrill/June 25, 1996
	subroutine find_index_v(iout,iflag,xi,x,n)
	real*4 x(1)
c
	call find_max(x,n,x_max,i_max)
	call find_min(x,n,x_min,i_min)
	if (xi .le. x_min) then
		iout=i_min
		if (i_max .gt. i_min) then
		  iflag=1.0
		else
		  iflag=-1.0
		endif
	elseif (xi .ge. x_max) then
		iout=i_max
                if (i_max .gt. i_min) then
                  iflag=-1.0
                else
                  iflag=1.0
                endif
 	else
	  iflag=0.0
	  ilow=min(i_min,i_max)
	  ihi=max(i_min,i_max)
	  imid=ilow+(ihi-ilow)/2
	  do while (imid .ne. ilow)
	    xlow=min(x(ilow),x(imid))
	    xmid=max(x(ilow),x(imid))
	    if ((xi .ge. xlow) .and. (xi .le. xmid) ) then
		ihi=imid
		imid=ilow+(ihi-ilow)/2
	    else
		ilow=imid
		imid=ilow+(ihi-ilow)/2
	    end if
	  end do
	  iout=ilow
	end if
	return
	end
c------------------------------------------------------------------------
c subroutine find_max(x,nx,x_max,i_max)
c find the max of a array and return the max value and index of the max.
c inputs:	x - input array dimensioned by nx
c		nx- length of the input array x
c outputs:	x_max - the max of x
c		i_max - the index of the max of x
c author: P. Wu/Anadrill/June, 1996
c
	subroutine find_max(x,nx,x_max,i_max)
c
	real*4 x(nx),x_max
	x_max=1.e-30	!enter a very small number for the application
	do i=1,nx
	  if(x(i) .gt. x_max) then	! > imply the the first occurance for multiple
		x_max=x(i)		! values of the same size number
		i_max=i
	  endif
	enddo
	return
	end
c subroutine find_min(x,nx,x_min,i_min)
c find the min of a array and return the min value and index of the min.
c inputs:       x - input array dimensioned by nx
c               nx- length of the input array x
c outputs:      x_min - the min of x
c               i_min - the index of the min of x
c author: P. Wu/Anadrill/June, 1996

        subroutine find_min(x,nx,x_min,i_min)
        real*4 x(nx),x_min
        x_min=1.e30    !enter a very large number for the application
        do i=1,nx
          if(x(i) .lt. x_min) then ! < imply the first occurance for multiple
                x_min=x(i)         ! values of the same size number
                i_min=i
          endif
        enddo
        return
        end
c
c
c
	subroutine l2norm(y_in,x_in,nx,errpc)
c
	parameter (max_m=20)
	real*4 y_in(max_m),x_in(max_m)
	real*4 errpc,sumsqr
c
c
	sumsqr=0.0
c
	do i=1,nx
		tmp=abs(y_in(i))
		if ((i .le. 5) ) then
			if (tmp .lt. 0.02 ) then
				tmp = 0.02
			end if
		else
			if (tmp .lt. 0.006) then
				tmp = 0.006
			end if
		end if
	    if(nx.GT.10) then 
	      if ((i .le. 15) ) then
			if (tmp .lt. 0.02 ) then
				tmp = 0.02
			end if
		  else
			if (tmp .lt. 0.006) then
				tmp = 0.006
			end if
		  end if
	    end if
	  sumsqr=sumsqr+
     1 x_in(i)*x_in(i)/tmp/tmp
	end do
c
	errpc=sqrt(sumsqr/float(nx))
c
	return
	end
c
	subroutine el2norm(y_in,w_in,m,idx_ps2,idx_ad2,idx_ps4,idx_ad4,cpsm2_m,cadm2_m,cpsm4_m,cadm4_m,el2_o)
c
	parameter (ntr=5,max_m=20)
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
	do i=1,max_m
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


c
	subroutine find_min_s(imin,xmin,x_in,nx)
c
	real*4 xmin,x_in(nx)
	integer nx,imin
	real*4, parameter:: abst=-999.25
c
	xmin=1.0e+30	!start with a very large number
	do i=1,nx
	 if(abs(x_in(i)-abst).gt.1.e-8) then
	    if( x_in(i) .lt. xmin ) then
	      xmin=x_in(i)
	      imin=i
	    end if
	 end if
	end do
	return
	end
c
c
        subroutine find_max_s(imax,xmax,x_in,nx)
c
        real*4 xmax,x_in(nx)
        integer nx,imax
        real*4, parameter:: abst=-999.25

c
        xmax=-1.0e+30     !start with a very small number
        do i=1,nx
          if(abs(x_in(i)-abst).gt.1.e-8) then
            
            if(x_in(i).gt.xmax) then
               xmax=x_in(i)
               imax=i
            end if
          end if
        end do
        return
        end
c
C subroutine find_sort(Xout, Xin, NRow, NCol) - Do bubble sort
C
C PURPOSE:
C	This subroutine sorts the input array to ascending order using a bubble sort.
C INPUTS:
C		Xin	Input data array to be sorted
C		Nrow	Number of rows in the input data array
C		Ncol	Number of columns in the input data array
C OUTPUTS:
C		Xout	Sorted data array
C
C Note:	This routine sorts only the first column of the input data array.
C	All of the other columns simply tag along with the first one.
C******************************* SWS-E-Sonic *****************************
	SUBROUTINE find_sort (Xout, Xin, Nrow, Ncol)
C
	real*4 Xin(Nrow,Ncol),Xout(Nrow,Ncol)
	LOGICAL*1 Done
c
c isolate output from input
c
	do i=1,ncol
	  do j=1,nrow
		xout(j,i)=xin(j,i)
	  end do
	end do
C
C...Perform the bubble sort
C
	I = 1
	Done = .FALSE.
	DO WHILE (I .LE. Nrow .AND. .NOT. Done)
	  Done = .TRUE.
	  DO J=1,Nrow-I
	    IF (Xout(J,1) .gt. Xout(J+1,1)) THEN
	      DO K=1,Ncol
		tmp	= xout(j,k)
	        Xout(j,k)   = Xout(J+1,k)
	        Xout(J+1,k) = tmp
	      ENDDO
	      Done = .FALSE.
	    ENDIF
	  ENDDO
	  I = I + 1
	ENDDO
C
C...Done
C
	RETURN
	END
c
c subroutine find_mean(x,n,xm)
c purpose:	compute xm as the mean of the input array x with dimension n
c P. WU/Anadrill/June 27, 1996
c
	subroutine find_mean(x,n,xm)
c
	real*4 x(n)
c
	xm=0.0
	do i=1,n
	  xm=xm+x(i)
	end do
	xm=xm/float(n)
	return
	end
c subroutine find_range(x,n,xr)
c purpose:	compute xr as the range of the input array x with dimension n
c P. WU/Anadrill/June 27, 1996
c
	subroutine find_range(x,n,xr)
c
	real*4 x(n)
c
	call find_max(x,n,x_max,i_max)
	call find_min(x,n,x_min,i_min)
	xr=x_max-x_min
	return
	end
c
c subroutine filtert(fx,x,nx,f,nf)
c purpose:	To filter a sequence x with a fitler f via convolution. 
c			nf, the length of the filter, is an odd number.
c			nx, the length of the sequence, is assumed bigger than nf.
c			nf/2 constant values are padded at the front and back of the
c			fitlered sequence to make a smooth output record of length nx
	subroutine filtert(fx,x,nx,f,nf)
c
	real*4 x(nx),fx(nx),f(nf)
	integer nx,nf,nf2,nf2p1
c
	nf2=nf/2
	nf2p1=nf2+1
c
	do i=nf2p1, nx-nf2 
		fx(i)=0.0
		k=0
		do j=i-nf2,i+nf2
			k=k+1
			fx(i)=fx(i)+x(j)*f(k)
		end do
	end do
c
	do i=1,nf2	
		fx(i)=fx(nf2p1)
	end do
	do i=nx-nf2+1,nx
		fx(i)=fx(nx-nf2)
	end do
	return
	end
c
c subroutine medfil01(ixf,ix,nx,nf)
c purpose:	do median filter of length nf on input sequence ix with 0 and 1 
c			as elements.  The nf is an odd number.  nx, the length of the
c			input sequence, is assumed larger then nf.  nf/2 constant 
c			values are padded at the front and back of the fitlered sequence
c			to make a smooth output record of length nx
c
	subroutine medfil01(ixf,ix,nx,nf)
c
	integer ix(nx),ixf(nx)
	integer nx,nf
c
	nf2=nf/2
	nf2p1=nf2+1
c
	isum=0
	do i=1,nf
		isum=isum+ix(i)
	end do
	if ( isum .ge. nf2p1 ) then
		ixf(nf2p1)=1
	else
		ixf(nf2p1)=0
	end if
	do i=nf2p1+1, nx-nf2
		isum=isum+ix(i+nf2)-ix(i-nf2p1)
		if ( isum .ge. nf2p1 ) then
			ixf(i)=1
		else
			ixf(i)=0
		end if
	end do
	do i=1,nf2
		ixf(i)=ixf(nf2p1)
	end do
	do i=nx-nf2+1,nx
		ixf(i)=ixf(nx-nf2)
	end do
	return
	end
c
c
c subroutine medfil010(ixf,ix,nx,nf,nc,nmodel,model_compute)
c purpose:	do counting filter of length nf on input sequence ix with 0 and 1 
c			as elements.  The nf is an odd number.  nx, the length of the
c			input sequence, is assumed larger then nf.  nf/2 constant 
c			values are padded at the front and back of the fitlered sequence
c			to make a smooth output record of length nx
cc             selection is done for multiple curves nc 
c              the purpose is to give one dominant flag over the nf interval out of
c              nc channels. 
c
        subroutine medfil010(ixf,ix,nx,nf,nc,nmodel,model_compute)
c
	integer ix(nx,nc),ixf(nx,nc)
	integer nx,nf,nfh,nsum,jselect,j,nc
        integer nmodel,isummul,isum,is 
        integer  model_compute(nmodel+2) 
c
        if(nc.gt.nmodel+2)then
        write(6,*)'nevent more than nmodel+2, stop'
        return
        end if 
 
        do i=1,nx
        do j=1,nc
        ixf(i,j)=0
        end do 
        end do 

c   
	nf2=nf/2
	nf2p1=nf2+1
c
c first the initial summation
        do j=1,nc 
        isum=0	
	do i=1,nf
		isum=isum+ix(i,j)
	end do
cc use the last entry temperorily to save the sum for later use 
        ixf(nx,j)=isum
        end do
c check the multieffect accumulation as well 
cc for the multiple effect, all flag needs to be zero 
        isummul=0
	do i=1,nf
        is=0 
        do j=1,nc 
        is=is+ix(i,j)
        end do
        if (is.eq.0)isummul=isummul+1 
	end do
c
cc use the last entry temperorily to save the sum for later use 

cc
        nsum=0
        jselect=0
        do j=1,nc 
        if (ixf(nx,j).gt.nsum.and.model_compute(j).ne.0)then
        nsum=ixf(nx,j)
        jselect=j
        end if
        end do
cccc the multieffect acculation needs to be smaller for the flag to set 
        if(jselect.ne.0.and.nsum.ge.isummul)ixf(nf2p1,jselect)=1

ccc Now go on for other points
	do i=nf2p1+1,nx-nf2
        nsum=0
        jselect=0
c       write(6,*)'iflagN=',(ix(i-nf2p1,j),j=1,nc)
c       write(6,*)'iflagP=',(ix(i+nf2,j),j=1,nc)
c
        do j=1,nc
        ixf(nx,j)=ixf(nx,j)+ix(i+nf2,j)-ix(i-nf2p1,j)
c       write(6,*)'j,nsum,ixf(nx,j),model_compute(j)=',
c    1                  j,nsum,ixf(nx,j),model_compute(j)
        if (ixf(nx,j).gt.nsum.and.model_compute(j).ne.0)then
        nsum=ixf(nx,j)
        jselect=j
        end if
        end do
c       write(6,*)'ixf(nx)=',(ixf(nx,j),j=1,nc)

        is=0 
        do j=1,nc 
        is=is+ix(i-nf2p1,j)
        end do
        if (is.eq.0)isummul = isummul -1
        is=0 
        do j=1,nc 
        is=is+ix(i+nf2,j)
        end do
        if (is.eq.0)isummul = isummul +1
c
c       write(6,*)'nsum,isummul,jselect',nsum,
c    1     isummul,jselect
        if(jselect.ne.0.and.nsum.ge.isummul)ixf(i,jselect)=1
c
        end do
c
cc now the two ends 

	do i=1,nf2
                do j=1,nc
		ixf(i,j)=ixf(nf2p1,j)
                end do 
	end do
	do i=nx-nf2+1,nx
                do j=1,nc
		ixf(i,j)=ixf(nx-nf2,j)
                end do 
	end do

cc
	return
        end 




c Subroutine array_log10(xlog,x,n)
c to convert an array x to log10(x)
c
	subroutine array_log10(xlog,x,n)
c
	integer n
	real*4 xlog(n),x(n)
c
	do i=1,n
		xlog(i)=log10(x(i))
	end do
	return
	end
c
c Subroutine array_10power(x10,x,n)
c to convert an array x to 10**x
c
	subroutine array_10power(x10,x,n)
c
	integer n
	real*4 x10(n),x(n)
c
	do i=1,n
		x10(i)=10**(x(i))
	end do
	return
	end
c
c=====================================================================
c common code for LM algorithm for LS inversion -SNLS1
c
c
      SUBROUTINE CHKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
      INTEGER M,N,LDFJAC,MODE
      REAL X(N),FVEC(M),FJAC(LDFJAC,N),XP(N),FVECP(M),ERR(M)
      INTEGER I,J
      REAL EPS,EPSF,EPSLOG,EPSMCH,FACTOR,ONE,TEMP,ZERO
      REAL R1MACH
      DATA FACTOR,ONE,ZERO /1.0E2,1.0E0,0.0E0/
C
C***FIRST EXECUTABLE STATEMENT  CHKDER
      EPSMCH = R1MACH(4)
C
      EPS = SQRT(EPSMCH)
C
      IF (MODE .EQ. 2) GO TO 20
C
C        MODE = 1.
C
         DO 10 J = 1, N
            TEMP = EPS*ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = EPS
            XP(J) = X(J) + TEMP
   10       CONTINUE
         GO TO 70
   20 CONTINUE
C
C        MODE = 2.
C
         EPSF = FACTOR*EPSMCH
         EPSLOG = ALOG10(EPS)
         DO 30 I = 1, M
            ERR(I) = ZERO
   30       CONTINUE
         DO 50 J = 1, N
            TEMP = ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = ONE
            DO 40 I = 1, M
               ERR(I) = ERR(I) + TEMP*FJAC(I,J)
   40          CONTINUE
   50       CONTINUE
         DO 60 I = 1, M
            TEMP = ONE
            IF (FVEC(I) .NE. ZERO .AND. FVECP(I) .NE. ZERO
     1          .AND. ABS(FVECP(I)-FVEC(I)) .GE. EPSF*ABS(FVEC(I)))
     2         TEMP = EPS*ABS((FVECP(I)-FVEC(I))/EPS-ERR(I))
     3                /(ABS(FVEC(I)) + ABS(FVECP(I)))
            ERR(I) = ONE
            IF (TEMP .GT. EPSMCH .AND. TEMP .LT. EPS)
     1         ERR(I) = (ALOG10(TEMP) - EPSLOG)/EPSLOG
            IF (TEMP .GE. EPS) ERR(I) = ZERO
   60       CONTINUE
   70 CONTINUE
C
      RETURN
C
C     LAST CARD OF SUBROUTINE CHKDER.
C
      END
      REAL FUNCTION ENORM(N,X)
      INTEGER N
      REAL X(N)
      INTEGER I
      REAL AGIANT,FLOATN,ONE,RDWARF,RGIANT,S1,S2,S3,XABS,X1MAX,X3MAX,
     1     ZERO
      DATA ONE,ZERO,RDWARF,RGIANT /1.0E0,0.0E0,3.834E-20,1.304E19/
C***FIRST EXECUTABLE STATEMENT  ENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = ABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
C
C              SUM FOR LARGE COMPONENTS.
C
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
C
C              SUM FOR SMALL COMPONENTS.
C
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
C
C           SUM FOR INTERMEDIATE COMPONENTS.
C
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
C
C     CALCULATION OF NORM.
C
      IF (S1 .EQ. ZERO) GO TO 100
         ENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX)
     1         ENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX)
     1         ENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            ENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
C
C     LAST CARD OF FUNCTION ENORM.
C
      END
      SUBROUTINE FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
      INTEGER M,N,LDFJAC,IFLAG
      REAL EPSFCN
      REAL X(N),FVEC(M),FJAC(LDFJAC,N),WA(M)
      INTEGER I,J
      REAL EPS,EPSMCH,H,TEMP,ZERO
      REAL R1MACH
      DATA ZERO /0.0E0/
C***FIRST EXECUTABLE STATEMENT  FDJAC3
      EPSMCH = R1MACH(4)
C
      EPS = SQRT(AMAX1(EPSFCN,EPSMCH))
ccc Qiming changed it 
      eps=0.01
      eps=0.0001
C      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
C           ARE TO BE RETURNED BY FCN.
      IFLAG = 1
      DO 20 J = 1, N
         TEMP = X(J)
         H = EPS*ABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
         CALL FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC)
         IF (IFLAG .LT. 0) GO TO 30
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE FDJAC3.
C
      END
      
      SUBROUTINE FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      
      INTEGER FUNCTION I1MACH(I)
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /      5 /
C     DATA IMACH( 2) /      6 /
C     DATA IMACH( 3) /      7 /
C     DATA IMACH( 4) /      6 /
C     DATA IMACH( 5) /     48 /
C     DATA IMACH( 6) /      6 /
C     DATA IMACH( 7) /      2 /
C     DATA IMACH( 8) /     39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /      8 /
C     DATA IMACH(11) /     13 /
C     DATA IMACH(12) /    -50 /
C     DATA IMACH(13) /     76 /
C     DATA IMACH(14) /     26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) / 6LOUTPUT /
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /  11 /
C     DATA IMACH( 2) /  12 /
C     DATA IMACH( 3) /   8 /
C     DATA IMACH( 4) /  10 /
C     DATA IMACH( 5) /  16 /
C     DATA IMACH( 6) /   2 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    0 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) /    3 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   23 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   38 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    4 /
C     DATA IMACH( 4) /    1 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   23 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   39 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    4 /
C     DATA IMACH( 4) /    1 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   23 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   55 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     6 /
C     DATA IMACH( 4) /     7 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -126 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1015 /
C     DATA IMACH(16) /  1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     1 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    36 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    27 /
C     DATA IMACH(12) /  -128 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    60 /
C     DATA IMACH(15) / -1024 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
c	DATA IMACH( 1) /    5 /
c	DATA IMACH( 2) /    6 /
c	DATA IMACH( 3) /    5 /
c	DATA IMACH( 4) /    6 /
c	DATA IMACH( 5) /   32 /
c	DATA IMACH( 6) /    4 /
c	DATA IMACH( 7) /    2 /
c	DATA IMACH( 8) /   31 /
c	DATA IMACH( 9) / 2147483647 /
c	DATA IMACH(10) /    2 /
c	DATA IMACH(11) /   24 /
c	DATA IMACH(12) / -127 /
c	DATA IMACH(13) /  127 /
c	DATA IMACH(14) /   56 /
c	DATA IMACH(15) / -127 /
c	DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /     5 /
C     DATA IMACH( 2) /     6 /
C     DATA IMACH( 3) /     6 /
C     DATA IMACH( 4) /     6 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -126 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1022 /
C     DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /    1 /
C     DATA IMACH( 2) /    1 /
C     DATA IMACH( 3) /    0 /
C     DATA IMACH( 4) /    1 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE IBM PC - MICROSOFT FORTRAN
C
      DATA IMACH( 1) /     5 /
      DATA IMACH( 2) /     6 /
      DATA IMACH( 3) /     6 /
      DATA IMACH( 4) /     0 /
      DATA IMACH( 5) /    32 /
      DATA IMACH( 6) /     4 /
      DATA IMACH( 7) /     2 /
      DATA IMACH( 8) /    31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /     2 /
      DATA IMACH(11) /    24 /
      DATA IMACH(12) /  -126 /
      DATA IMACH(13) /   127 /
      DATA IMACH(14) /    53 /
      DATA IMACH(15) / -1022 /
      DATA IMACH(16) /  1023 /
C
C     MACHINE CONSTANTS FOR THE IBM PC - PROFESSIONAL FORTRAN
C     AND LAHEY FORTRAN
C
C     DATA IMACH( 1) /     4 /
C     DATA IMACH( 2) /     7 /
C     DATA IMACH( 3) /     7 /
C     DATA IMACH( 4) /     0 /
C     DATA IMACH( 5) /    32 /
C     DATA IMACH( 6) /     4 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    24 /
C     DATA IMACH(12) /  -126 /
C     DATA IMACH(13) /   127 /
C     DATA IMACH(14) /    53 /
C     DATA IMACH(15) / -1022 /
C     DATA IMACH(16) /  1023 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH=IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE(OUTPUT,9000)
9000  FORMAT('1ERROR    1 IN I1MACH - I OUT OF BOUNDS ')
C
C     CALL FDUMP
C
C
      STOP
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,WA1,WA2)
      INTEGER N,LDR
      INTEGER IPVT(N)
      REAL DELTA,PAR
      REAL R(LDR,N),DIAG(N),QTB(N),X(N),SIGMA(N),WA1(N),WA2(N)
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      REAL DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001,SUM,TEMP,ZERO
      REAL R1MACH,ENORM
      DATA P1,P001,ZERO /1.0E-1,1.0E-3,0.0E0/
C***FIRST EXECUTABLE STATEMENT  LMPAR
      DWARF = R1MACH(1)
C
C     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
C     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 10 J = 1, N
         WA1(J) = QTB(J)
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA1(J) = ZERO
   10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
   60    CONTINUE
C
C     INITIALIZE THE ITERATION COUNTER.
C     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
C     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
C
      ITER = 0
      DO 70 J = 1, N
         WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
      DXNORM = ENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
C
C     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
C     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
C     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
C
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
         L = IPVT(J)
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
      DO 110 J = 1, N
         SUM = ZERO
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 100
         DO 90 I = 1, JM1
            SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
      TEMP = ENORM(N,WA1)
      PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
C
C     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
C
      DO 140 J = 1, N
         SUM = ZERO
         DO 130 I = 1, J
            SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
         L = IPVT(J)
         WA1(J) = SUM/DIAG(L)
  140    CONTINUE
      GNORM = ENORM(N,WA1)
      PARU = GNORM/DELTA
      IF (PARU .EQ. ZERO) PARU = DWARF/AMIN1(DELTA,P1)
C
C     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
C     SET PAR TO THE CLOSER ENDPOINT.
C
      PAR = AMAX1(PAR,PARL)
      PAR = AMIN1(PAR,PARU)
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM
C
C     BEGINNING OF AN ITERATION.
C
  150 CONTINUE
         ITER = ITER + 1
C
C        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
C
         IF (PAR .EQ. ZERO) PAR = AMAX1(DWARF,P001*PARU)
         TEMP = SQRT(PAR)
         DO 160 J = 1, N
            WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2)
         DO 170 J = 1, N
            WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
         DXNORM = ENORM(N,WA2)
         TEMP = FP
         FP = DXNORM - DELTA
C
C        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
C        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
C        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
C
         IF (ABS(FP) .LE. P1*DELTA
     1       .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP
     2            .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220
C
C        COMPUTE THE NEWTON CORRECTION.
C
         DO 180 J = 1, N
            L = IPVT(J)
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
         DO 210 J = 1, N
            WA1(J) = WA1(J)/SIGMA(J)
            TEMP = WA1(J)
            JP1 = J + 1
            IF (N .LT. JP1) GO TO 200
            DO 190 I = JP1, N
               WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
         TEMP = ENORM(N,WA1)
         PARC = ((FP/DELTA)/TEMP)/TEMP
C
C        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
C
         IF (FP .GT. ZERO) PARL = AMAX1(PARL,PAR)
         IF (FP .LT. ZERO) PARU = AMIN1(PARU,PAR)
C
C        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
C
         PAR = AMAX1(PARL,PAR+PARC)
C
C        END OF AN ITERATION.
C
         GO TO 150
  220 CONTINUE
C
C     TERMINATION.
C
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
C
C     LAST CARD OF SUBROUTINE LMPAR.
C
      END
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      REAL A(LDA,N),SIGMA(N),ACNORM(N),WA(N)
      INTEGER I,J,JP1,K,KMAX,MINMN
      REAL AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      REAL R1MACH,ENORM
      DATA ONE,P05,ZERO /1.0E0,5.0E-2,0.0E0/
C***FIRST EXECUTABLE STATEMENT  QRFAC
      EPSMCH = R1MACH(4)
C
C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
         ACNORM(J) = ENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
         AJNORM = ENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
C
C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C        AND UPDATE THE NORMS.
C
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*SQRT(AMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            SIGMA(K) = ENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE QRFAC.
C
      END
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
      INTEGER N,LDR
      INTEGER IPVT(N)
      REAL R(LDR,N),DIAG(N),QTB(N),X(N),SIGMA(N),WA(N)
      INTEGER I,J,JP1,K,KP1,L,NSING
      REAL COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      DATA P5,P25,ZERO /5.0E-1,2.5E-1,0.0E0/
C***FIRST EXECUTABLE STATEMENT  QRSOLV
      DO 20 J = 1, N
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10       CONTINUE
         X(J) = R(J,J)
         WA(J) = QTB(J)
   20    CONTINUE
C
C     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
C
      DO 100 J = 1, N
C
C        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
C        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
C
         L = IPVT(J)
         IF (DIAG(L) .EQ. ZERO) GO TO 90
         DO 30 K = J, N
            SIGMA(K) = ZERO
   30       CONTINUE
         SIGMA(J) = DIAG(L)
C
C        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
C        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
C        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
C
         QTBPJ = ZERO
         DO 80 K = J, N
C
C           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
C
            IF (SIGMA(K) .EQ. ZERO) GO TO 70
            IF (ABS(R(K,K)) .GE. ABS(SIGMA(K))) GO TO 40
               COTAN = R(K,K)/SIGMA(K)
               SIN = P5/SQRT(P25+P25*COTAN**2)
               COS = SIN*COTAN
               GO TO 50
   40       CONTINUE
               TAN = SIGMA(K)/R(K,K)
               COS = P5/SQRT(P25+P25*TAN**2)
               SIN = COS*TAN
   50       CONTINUE
C
C           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
C           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
C
            R(K,K) = COS*R(K,K) + SIN*SIGMA(K)
            TEMP = COS*WA(K) + SIN*QTBPJ
            QTBPJ = -SIN*WA(K) + COS*QTBPJ
            WA(K) = TEMP
C
C           ACCUMULATE THE TRANFORMATION IN THE ROW OF S.
C
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 70
            DO 60 I = KP1, N
               TEMP = COS*R(I,K) + SIN*SIGMA(I)
               SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I)
               R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
C
C        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
C        THE CORRESPONDING DIAGONAL ELEMENT OF R.
C
         SIGMA(J) = R(J,J)
         R(J,J) = X(J)
  100    CONTINUE
C
C     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
C     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 110 J = 1, N
         IF (SIGMA(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA(J) = ZERO
  110    CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
         J = NSING - K + 1
         SUM = ZERO
         JP1 = J + 1
         IF (NSING .LT. JP1) GO TO 130
         DO 120 I = JP1, NSING
            SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
         WA(J) = (WA(J) - SUM)/SIGMA(J)
  140    CONTINUE
  150 CONTINUE
C
C     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
C
      DO 160 J = 1, N
         L = IPVT(J)
         X(L) = WA(J)
  160    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE QRSOLV.
C
      END
      REAL FUNCTION R1MACH(I)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL/20K,0/,LARGE/77777K,177777K/
C     DATA RIGHT/35420K,0/,DIVER/36020K,0/
C     DATA LOG10/40423K,42023K/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE91), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     R1MACH(1) = 1.17549435E-38
C     R1MACH(2) = 1.70141163E+38
C     R1MACH(3) = 5.960464478E-8
C     R1MACH(4) = 1.119209290E-7
C     R1MACH(5) = 3.01030010E-1
C
C     DATA SMALL(1) / 00040000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
C     *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
 	DATA SMALL(1) /       128 /			!for unix system, uncommented this line
 	DATA LARGE(1) /    -32769 /			!for unix system, uncommented this line
 	DATA RIGHT(1) /     13440 /			!for unix system, uncommented this line
 	DATA DIVER(1) /     13568 /			!for unix system, uncommented this line
 	DATA LOG10(1) / 547045274 /			!for unix system, uncommented this line
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     ASSUMING REAL*4 IS THE DEFAULT REAL
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256 /
C     DATA LARGE(1), LARGE(2) /    -1,   -129 /
C     DATA RIGHT(1), RIGHT(2) /     0,  26880 /
C     DATA DIVER(1), DIVER(2) /     0,  27136 /
C     DATA LOG10(1), LOG10(2) /  8347,  32538 /
C
C     MACHINE CONSTANTS FOR THE IBM PC - MICROSOFT FORTRAN
C
c     DATA SMALL(1) / #00800000 /		!for PC sytem, uncommented this line
c     DATA LARGE(1) / #7F7FFFFF /		!for PC sytem, uncommented this line
c     DATA RIGHT(1) / #33800000 /		!for PC sytem, uncommented this line
c     DATA DIVER(1) / #34000000 /		!for PC sytem, uncommented this line
c     DATA LOG10(1) / #3E9A209A /		!for PC sytem, uncommented this line
C
C     MACHINE CONSTANTS FOR THE IBM PC - PROFESSIONAL FORTRAN
C     AND LAHEY FORTRAN
C
C     DATA SMALL(1)/ Z'00800000' /
C     DATA LARGE(1)/ Z'7F7FFFFF' /
C     DATA RIGHT(1)/ Z'33800000' /
C     DATA DIVER(1)/ Z'34000000' /
C     DATA LOG10(1)/ Z'3E9A209A' /
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ( 'R1MACH -- I OUT OF BOUNDS',25,1,2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END
      SUBROUTINE RWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN)
      INTEGER N,LDR
      REAL ALPHA
      REAL R(LDR,N),W(N),B(N),COS(N),SIN(N)
      INTEGER I,J,JM1
      REAL COTAN,ONE,P5,P25,ROWJ,TAN,TEMP,ZERO
      DATA ONE,P5,P25,ZERO /1.0E0,5.0E-1,2.5E-1,0.0E0/
C***FIRST EXECUTABLE STATEMENT  RWUPDT
      DO 60 J = 1, N
         ROWJ = W(J)
         JM1 = J - 1
C
C        APPLY THE PREVIOUS TRANSFORMATIONS TO
C        R(I,J), I=1,2,...,J-1, AND TO W(J).
C
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            TEMP = COS(I)*R(I,J) + SIN(I)*ROWJ
            ROWJ = -SIN(I)*R(I,J) + COS(I)*ROWJ
            R(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J).
C
         COS(J) = ONE
         SIN(J) = ZERO
         IF (ROWJ .EQ. ZERO) GO TO 50
         IF (ABS(R(J,J)) .GE. ABS(ROWJ)) GO TO 30
            COTAN = R(J,J)/ROWJ
            SIN(J) = P5/SQRT(P25+P25*COTAN**2)
            COS(J) = SIN(J)*COTAN
            GO TO 40
   30    CONTINUE
            TAN = ROWJ/R(J,J)
            COS(J) = P5/SQRT(P25+P25*TAN**2)
            SIN(J) = COS(J)*TAN
   40    CONTINUE
C
C        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA.
C
         R(J,J) = COS(J)*R(J,J) + SIN(J)*ROWJ
         TEMP = COS(J)*B(J) + SIN(J)*ALPHA
         ALPHA = -SIN(J)*B(J) + COS(J)*ALPHA
         B(J) = TEMP
   50    CONTINUE
   60    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE RWUPDT.
C
      END
      SUBROUTINE SNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,GTOL,
     1   MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,IPVT,QTF,
     2   WA1,WA2,WA3,WA4)
      INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
      INTEGER IJUNK,NROW,IPVT(N)
      REAL FTOL,XTOL,GTOL,FACTOR,EPSFCN
      REAL X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),WA1(N),WA2(N),
     1     WA3(N),WA4(M)
      LOGICAL SING
      EXTERNAL FCN
      INTEGER I,IFLAG,ITER,J,L,MODECH
      REAL ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM,ONE,PAR,
     1     PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,SUM,TEMP,TEMP1,
     2     TEMP2,XNORM,ZERO
      REAL R1MACH,ENORM,ERR(M),CHKLIM
      DATA CHKLIM/.1E0/
      DATA ONE,P1,P5,P25,P75,P0001,ZERO
     1     /1.0E0,1.0E-1,5.0E-1,2.5E-1,7.5E-1,1.0E-4,0.0E0/
c     write(6,*)'x=',x(1),x(2),x(3),x(4),x(5)
C
C***FIRST EXECUTABLE STATEMENT  SNLS1
      EPSMCH = R1MACH(4)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NJEV = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR. N .LE. 0 .OR.
     1    M .LT. N .OR. LDFJAC .LT. N .OR. FTOL .LT. ZERO
     2    .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO
     3    .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (IOPT .LT. 3 .AND. LDFJAC .LT. M) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
C
C     EVALUATE THE FUNCTION AT THE STARTING POINT
C     AND CALCULATE ITS NORM.
C
      IFLAG = 1
      IJUNK = 1
      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = ENORM(M,FVEC)
C
C     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
C
      PAR = ZERO
      ITER = 1
C
C     BEGINNING OF THE OUTER LOOP.
C
   30 CONTINUE
C
C        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
         IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
C
C        CALCULATE THE JACOBIAN MATRIX.
C
      IF (IOPT .EQ. 3) GO TO 475
C
C     STORE THE FULL JACOBIAN USING M*N STORAGE
C
      IF (IOPT .EQ. 1) GO TO 410
C
C     THE USER SUPPLIES THE JACOBIAN
C
         IFLAG = 2
         CALL FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
         NJEV = NJEV + 1
C
C             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
C
         IF(ITER .GT. 1) GO TO 355
         IF(IFLAG .LT. 0) GO TO 300
C
C             GET THE INCREMENTED X-VALUES INTO WA1(*).
C
         MODECH = 1
         CALL CHKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C             EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
C
         IFLAG = 1
         CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC)
         NFEV = NFEV + 1
         IF(IFLAG .LT. 0) GO TO 300
         DO 350 I = 1, M
           MODECH = 2
           CALL CHKDER(1,N,X,FVEC(I),FJAC(I,1),LDFJAC,WA1,
     1          WA4(I),MODECH,ERR)
           IF(ERR(I) .GE. CHKLIM) GO TO 350
           CALL XERRWV( ' SNLS1--DERIVATIVE OF FUNCTION I1 MAY BE WRONG,
     1 ERR=R1 TOO CLOSE TO 0.' ,70,7,0,1,I,0,1,ERR(I),ERR(I))
350      CONTINUE
355      CONTINUE
C
         GO TO 420
C
C     THE CODE APPROXIMATES THE JACOBIAN
C
410      IFLAG = 1
         CALL FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
         NFEV = NFEV + N
  420    IF (IFLAG .LT. 0) GO TO 300
C
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
         CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
C
C        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
C        QTF.
C
         DO 430 I = 1, M
            WA4(I) = FVEC(I)
  430         CONTINUE
         DO 470 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 460
            SUM = ZERO
            DO 440 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  440          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 450 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  450          CONTINUE
  460       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  470       CONTINUE
         GO TO 560
C
C        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
C        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
C        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
C        N COMPONENTS IN QTF.
C
  475    DO 490 J = 1, N
            QTF(J) = ZERO
            DO 480 I = 1, N
               FJAC(I,J) = ZERO
  480          CONTINUE
  490        CONTINUE
         DO 500 I = 1, M
            NROW = I
            IFLAG = 3
            CALL FCN(IFLAG,M,N,X,FVEC,WA3,NROW)
            IF (IFLAG .LT. 0) GO TO 300
C
C            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
C
            IF(ITER .GT. 1) GO TO 498
C
C            GET THE INCREMENTED X-VALUES INTO WA1(*).
C
            MODECH = 1
            CALL CHKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED.
C
            IF(I .NE. 1) GO TO 495
C
C            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
495         CONTINUE
            MODECH = 2
            CALL CHKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR)
            IF(ERR(I) .GE. CHKLIM) GO TO 498
            CALL XERRWV( ' SNLS1--DERIVATIVE OF FUNCTION I1 MAY BE WRONG
     1, ERR=R1 TOO CLOSE TO 0.' ,70,7,0,1,I,0,1,ERR(I),ERR(I))
498         CONTINUE
C
            TEMP = FVEC(I)
            CALL RWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2)
  500       CONTINUE
         NJEV = NJEV + 1
C
C        IF THE JACOBIAN IS RANK DEFICIENT, CALL QRFAC TO
C        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
C
         SING = .FALSE.
         DO 510 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) SING = .TRUE.
            IPVT(J) = J
            WA2(J) = ENORM(J,FJAC(1,J))
  510       CONTINUE
         IF (.NOT.SING) GO TO 560
         CALL QRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
         DO 550 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 540
            SUM = ZERO
            DO 520 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  520         CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 530 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  530          CONTINUE
  540       CONTINUE
            FJAC(J,J) = WA1(J)
  550       CONTINUE
  560    CONTINUE
C
C        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
C
C        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
C        AND INITIALIZE THE STEP BOUND DELTA.
C
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = ENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
C
C        COMPUTE THE NORM OF THE SCALED GRADIENT.
C
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = AMAX1(GNORM,ABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
C
C        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
C
         IF (GNORM .LE. GTOL) INFO = 4
         IF (INFO .NE. 0) GO TO 300
C
C        RESCALE IF NECESSARY.
C
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = AMAX1(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
C
C        BEGINNING OF THE INNER LOOP.
C
  200    CONTINUE
C
C           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
C
            CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2,
     1                 WA3,WA4)
C
C           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = ENORM(N,WA3)
C
C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
            IF (ITER .EQ. 1) DELTA = AMIN1(DELTA,PNORM)
C
C           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = ENORM(M,WA4)
C
C           COMPUTE THE SCALED ACTUAL REDUCTION.
C
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C           COMPUTE THE SCALED PREDICTED REDUCTION AND
C           THE SCALED DIRECTIONAL DERIVATIVE.
C
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = ENORM(N,WA3)/FNORM
            TEMP2 = (SQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
C
C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C           REDUCTION.
C
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
C
C           UPDATE THE STEP BOUND.
C
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO)
     1            TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*AMIN1(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
C
C           TEST FOR SUCCESSFUL ITERATION.
C
            IF (RATIO .LT. P0001) GO TO 290
C
C           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = ENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
C
C           TESTS FOR CONVERGENCE.
C
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE) INFO = 1
            IF (DELTA .LE. XTOL*XNORM) INFO = 2
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (INFO .NE. 0) GO TO 300
C
C           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (ABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH
     1          .AND. P5*RATIO .LE. ONE) INFO = 6
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
C
C           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
C
            IF (RATIO .LT. P0001) GO TO 200
C
C        END OF THE OUTER LOOP.
C
         GO TO 30
  300 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      IF (INFO .LT. 0) CALL XERROR( 'SNLS1  -- EXECUTION TERMINATED BECA
     1USE USER SET IFLAG NEGATIVE.',63,1,1)
      IF (INFO .EQ. 0) CALL XERROR( 'SNLS1  -- INVALID INPUT PARAMETER.'
     1,34,2,1)
      IF (INFO .EQ. 4) CALL XERROR( 'SNLS1  -- THIRD CONVERGENCE CONDITI
     1ON, CHECK RESULTS BEFORE ACCEPTING.',70,1,1)
      IF (INFO .EQ. 5) CALL XERROR( 'SNLS1  -- TOO MANY FUNCTION EVALUAT
     1IONS.',40,9,1)
      IF (INFO .GE. 6) CALL XERROR(  'SNLS1  -- TOLERANCES TOO SMALL, NO
     1 FURTHER IMPROVEMENT POSSIBLE.',64,3,1)
c      write(6,*)'FVEC=',fvec(1),fvec(2),fvec(3),fvec(4),fvec(5)
      RETURN
C
C     LAST CARD OF SUBROUTINE SNLS1.
C
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      return !change to return to stop the hard crash with matlab
c     STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
c            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
c         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
c               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
c               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
c               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
c               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
c               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
c         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ( (LLEVEL.EQ.2).OR.( (LLEVEL.EQ.1).AND.(MKNTRL.EQ.2).and.
	1	(lerr.ne.9)) )
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
c            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
c               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
c            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
c            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
C     END OF SOURCE FOR SNLS1
