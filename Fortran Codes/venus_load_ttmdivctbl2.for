c       subroutine arc_load_ttmbhctbl(id_tool)
c       purpose: to load the ARC borehole correction tables 
c	  from TTM model into the common block 
c       Author:	P.Wu/SPC Oct. 24, 2001
c
	subroutine venus_load_ttmdivtbl(id_tool,Data_FilePath,path_size)
c
	character*300 filename
	integer path_size
      real*4  hd5(19),hd6(19)
	character Data_FilePath*path_size
	include 'venus_ttmdivtab_def.inc'


! assign Rt grid
      data (rt(i),i=1,26)/ 0.10,0.15,0.20,0.30,0.50,0.70,1.0,1.5,2.0,
     &  3.0,5.0,7.0,10.0,15.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,
     &  300.0,500.0,700.0,1000.0,1500.0/

! assign Rm grid
      data (rm(i),i=1,19)/ 0.015,0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3,
     &                0.5,0.7,1.0,1.5,2.0,3.0,5.0,7.0,10.0,1000.0/
    

     

! assign hd grid
        ! V475
      data (hd5(i),i=1,19)/5.25,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,
     &    11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0/

      ! for v675
      data (hd6(i),i=1,19)/7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12.0,13.0,14.0,
     1          15.0,16.0,17.0,18.0,19.0,20.0,22.0,24.0/

! Assign epsr grid

     
      data (epsr(i),i=1,27)/1.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,
     &  90.0,100.0,125.0,150.0,175.0,200.0,225.0,250.0,275.0,300.0,
     &  500.0, 1000.0,2000.0,5000.0,10000.0,20000.0,50000.0,100000.0/


      
      nrt_t=nrt
      nrm_t=nrm
      nhd_t=nhd
      nepsr_t=nepsr
      
!     Initialize the table	
	ps_a=-999.25
	ps_b=-999.25
	ad_a=-999.25
	ad_b=-999.25

	
! load in the table for Venus tool
     
      if (id_tool .eq. 4) then
	
      filename='a5_di_db2_rt26.dat'   ! arc475 table
      do k=1,19
         hd(k)=hd5(k)
      enddo
      
      else if (id_tool .eq. 5) then
  
          filename='v5_di_db2.dat'   ! v475 table
        do k=1,19
            hd(k)=hd5(k)
        enddo
        
        else if (id_tool .eq. 6) then
     
          filename='v6_di_db2_hd19.dat'   ! v675 table
        do k=1,19
            hd(k)=hd6(k)
       enddo
      
      else if (id_tool .eq. 7) then
      
          filename='a6_di_db2.dat'   ! arc675 table
        do k=1,19
            hd(k)=hd6(k)
       enddo
        
      endif
	filename=trim(Data_FilePath)//filename
!!	write(*,*) filename
   
       
  
      !read in table Rt,Rm,HD and epsr <---> unmixed ARC channels (in db/deg)

      call read_venus_diinv_data(filename,nrt,nrm,nhd,nepsr,ntr,ad_b,ps_b,ad_a,ps_a) 
      
      !  test if read in the correct data
c       print *," first example is", ad_b(4,5,7,4,4),ps_b(4,5,7,4,4)     ! ad_b(irt,iHD,irm,itr,iepr)
c       print *," second example is", ad_a(15,12,17,1,23),ps_a(15,12,17,1,23)  ! ad_b(irt,iHD,irm,itr,iepr)
       
      
      
! fill in the blank spots with unrealistic eccentricity

!      do itr=1,ntr
!            do irt=1,nrt
!               do irm=1,nrm
!                  do ihd=1,nhd 
!      !               ps_a(nrt,nhd,nrm,ntr,necc)
!                      ii=necc_hd(ihd)
!                      ps_a(irt,ihd,irm,itr,ii:necc)=ps_a(irt,ihd,irm,itr,ii-1)
!                      ad_a(irt,ihd,irm,itr,ii:necc)=ad_a(irt,ihd,irm,itr,ii-1)
!                  end do
!               end do
!            enddo
!      enddo
! END FILLING




! **********************************************************************************
!   reload in the table from the tubtumud

!      if(0) then
!!     load in tubtomud table  
!            open(100,file='ttmtable.txt')
!!            write(100,'(1x,A, 100e20.8)') '%rt=', rt
!            read(100,*)
!            read(100,*)
!            read(100,*)
!            read(100,*)
!            read(100,*) 

!            do ihd=1,nhd
!               do  irt=1,nrt
!                   do irm=1,nrm
!                      read(100,*)  ((ad_b(irt,ihd,irm,itr,1),ps_b(irt,ihd,irm,itr,1)),itr=1,ntr),
!     &                             ((ad_a(irt,ihd,irm,itr,1),ps_a(irt,ihd,irm,itr,1)),itr=1,ntr)  
!                   enddo
!               enddo
!            enddo

!            close(100)

!       endif
!**********************************************************
!!  Testing the table
!      if(0) then
!!          Write table in the format of tiltecc table 
!           open(100,file='ttmtablefull_test.txt')
!!          write(100,'(1x,A, 100e20.8)') '%rt=', rt
!           write(100,'(A, 100e20.8)') '%rt=', rt
!           write(100,'(A, 100e20.8)') '%rm=',rm
!           write(100,'(A, 100e20.8)') '%hd=',hd
!           write(100,'(A)')  '%orders (((irm=1,nrm), irt=1,nrt), ihd=1,nhd)'
!           write(100,'(A100)') '%A16L P16L A22L P22L A28L P28L A34L P34L A40L P40L A16H P16H A22H P22H
!     & A28H P28H A34H P34H A40H P40H' 

c            do ihd=1,nhd
c               do  irt=1,nrt
c                   do irm=1,nrm
c                      write(100,'(20e20.8)')  ((ad_b(irt,ihd,irm,itr,1),ps_b(irt,ihd,irm,itr,1)),itr=1,ntr),
c     &                                        ((ad_a(irt,ihd,irm,itr,1),ps_a(irt,ihd,irm,itr,1)),itr=1,ntr)  
c                   enddo
c               enddo
c            enddo
          
c              !  write back the whole lookup table, and compare to v5_di_db2.dat, almost identical, accuracy differs
c            do ihd=1,nhd
c               do iepsr=1,nepsr
                    
c                 do irm=1,nrm
c                   do  irt=1,nrt
                  
c                  write(100,'(20f20.8)')  ((ad_b(irt,ihd,irm,itr,iepsr),ps_b(irt,ihd,irm,itr,iepsr)),itr=1,ntr),
c     &                                  ((ad_a(irt,ihd,irm,itr,iepsr),ps_a(irt,ihd,irm,itr,iepsr)),itr=1,ntr)  
c                   enddo
c                enddo
c               enddo
c           enddo 

c            close(100)
!      endif


	return
	end

