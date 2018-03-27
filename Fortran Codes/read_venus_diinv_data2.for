! read in the table for the given tool size
!      program  read_data
!      implicit none
	subroutine 	read_venus_diinv_data(filename,nrt,nrm,nhd,nepsr,ntr,ad_b,ps_b,ad_a,ps_a) 
 !read in table Rt,Rm,HD and epr <---> unmixed ARC channels (in db/deg)

! Input grid sizes
	integer nrt,nrm,nhd,nepsr
!	integer necc_hd(nhd)
! output -- unmixed, non-aircalled raw tool response loaded from the bhc table
!      real*4   ps_a(26,19,19,5,27),ad_a(26,19,19,5,27)
!      real*4   ps_b(26,19,19,5,27),ad_b(26,19,19,5,27)

!      real*4, allocatable:: ps_a(:,:,:,:,:),ad_a(:,:,:,:,:)
      real*4 ps_a(nrt,nhd,nrm,ntr,nepsr),ad_a(nrt,nhd,nrm,ntr,nepsr)
      real*4 ps_b(nrt,nhd,nrm,ntr,nepsr),ad_b(nrt,nhd,nrm,ntr,nepsr)
!
	character*300 filename,tmp
	integer iHD, iepsr,irt,irm,icase, itr,i
      
!      nrt= 26
!      nepr=27
!      nhd=19
!      nrm=19

!      filename='v5_di_db2.dat'
      open(unit=11,status='unknown', file=filename)
      
      
!      do i=1,9
!        read(11,*)
!      enddo
      
      icase=0;
      do iHD=1,nhd
    
          do iepsr=1,nepsr
!                do i=1,4
!                  read(11,*)
!                enddo                
               do irm=1,nrm ! loop over R_m
                  do irt=1,nrt ! loop over R_t
                  icase=icase+1;
            read(11,*) ((ad_b(irt,iHD,irm,itr,iepsr),ps_b(irt, iHD,irm,itr,iepsr)),itr=1,5),
     &      ((ad_a(irt,iHD,irm,itr,iepsr),ps_a(irt,iHD,irm,itr,iepsr)),itr=1,5)
                   
                     enddo
                  
                enddo
              
          enddo
   
      enddo
      
      

      close (11)
    
      

	return
	end
