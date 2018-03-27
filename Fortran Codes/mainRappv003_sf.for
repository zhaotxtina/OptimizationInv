!default order of channels: P16H,P22H,P28H,P34H,P40H, A16H,A22H,A28H,A34H,A40H,    P16L,P22L,P28L,P34L,P40L,A16L,A22L,A28L,A34L,A40L
! the inversion is for two frequencies, for each frequency 2MHz,400kHz, there is one set of output, so the channels of H and L should be separate for inversion input
! Start developing Jan 19, 2015, based on ecaliper codes and subrountines from Arcwizard as well
! By Tina Zhao
	program vdiinv_2015
	implicit none
	integer i,j,k,l,m,n,irt
	integer, parameter:: max_m=20,large_number=100000,ntr=5,nd=100000
	real*4, parameter:: abst = -999.25 !absent value
      integer id_tool
      real*4 rm,bit_size
      integer total_number
      character*100 Input_filename
      CHARACTER Data_FilePath*200
      integer path_size
      real*4 t_start,t_end,tstart,t1
      character*300 filename
!      integer nd ! number of data points
      integer nch_desire  !selecting 10 channels for inversion 
      integer i_auto

      integer nch
      integer, allocatable:: index_channel(:),index_channel_tmp(:)
      real*4, allocatable::  weight_channel(:),weight_channel_tmp(:)
 
 !     real*4 depth(nd),dh_input(nd),rm_input(nd),rt_input(nd),ecc_input(nd)
 !     real*4 dh_bhc(nd),rm_bhc(nd), rt_bhc(nd), ecc_bhc(nd),err_bhc(nd),QC(nd)
 !     real*4 meas(nd,max_m)    

      real*4,allocatable:: depth(:),dh_input(:),rm_input(:),rt_input(:),epsr_input(:)
   !  here have two groups of output for hd,rt,epsr,rm for two frequencies   
      real*4,allocatable:: dh_bhc2(:),rm_bhc2(:), rt_bhc2(:), epsr_bhc2(:),err_bhc2(:),QC2(:)
      real*4,allocatable:: dh_bhc4(:),rm_bhc4(:), rt_bhc4(:), epsr_bhc4(:),err_bhc4(:),QC4(:)
!      real*4,allocatable:: std_hd2(:),std_rm2(:),std_rt2(:),std_epsr2(:),std_hd4(:),std_rm4(:),std_rt4(:),std_epsr4(:)
      real*4,allocatable:: meas(:,:),meas_tmp(:,:)  !input measurement channels
 
      integer index_init_rt(max_m)
      data index_init_rt /20,15,19,14,18,13,10,5,9,4,8,3,17,12,16,11,7,2,6,1/ !preferred channels to assign initial Rt in cases it is not provided

!     channel order from las file: A16H A16L A22H A22L A28H A28L A34H A34L A40H A40L    P16H P16L P22H P22L P28H P28L P34H P34L P40H P40L      
!     integer index_las(20)
!     data index_las /6,16,7,17,8,18,9,19,10,20,  1,11,2,12,3,13,4,14,5,15/    !Default channel orders from the las file
!     index corresponds to index 1-10 is for 2MHz: P16-P40H, A16-A40H; index 11-20 is for 400kHz: P16-P40L,A16-A40L

!  the above is for Rapp input, now the PSAD input channels are: 
!  index 1-10 is for 2MHz: P16-P40H, A16-A40H; index 11-20 is for 400kHz: P16-P40L,A16-A40L

      CHARACTER*20000  Stmp
      integer mode  !mode of inversion  still under inverstigation, main purpose is for Rt and epsr
                    !what is the OBM impact on mode selection
                    !mode=1, invert for Rt and epsr
                    !mode=2, invert for Rt, epsr, dh, rm
                    
	integer is_OBM ! OBM indicator
	real*4  misfit_epsr !Ecc inversion threshhold  ! need to change, does epsr needs threshold check?
      
!      real*4 atmp(nd)  !Borehole temperature, used for estimating (tmeperature dependent) mud resistivity
      real*4 test(21),tmpv(21)  ! what are these variables?

      path_size=200
      ! id_tool = 4
      !id_tool = 7
      
      ! initialize weight_channel, and index channel
      weight_channel=0
      index_channel=0
      
        open(111,file='parameter_rt_v675_040115wbm.inp')
        ! open(111,file='parameter_v5wbmofgm2_051815.txt')
        ! open(111,file='parameter_v6wbmofgm2_062215.txt')
        ! open(111,file='parameter_expert_eco_dual051315.txt')
         read(111,*) is_OBM, mode, id_tool
         read(111,*) Input_filename
         read(111,*) misfit_epsr  !misfit_ecc=1.0e-3
         read(111,*) bit_size    !bitsize can be actual bit size or the reamer size, so it can be different from id_tool
         read(111,*) rm          !user provided parameter rm, used for mud ressitivity in case mud is not inverted, 
                                 !between 0.015 to 10 ohm.m for WBM, fixed at 10000ohm.m for OBM
         read(111,*) nch         ! normally 20 channels
         
         allocate (index_channel(nch))        ! allocate index_channels and weight_channels for 20 channels
         allocate (index_channel_tmp(nch))  !not defined yet
         allocate (weight_channel_tmp(nch)) ! not defined yet

         allocate (weight_channel(nch))

         read(111,*) index_channel(1:nch)
         read(111,*) weight_channel(1:nch)
      close(111)          ! finish reading parameter_expert.inp
      
      
      open(111,file='auto.inp')
         read(111,*) nch_desire,i_auto    ! noe set to 20, 0, so what is i_auto?
      close(111)
      
       ! check if it's OBM, then mode with Rm output will be given warnnings 

      if (is_OBM.eq.1) rm=1000.0       ! is is_OBM=1, is OBM, if is 0, it's WBM

      if(is_OBM.eq.1) then      ! not sure that many modes yet, but if it's OBM, cannot be inverted if mode asks to invert Rm, so far mode=2 in my
       !   if(mode.eq.2.or.mode.eq.3.or.mode.eq.5) then
          if(mode.eq.2 ) then
             write(*,*) 'Rm cannot be inverted for OBM'
             stop
          endif
      endif



      if(is_OBM) rm_input=1000.0
      
         ! find how many inputs are there

      OPEN(8, file=Input_filename,SHARE='DENYNONE')
      Stmp='   '
      
	DO  while (Stmp(1:2).NE.'~A')
        READ(8,'(A)')  Stmp
	ENDDO
   
      
!      OPEN(8, file=Input_filename,SHARE='DENYNONE')
      
      !finding out the total sample numbre total_number
      total_number=1
	do i=1,large_number
	  if (EOF(8)) goto 300
        READ(8,*) test(1)
	  total_number=total_number+1;
	enddo

300   CONTINUE 
	close(8)
	
	total_number=total_number-1    ! don't know about the difference of imput and bhc yet, later in other subrountines...
!	total_number= 27
      allocate( depth(total_number),dh_input(total_number),rm_input(total_number),rt_input(total_number),epsr_input(total_number))
      allocate( dh_bhc2(total_number),rm_bhc2(total_number), rt_bhc2(total_number), epsr_bhc2(total_number),
     &           err_bhc2(total_number),QC2(total_number))             ! for 2MHz output
!       allocate( std_hd2(total_number),std_rm2(total_number), std_rt2(total_number), std_epsr2(total_number))
                          
      allocate( dh_bhc4(total_number),rm_bhc4(total_number), rt_bhc4(total_number), epsr_bhc4(total_number),
     &           err_bhc4(total_number),QC4(total_number))             ! for 400kHz output
!      allocate( std_hd4(total_number),std_rm4(total_number),std_rt4(total_number),std_epsr4(total_number))
      
      allocate( meas(total_number,nch))
      allocate( meas_tmp(total_number,nch))

        ! set default or absent values for the four parameters
      dh_input=abst
      rt_input=abst
      rm_input=abst
      epsr_input=abst


      meas=abst
      
      
         ! input from LAS file the measured data
      OPEN(8, file=Input_filename,SHARE='DENYNONE') 
      Stmp='   '
     
	DO  while (Stmp(1:2).NE.'~A')
      READ(8,'(A)')  Stmp
	ENDDO
   
	do i=1,total_number

          READ(8,*) depth(i),meas(i,1:nch) ! ,dh_input(i),rm_input(i) ! , dh_input(i),rm_input(i) !read in measurement channels, data from Hong doesn't have Rm at the end 
! the rt_input is supposed to be computed later in the subroutine Vdie_inv_v006.for, because most input data won't have it
!          ! Convert to default order P16H,P22H,P28H,P34H,P40H, A16H,A22H,A28H,A34H,A40H,  P16L,P22L,P28L,P34L,P40L,A16L,A22L,A28L,A34L,A40L
!          do j=1,5
!             test(15+j)=tmpv(2*j)        !A16L,A22L,...,A40L
!             test(10+j)=tmpv(10+2*j)     !P16L,P22L,...,P40L
!             test(5+j) =tmpv(2*j-1)      !A16H,A22H,...,A40H
!             test(j)=tmpv(10+2*j-1)      !P16H,P22H,...,P40H
!          enddo
!      
!          ! Save to meas, which only contains selected channels
!          do j=1,nch
!             meas(i,j)=test(index_channel(j))
!          enddo

	enddo
30    CONTINUE 
	close(8)


! assigning initial values in cases not provided

       do i=1,total_number 
          if(abs(dh_input(i)-abst).lt.1.0e-4)  dh_input(i)=bit_size
          if(abs(rm_input(i)-abst).lt.1.0e-4)  rm_input(i)=rm  ! where rm is from?
          if(abs(rt_input(i)-abst).lt.1.0e-4)  rt_input(i)=100  ! old flow
           rt_input(i)=100  ! make it constant now, and replaced later in subroutine Vdie_inv_v009.for
         if(abs(epsr_input(i)-abst).lt.1.0e-4)   epsr_input=10.0   ! just assing a random guess for epsr
          
          if(is_OBM) rm_input=1000.0
       enddo

! Assigning initial Rt 

!       irt=9
       
!       do i=1,20
!         do j=1,nch
!           if (index_init_rt(i).eq.index_channel(j)) then
!              irt=j 
!              go to 40
!           endif
!         enddo
!       enddo
!
!40     continue

!       rt_input=meas(:,irt)       


! reordering the measuremetn channels

    
  
      do m=1,nch
         do j=1,nch
            if (index_channel(j).eq.m) then
               index_channel_tmp(m)=m
!               meas_tmp(:,i)=meas(:,j)               
               weight_channel_tmp(m)=weight_channel(j)
            endif
         enddo
      enddo    

      rt_input(:)=meas(:,10)    !A40H,     this is included in input last column

c      Data_FilePath = 'D:\work\eCaliper\tables\'
c       Data_FilePath = 'D:\table_venus_diinv\'
      Data_FilePath = 'D:\Documents\Venus_AP_study\dielectric\fulltablecode_2015_sf\'

      CALL CPU_TIME(t1)
      t_start=t1

      open(1000,file='before_passing.dat')
       write(1000,*) id_tool,bit_size,mode,   ! tool ID, bit_size, inversion mode
     &     nch,index_channel,weight_channel,          ! number of input channels, channels index and weighting 
     &     total_number,depth,meas,nch_desire,i_auto, ! No. measurement points, depth and measurement channels
     &     dh_input,rm_input,rt_input,epsr_input,      ! model parameters in case they are not inverted
     &     Data_FilePath,path_size, is_OBM,misfit_epsr, ! location of data files;  OBM indicator, misfit threshhold for Epsr inversion
     &     dh_bhc2,rm_bhc2,rt_bhc2,epsr_bhc2,err_bhc2,QC2,  ! inversion results                        
     &     dh_bhc4,rm_bhc4,rt_bhc4,epsr_bhc4,err_bhc4,QC4  ! inversion results
      close(1000)

      
      call Vdie_inv_v003(id_tool,bit_size,mode,         ! tool ID, bit_size, inversion mode
     &     nch,index_channel_tmp,weight_channel_tmp,    ! number of input channels, channels index and weighting 
     &     total_number,depth,meas,nch_desire,i_auto, ! No. measurement points, depth and measurement channels
     &     dh_input,rm_input,rt_input,epsr_input,      ! model parameters in case they are not inverted
     &     Data_FilePath,path_size, is_OBM,misfit_epsr,  ! location of data files; OBM indicator, misfit threshhold for Epsr inversion?
     &     dh_bhc2,rm_bhc2,rt_bhc2,epsr_bhc2,err_bhc2,QC2,                         
     &     dh_bhc4,rm_bhc4,rt_bhc4,epsr_bhc4,err_bhc4,QC4)   ! inversion results
!     &     std_hd2,std_rm2,std_rt2,std_epsr2,   ! uncertainty results
!     &     std_hd4,std_rm4,std_rt4,std_epsr4)   ! uncertainty results


!      call ARCEcaliper(id_tool,bit_size,mode,         ! tool ID, bit_size, inversion mode
!     &     nch,index_channel,weight_channel,          ! number of input channels, channels index and weighting 
!     &     total_number,depth,meas,nch_desire,i_auto, ! No. measurement points, depth and measurement channels
!     &     dh_input,rm_input,rt_input,ecc_input,      ! model parameters in case they are not inverted
!     &     Data_FilePath,path_size,                   ! location of data files
!     &     is_OBM,misfit_ecc,                         ! OBM indicator, misfit threshhold for Ecc inversion
!     &     dh_bhc,rm_bhc,rt_bhc,ecc_bhc,err_bhc,QC)   ! inversion results



      CALL CPU_TIME(t1)
      t_end=t1

      write(*,*) t_end,t_start,t_end-t_start
       filename='wizard2M_v6wbm_040115rtv3.txt'
      !filename='wizard2M_v6wbmofgm2_080715.txt'
      ! filename='wizard2M_v6wbmofgm2_080615.txt'
      ! filename='wizard2M_eco_psad_owbm_0623.txt'
       open(unit=11,file=filename,status='unknown') 
       
       filename='wizard400k_v6wbm_040115rtv3.txt'
      ! filename='wizard400k_v6wbmofgm2_080715.txt'
      ! filename='wizard400k_v6wbmofgm2_080615.txt'
      ! filename='wizard400k_eco_psad_owbm_0623.txt'
       open(unit=13,file=filename,status='unknown') 
      do i=1,total_number

!          write(11,1001) depth(i),rt_bhc2(i),dh_bhc2(i),rm_bhc2(i),epsr_bhc2(i),err_bhc2(i),
!     1                QC2(i),std_hd2(i),std_rm2(i),std_rt2(i),std_epsr2(i)
       
!        write(13,1003) depth(i),rt_bhc4(i),dh_bhc4(i),rm_bhc4(i),epsr_bhc4(i),err_bhc4(i),
!     1                QC4(i),std_hd4(i),std_rm4(i),std_rt4(i),std_epsr4(i)
        
         write(11,1001) depth(i),rt_bhc2(i),dh_bhc2(i),rm_bhc2(i),epsr_bhc2(i),err_bhc2(i),
     1                QC2(i)
       
        write(13,1003) depth(i),rt_bhc4(i),dh_bhc4(i),rm_bhc4(i),epsr_bhc4(i),err_bhc4(i),
     1                QC4(i)
       
       
      enddo 
      
       close(11)      
       close(13)
      

1001  format(7(1pe16.7))  
1003  format(7(1pe16.7)) 
      

      
!      write(100,'(A)') '         depth,                dh_bhc,             rm_bhc,            rt_bhc,            
!     & ecc_bhc,             err_bhc,              QC '
!      do i=1, nd
!      write(100,'(100e20.8)') depth(i), dh_bhc(i),rm_bhc(i), rt_bhc(i),ecc_bhc(i), err_bhc(i),QC(i) 
!      enddo
      
      deallocate( depth,dh_input,rm_input,rt_input,epsr_input)
!      deallocate( dh_bhc2,rm_bhc2, rt_bhc2, epsr_bhc2,err_bhc2,QC2,std_hd2,std_rm2,std_rt2,std_epsr2) 
!      deallocate( dh_bhc4,rm_bhc4, rt_bhc4, epsr_bhc4,err_bhc4,QC4,std_hd4,std_rm4,std_rt4,std_epsr4) 
      deallocate( dh_bhc2,rm_bhc2, rt_bhc2, epsr_bhc2,err_bhc2,QC2) 
      deallocate( dh_bhc4,rm_bhc4, rt_bhc4, epsr_bhc4,err_bhc4,QC4)      
      deallocate (meas,meas_tmp)

     
      deallocate (index_channel)
      deallocate (weight_channel)

      deallocate (index_channel_tmp)
      deallocate (weight_channel_tmp)

      
      stop
      end