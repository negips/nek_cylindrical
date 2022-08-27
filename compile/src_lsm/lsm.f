!======================================================================
!     Level Set Method
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine lsm_march_phi

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      include 'TEST'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /scruz/ ta1,ta2,ta3,ta4

      real SSign, dummy
      common /scrch/ SSign(lv), dummy(lv)

      real distn, eps
      common /levsetm/ distn, eps
      real lsm_gd(lv,3)
      real smoothsign(lv)
      common /levsetmv/ lsm_gd,smoothsign


      integer i,j
      integer n

      real Seps         ! smoothening for Sign
      real gdn
      real steadyn
      real gl2norm,op_glsc2_wt,glmax
      logical ifds
      real gmax,gmaxdt
      real s1


      n = lx1*ly1*lz1*nelv

      call rzero3(lsm_gd(1,1),lsm_gd(1,2),lsm_gd(1,3),n)
      call copy(ta4,t,n)

      call gradm1(lsm_gd(1,1),lsm_gd(1,2),
     $            lsm_gd(1,3),t) 

      distn = sqrt(op_glsc2_wt(lsm_gd(1,1),lsm_gd(1,2),lsm_gd(1,3),
     $              lsm_gd(1,1),lsm_gd(1,2),lsm_gd(1,3),bm1)/volvm1)
      if (nid.eq.0) write(6,*) '|Distance|: ', istep, time, distn

      ifds = .false.
      if (ifds) then
        call dsavg(lsm_gd(1,1))
        call dsavg(lsm_gd(1,2))
        call dsavg(lsm_gd(1,3))
      endif

!     Create Smooth sign fuction
      Seps = 1.0e-2
      do i=1,n
        gdn = sqrt(lsm_gd(i,1)**2 + lsm_gd(i,2)**2 + lsm_gd(i,3)**2)
        SSign(i) = t(i,1,1,1,1)/(sqrt(t(i,1,1,1,1)**2 
     $                           + (gdn*Seps)**2))    ! sign function
      enddo
      
      call copy(Smoothsign,SSign,n)

      if (abs(distn-1.0).lt.1.0e-6) return

      do i=1,2000
        call copy(dummy,t,n)
!        call lsm_rk2(t,SSign)
!        call lsm_rk4(t,SSign)
        call lsm_ssprk3(t,SSign)
!        call lsm_ExpEuler(t,SSign)

        call gradm1(ta1,ta2,ta3,t) 
        if (ifds) then
          call dsavg(ta1)
          call dsavg(ta2)
          call dsavg(ta3)
        endif  
      
        distn = sqrt(op_glsc2_wt(ta1,ta2,ta3,ta1,ta2,ta3,bm1)/volvm1)

        do j=1,n
          s1       = ta1(i)**2 + ta2(i)**2
          if (ndim.eq.3) s1 = s1 + ta3(i)**2
          ta4(j) = (1.0-sqrt(s1))
        enddo
        call copy(tmp1,ta4,n)
        gmax = glmax(ta4,n)
        gmaxdt = gmax*dt
       
!        if (nid.eq.0) write(6,*) '|Distance-2|: ', istep, i,
!     $                           distn,gmax,gmaxdt

        if (abs(distn-1.0).lt.1.0e-08) exit

      enddo  

      return
      end subroutine lsm_march_phi
!---------------------------------------------------------------------- 
      subroutine lsm_rk2(d0,SSign)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'

      include 'TEST'

      real d0(lx1*ly1*lz1*lelv)
      real SSign(lx1*ly1*lz1*lelv)
      real dt2

      real ddt(lx1*ly1*lz1*lelt,4)
      common /scrmg/ ddt

      real ta4(lx1*ly1*lz1*lelv,4)
      common /scruz/ ta4

      real fil 
      common /scrcg/ fil(lx1,ly1,lz1,lelt)

      real tmpsol(lx1*ly1*lz1*lelt,3)
      common /scrsf/ tmpsol

      real gdn

      integer i,n
      logical ifds,iffil

      real distn
      real op_glsc2_wt,glmax,glamax
      real gmax,gmaxdt

      ifds  = .false.
      iffil = .false.
      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt

      ifield = 2

      call rzero3(ta4(1,1),ta4(1,2),ta4(1,3),n)
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),d0)
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        ta4(i,4) = (sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)-1.0)
      enddo
      call copy(tmp1,ta4(1,4),n)

      gmax   = glamax(ta4(1,4),n)
      gmaxdt = gmax*dt

      if (abs(gmax).lt.1.0e-8) return

      dt2 = 1.0e-3 ! 1.0e-5/gmax

!      if (nid.eq.0) write(6,*) '|Distance-2|: ', istep, i,
!     $                           distn,gmax,dt2

      do i=1,n
        gdn = sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
        ddt(i,1) = SSign(i)*(1.0 - gdn)
      enddo
      if (iffil) then
        call hpf_dist(fil,d0)
        call add2(ddt(1,1),fil,n)           ! already has negative sign
      endif  
      call col3(fil,ddt(1,1),bm1,n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),d0,dt2,n)

!     2nd step      
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),ta4(1,4))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        gdn       = sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
        ddt(i,2)  = SSign(i)*(1.0 - gdn)
        ddt(i,2)  = ddt(i,2) + ddt(i,1)
      enddo
      if (iffil) then
        call hpf_dist(fil,ta4(1,4)) 
        call add2(ddt(1,2),fil,n)
      endif        
      call col3(fil,ddt(1,2),bm1,n)
      call invertB(ta4(1,4),fil)
      call add2s2(d0,ta4(1,4),dt2/2.0,n)

      return
      end subroutine lsm_rk2
!---------------------------------------------------------------------- 

      subroutine lsm_rk4(d0,SSign)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'

      include 'TEST'

      real d0(lx1*ly1*lz1*lelv)
      real SSign(lx1*ly1*lz1*lelv)
      real dt2

      real ddt(lx1*ly1*lz1*lelt,4)
!      common /scrmg/ ddt

      real ta4(lx1*ly1*lz1*lelv,4)
!      common /scruz/ ta4

      real fil(lx1*ly1*lz1*lelt) 
!      common /scrcg/ fil

      real tmpsol(lx1*ly1*lz1*lelt,3)
!      common /scrsf/ tmpsol

      real gdn

      integer i,n
      logical ifds,iffil,ifdeal

      real distn
      real op_glsc2_wt,glmax,glamax
      real gmax,gmaxdt

      ifds = .false.
      iffil = .false.
      ifdeal = .false.

      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt

      ifield = 2

      call rzero3(ta4(1,1),ta4(1,2),ta4(1,3),n)
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),d0)
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        ta4(i,4) = (1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2))
      enddo
!      call glmax(ta4(1,4),n)
      call copy(tmp1,ta4(1,4),n)
      gmax   = glamax(ta4(1,4),n)
      gmaxdt = 1.0e-5

      if (abs(gmax).lt.1.0e-8) return

!      dt2 = 1.0e-3 !gmaxdt/gmax

!      if (nid.eq.0) write(6,*) '|Distance-3|: ', istep, 
!     $                           gmax,dt2

!     1st step      
      do i=1,n
        tmpsol(i,1)=1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,1),tmpsol,SSign)
      else
        call col3(ddt(1,1),tmpsol,SSign,n)
        call col2(ddt(1,1),bm1,n)
      endif  
      if (iffil) then
        call hpf_dist(fil,d0)
        call add2col2(ddt(1,1),fil,bm1,n)           ! already has negative sign
      endif        
      call copy(fil,ddt(1,1),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),d0,dt2/2.0,n)

!     2nd step
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),ta4(1,4))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        tmpsol(i,1)=1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,2),tmpsol,SSign)
      else
        call col3(ddt(1,2),tmpsol,SSign,n)
        call col2(ddt(1,2),bm1,n)
      endif  

      if (iffil) then
        call hpf_dist(fil,ta4(1,4)) 
        call add2col2(ddt(1,2),fil,bm1,n)
      endif 
      call copy(fil,ddt(1,2),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),d0,dt2/2.0,n)

!     3rd step      
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),ta4(1,4))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        tmpsol(i,1)=1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,3),tmpsol,SSign)
      else
        call col3(ddt(1,3),tmpsol,SSign,n)
        call col2(ddt(1,3),bm1,n)
      endif  

      if (iffil) then 
        call hpf_dist(fil,ta4(1,4))
        call add2col2(ddt(1,3),fil,bm1,n)
      endif  
      call copy(fil,ddt(1,3),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),d0,dt2,n)

!     4th step
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),ta4(1,4))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        tmpsol(i,1)=1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,4),tmpsol,SSign)
      else
        call col3(ddt(1,4),tmpsol,SSign,n)
        call col2(ddt(1,4),bm1,n)
      endif  
      if (iffil) then
        call hpf_dist(fil,ta4(1,4)) 
        call add2col2(ddt(1,4),fil,bm1,n)
      endif  
      do i=1,n
        fil(i) = (ddt(i,1) + 2.0*ddt(i,2) 
     $              + 2.0*ddt(i,3) + ddt(i,4))
      enddo  
      call invertB(ta4(1,4),fil)
      call add2s2(d0,ta4(1,4),dt2/6.0,n)


      return
      end subroutine lsm_rk4
!---------------------------------------------------------------------- 
      subroutine lsm_ssprk3(d0,SSign)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'

      include 'TEST'

      real d0(lx1*ly1*lz1*lelv)
      real SSign(lx1*ly1*lz1*lelv)
      real dt2

      real ddt(lx1*ly1*lz1*lelt,4)
      common /scrmg/ ddt

      real ta4(lx1*ly1*lz1*lelv,4)
      common /scruz/ ta4

      real fil(lx1*ly1*lz1*lelt) 
      common /scrcg/ fil

      real tmpsol(lx1*ly1*lz1*lelt,3)
      common /scrsf/ tmpsol

      real gdn

      integer i,n
      logical ifds,iffil,ifdeal

      real distn
      real op_glsc2_wt,glmax,glamax
      real gmax,gmaxdt

      real s1,s2

      ifds = .true.
      iffil = .false.
      ifdeal = .false.

      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt

      ifield = 2

      call rzero3(ta4(1,1),ta4(1,2),ta4(1,3),n)
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),d0)
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        if (ndim.eq.3) call dsavg(ta4(1,3))
      endif  
      do i=1,n
        s1       = ta4(i,1)**2 + ta4(i,2)**2
        if (ndim.eq.3) s1 = s1 + ta4(i,3)**2
        ta4(i,4) = (1.0-sqrt(s1))
      enddo
      call copy(tmp1,ta4(1,4),n)
      gmax   = glamax(ta4(1,4),n)
      gmaxdt = 1.0e-5

      if (abs(gmax).lt.1.0e-8) return

!      dt2 = 1.0e-3 !gmaxdt/gmax

!     1st step      
      do i=1,n
        s1       = ta4(i,1)**2 + ta4(i,2)**2
        if (ndim.eq.3) s1 = s1 + ta4(i,3)**2
        tmpsol(i,1) = (1.0-sqrt(s1))
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,1),tmpsol,SSign)
      else
        call col3(ddt(1,1),tmpsol,SSign,n)
        call col2(ddt(1,1),bm1,n)
      endif  
      if (iffil) then
        call hpf_dist(fil,d0)
        call add2col2(ddt(1,1),fil,bm1,n)           ! already has negative sign
      endif        
      call copy(fil,ddt(1,1),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),d0,dt2,n)
      call copy(tmpsol(1,1),ta4(1,4),n)             ! 1st soln: f1 = Euler(f0)

!     2nd step
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),tmpsol(1,1))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        if (ndim.eq.3) call dsavg(ta4(1,3))
      endif  
      do i=1,n
        s1       = ta4(i,1)**2 + ta4(i,2)**2
        if (ndim.eq.3) s1 = s1 + ta4(i,3)**2
        tmpsol(i,2) = (1.0-sqrt(s1))
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,2),tmpsol(1,2),SSign)
      else
        call col3(ddt(1,2),tmpsol(1,2),SSign,n)
        call col2(ddt(1,2),bm1,n)
      endif  

      if (iffil) then
        call hpf_dist(fil,tmpsol(1,1)) 
        call add2col2(ddt(1,2),fil,bm1,n)
      endif 
      call copy(fil,ddt(1,2),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),tmpsol(1,1),dt2,n)
      s2    = 1.0/4.0
      call cmult(ta4(1,4),s2,n)
      call copy(tmpsol(1,2),ta4(1,4),n)
      s2    = 3.0/4.0
      call add2s2(tmpsol(1,2),d0,s2,n)    ! 2nd soln: f2 = 3/4*f0 + 1/4*Euler(f1)

!     3rd step      
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),tmpsol(1,2))
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        if (ndim.eq.3) call dsavg(ta4(1,3))
      endif  
      do i=1,n
        s1       = ta4(i,1)**2 + ta4(i,2)**2
        if (ndim.eq.3) s1 = s1 + ta4(i,3)**2
        tmpsol(i,3) = (1.0-sqrt(s1))
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,3),tmpsol(1,3),SSign)
      else
        call col3(ddt(1,3),tmpsol(1,3),SSign,n)
        call col2(ddt(1,3),bm1,n)
      endif  

      if (iffil) then 
        call hpf_dist(fil,tmpsol(1,2))
        call add2col2(ddt(1,3),fil,bm1,n)
      endif  
      call copy(fil,ddt(1,3),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),tmpsol(1,2),dt2,n)
      s2    = 2.0/3.0
      call cmult(ta4(1,4),s2,n)
      call copy(tmpsol(1,3),ta4(1,4),n)
      s2    = 1.0/3.0
      call add2s2(tmpsol(1,3),d0,s2,n)    ! Final soln: fn+1 = 1/3*f0 + 2/3*Euler(f2)

      call copy(d0,tmpsol(1,3),n)


      return
      end subroutine lsm_ssprk3
!---------------------------------------------------------------------- 

      subroutine invertB(out1,in1)

      implicit none
      
      include 'SIZE'
      include 'MASS'

      real out1(1)
      real in1(1)

      integer i,n

      n=lx1*ly1*lz1*nelv

      call dssum(in1,nx1,ny1,nz1)
      call copy(out1,in1,n)
      call col2(out1,binvm1,n)

      return
      end subroutine invertB
!---------------------------------------------------------------------- 

      subroutine lsm_ExpEuler(d0,SSign)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'

      real d0(lx1*ly1*lz1*lelv)
      real SSign(lx1*ly1*lz1*lelv)
      real dt2

      real ddt(lx1*ly1*lz1*lelt,4)
!      common /scrmg/ ddt

      real ta4(lx1*ly1*lz1*lelv,4)
!      common /scruz/ ta4

      real tmpsol(lx1*ly1*lz1*lelt,3)
!      common /scrsf/ tmpsol

      real fil(lx1*ly1*lz1*lelv) 
!      common /scrcg/ fil

      real gdn

      integer i,n
      logical ifds,iffil,ifdeal

      real distn
      real op_glsc2_wt,glmax,glamax
      real gmax,gmaxdt

      ifds = .false.
      iffil = .false.
      ifdeal = .false.

      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt

      ifield = 2

      call rzero3(ta4(1,1),ta4(1,2),ta4(1,3),n)
      call gradm1(ta4(1,1),ta4(1,2),ta4(1,3),d0)
      if (ifds) then
        call dsavg(ta4(1,1))
        call dsavg(ta4(1,2))
        call dsavg(ta4(1,3))
      endif  
      do i=1,n
        ta4(i,4) = (sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)-1.0)
      enddo
      gmax   = glamax(ta4(1,4),n)
      gmaxdt = 1.0e-5

      if (abs(gmax).lt.1.0e-8) return

      dt2 = 1.0e-4 !gmaxdt/gmax

      if (nid.eq.0) write(6,*) '|Distance-3|: ', istep, 
     $                           gmax,dt2

!     Explicit Euler 
      do i=1,n
        tmpsol(i,1)=1.0-sqrt(ta4(i,1)**2 + ta4(i,2)**2 + ta4(i,3)**2)
      enddo
      if (ifdeal) then
        call dealias_u_c(ddt(1,1),tmpsol,SSign)
      else
        call col3(ddt(1,1),tmpsol,SSign,n)
        call col2(ddt(1,1),bm1,n)
      endif  
      if (iffil) then
        call hpf_dist(fil,d0)
        call add2col2(ddt(1,1),fil,bm1,n)           ! already has negative sign
      endif        
      call copy(fil,ddt(1,1),n)
      call invertB(ta4(1,4),fil)
      call add2s2(d0,ta4(1,4),dt2,n)


      return
      end subroutine lsm_ExpEuler
!---------------------------------------------------------------------- 

      subroutine hpf_dist(fil,fld)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'         ! param(110),(111)
      include 'TSTEP'         ! ifield
      include 'MASS'          ! BM1

      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      integer n

      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)

      integer hpf_kut
      real hpf_chi
      logical hpf_ifboyd

      integer nel

      integer icalld
      save icalld
      data icalld /0/

      real hpf_op(lx1,lx1)
      save hpf_op

      real fld(lx1,ly1,lz1,lelv)    ! input field
      real fil(lx1,ly1,lz1,lelv)    ! filtered field

      integer hpf_ks,hpf_ke

c---------------------------------------- 
!      if(.not. iffilter(ifield)) return

      hpf_kut = int(param(101))+1
      hpf_chi = -1.0*abs(param(103))
c     Boyd transform to preserve element boundary values is 
c     linearly unstable when used as forcing.
      hpf_ifboyd = .false.      

      hpf_ks = lx1-hpf_kut
      hpf_ke = hpf_ks+5

      nel = nelv
      n = nxyz*nel

!      if (hpf_chi.eq.0) return
!      if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'apply hpf ',
!     $                                 ifield, hpf_kut, hpf_chi

      if (icalld.eq.0) then
c       Create the filter transfer function
        call hpf_trns_fcn(hpf_filter,hpf_kut)
!        call hpf_trns_fcn2(hpf_filter,hpf_ks,hpf_ke)

c       Build the matrix to apply the filter function
c       to an input field
        call build_hpf_mat(hpf_op,hpf_filter,hpf_ifboyd)

c       Only initialize once    
        icalld=icalld+1 
      endif


c     Apply filter to temp/passive scalar fields
      call build_hpf_fld(fil,fld,
     $     hpf_op,lx1,lz1)

c     Multiply by filter weight (chi)
      call cmult(fil,hpf_chi,n)    

c     Multiply by Mass matrix    
c     and add to source term
!      call addcol3(bq(1,1,1,1,ifield-1),ta1,bm1,n)

      return
      end subroutine hpf_dist

c----------------------------------------------------------------------

      subroutine hpf_trns_fcn2(diag,ks,ke)

c      implicit none

      include 'SIZE'
      include 'PARALLEL'

      real diag(lx1*lx1)
      integer nx,ks,ke
      integer k0,kk,k
      integer kut

      real amp

c     Set up transfer function
c
      nx = lx1
      call ident   (diag,nx)
c
      k0 = max(1,ks)
      k1 = min(lx1,ke)
      kut=k1-k0                                  ! kut=additional modes
      do k=k0,lx1
        kk = k+nx*(k-1)
        amp = ((k-k0)**2 + 0.)/(kut**2 + 0.)     ! Normalized amplitude. quadratic growth
        diag(kk) = 1.-amp
        if (diag(kk).lt.0.0) diag(kk) = 0.0
      enddo

c     Output normalized transfer function
      k0 = lx1+1
      if (nio.eq.0) then
        write(6,6) 'HPF :',((1.-diag(k)), k=1,lx1*lx1,k0)
   6    format(a8,16f9.6,6(/,8x,16f9.6))
      endif

      return
      end subroutine hpf_trns_fcn2

c---------------------------------------------------------------------- 
      subroutine dealias_u_c(cku,u,Cz)

!     Compute dealiased form:  J^T Bf *JCz .Ju w/ correct Jacobians

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
!      include 'TOTAL'

      real cku(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf,wd2,jacm1d
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , jacm1d(ltd),tr(ltd)
     $             , wd2(ltd),uf(ltd)

      integer e,k
      integer iu,ic,ijc,ick

      integer nxyz1,nxyzc,nxyzd,nxyzu,nxyzj

      real zd,wd
      common /dealias1/ zd(lxd),wd(lxd)

      integer i,j,l

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
!      if (ifcf) nxyzc = nxyzd

      nxyzj = nxyz1


      iu  = 1    ! pointer to scalar field u
      ic  = 1    ! pointer to vector field C
      ijc = 1    ! pointer to scalar JACM1
      ick = 1    ! pointer to scalar cku 

      call zwgl (zd,wd,lxd)  ! zwgl -- NOT zwgll!

      if (if3d) then
        do k=1,lzd
        do j=1,lyd
        do i=1,lxd
           l = (k-1)*lyd*lxd + (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)*wd(k)
        enddo
        enddo
        enddo
      else
        do j=1,lyd
        do i=1,lxd
           l = (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)
        enddo
        enddo

      endif


      do e=1,nelv

!       Interpolate Convecting Field   
        call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Convected Field   
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Jacobian (Probably only needs to be done once) 
        call intp_rstd(jacm1d,jacm1(ijc,1,1,1),lx1,lxd,if3d,0) ! 0 --> forward

        do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
           tr(i) = wd2(i)*jacm1d(i)*uf(i)*fz(i)
        enddo
        call intp_rstd(cku(ick),tr,lx1,lxd,if3d,1) ! Project back to coarse

        ic  = ic  + nxyzc
        iu  = iu  + nxyzu
        ijc = ijc + nxyzj
        ick = ick + nxyz1

      enddo

      return
      end subroutine dealias_u_c
!-----------------------------------------------------------------------

      subroutine lsm_jacobian_action2(SSign,px,py,pz)

      implicit none

      include 'SIZE'
      include 'MASS'

!     Action of the Jacobian
      real Ax(lx1*ly1*lz1*lelv)

!     Smooth Sign
      real SSign(lx1*ly1*lz1*lelv)

!     gradients of perturbations      
      real px(lx1*ly1*lz1*lelv)
      real py(lx1*ly1*lz1*lelv)
      real pz(lx1*ly1*lz1*lelv)

      real phip(lx1*ly1*lz1*lelv)         ! perturbation

!     pert gradients      
      real ppx(lx1*ly1*lz1*lelv)
      real ppy(lx1*ly1*lz1*lelv)
      real ppz(lx1*ly1*lz1*lelv)
      real ppt(lx1*ly1*lz1*lelv)
      common /scruz/ ppx,ppy,ppz,ppt

!     Work arrays      
      real tb1(lx1*ly1*lz1*lelt)
      real tb2(lx1*ly1*lz1*lelt)
      real tb3(lx1*ly1*lz1*lelt)
      real tb4(lx1*ly1*lz1*lelt)
      common /scrmg/ tb1,tb2,tb3,tb4

      integer i,n
     
      n = lx1*ly1*lz1*nelv 

      call gradm1(ppx,ppy,ppz,phip) 
      call rzero(ppt,n)

      do i=1,n
        tb1(i) = px(i)**2 + py(i)**2
        if (ndim.eq.3) tb1(i) = tb1(i)+pz(i)**2

        px(i) = SSign(i)*px(i)/tb1(i)
        py(i) = SSign(i)*py(i)/tb1(i)
        pz(i) = SSign(i)*pz(i)/tb1(i)
      enddo


      return
      end subroutine lsm_jacobian_action2 
!---------------------------------------------------------------------- 

      subroutine lsm_jacobian_action3(Ax,px,py,pz,phip)

      implicit none

      include 'SIZE'
      include 'MASS'

!     Action of the Jacobian
      real Ax(lx1*ly1*lz1*lelv)

!     Smooth Sign
      real SSign(lx1*ly1*lz1*lelv)

!     gradients of perturbations      
      real px(lx1*ly1*lz1*lelv)
      real py(lx1*ly1*lz1*lelv)
      real pz(lx1*ly1*lz1*lelv)

      real phip(lx1*ly1*lz1*lelv)         ! perturbation

!     pert gradients      
      real ppx(lx1*ly1*lz1*lelv)
      real ppy(lx1*ly1*lz1*lelv)
      real ppz(lx1*ly1*lz1*lelv)
      real ppt(lx1*ly1*lz1*lelv)
      common /scruz11/ ppx,ppy,ppz,ppt

!     Work arrays      
      real tb1(lx1*ly1*lz1*lelt)
      real tb2(lx1*ly1*lz1*lelt)
      real tb3(lx1*ly1*lz1*lelt)
      real tb4(lx1*ly1*lz1*lelt)
      common /scrmg/ tb1,tb2,tb3,tb4

      integer i,n
     
      n = lx1*ly1*lz1*nelv 

      call gradm1(ppx,ppy,ppz,phip)
      call rzero(ppt,n)

      do i=1,n
!       grad\phi.grad\phi'
        ppt(i) = px(i)*ppx(i) + py(i)*ppy(i)
        if (ndim.eq.3) ppt(i) = ppt(i) + pz(i)*ppz(i)
      enddo

      call chsign(ppt,n)
      call col2(ppt,bm1,n)
      call invertB(Ax,ppt)
!      call copy(Ax,ppt,n)

!     Ax is now: Ax = -S0*(\grad\phi.\grad\phi')/(\grad\phi.\grad\phi)

      return
      end subroutine lsm_jacobian_action3 
!---------------------------------------------------------------------- 
