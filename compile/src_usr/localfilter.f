!====================================================================== 
!     Modified hpf.f to use an adhoc local filter
!          
!
!======================================================================       
      subroutine localfilterfld(ux)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'         ! param(110),(111)
      include 'TSTEP'         ! ifield
      include 'MASS'          ! BM1

      real ux(lx1,ly1,lz1,lelv)

      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      integer n

      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real filter(lm2)

      integer kut
      logical ifboyd
      real wght

      integer nel

      integer icalld
      save icalld
      data icalld /0/

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      real f_op(lx1,lx1)
      save f_op

      kut  =  int(uparam(4))
      wght = -abs(uparam(5))

!     Boyd transform to preserve element boundary values is 
!     linearly unstable when used as forcing.
      ifboyd = .false. 

      nel = nelv
      n = nxyz*nel

      if (icalld.eq.0) then
!       Create the filter transfer function
        call filter_trns_fcn(filter,kut)

!       Build the matrix to apply the filter function
!       to an input field
        call build_hpf_mat(f_op,filter,ifboyd)

!       Only initialize once    
        icalld=icalld+1 
      endif

!     Apply the filter.
!     ta1 has the high pass filtered field
      call build_hpf_fld(ta1,ux,f_op,lx1,lz1)
!      call copy(ux,ta1,n)
      call add2s2(ux,ta1,wght,n)     ! if we want the low pass filtered field

      return
      end subroutine localfilterfld

!----------------------------------------------------------------------
      subroutine filter_trns_fcn(diag,kut)

c      implicit none

      include 'SIZE'
      include 'PARALLEL'

      real diag(lx1*lx1)
      integer nx,k0,kut,kk,k

      real amp

c     Set up transfer function
c
      nx = lx1
      call ident   (diag,nx)
c
      kut=kut                                     ! kut=additional modes
      k0 = nx-kut
      do k=k0+1,nx
        kk = k+nx*(k-1)
        amp = ((k-k0)*(k-k0)+0.)/(kut*kut+0.)     ! Normalized amplitude. quadratic growth
        diag(kk) = 1.-amp
      enddo

c     Output normalized transfer function
      k0 = lx1+1
      if (nio.eq.0) then
        write(6,6) 'FILT :',((1.-diag(k)), k=1,lx1*lx1,k0)
   6    format(a8,16f9.6,6(/,8x,16f9.6))
      endif

      return
      end subroutine filter_trns_fcn

!---------------------------------------------------------------------- 


