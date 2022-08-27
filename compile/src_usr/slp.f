!======================================================================
!     Routines for introducing slip velocities across element boundaries
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine slp_mark_faces

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'

      integer icalld
      save    icalld
      data    icalld /0/

      integer i,j,n,m
      integer nfaces

      real ut1(lx1,ly1,lz1,lelt)
      real ut2(lx1,ly1,lz1,lelt)
      real ut3(lx1,ly1,lz1,lelt)
      integer flowtype(lelt)

      common /testvel1/ ut1,ut2,ut3,flowtype

      real ut4(lx1,ly1,lz1,lelt)
      real ut5(lx1,ly1,lz1,lelt)
      real ut6(lx1,ly1,lz1,lelt)

      common /testvel2/ ut4,ut5,ut6

      real ymean,xmean

      real vlsum        ! function
      real facevals(lx1,ly1)
      integer intfaces(2**ldim,lelt)
      integer intels(lelt)

      real eps

      real offst        ! Slip velocity across interfaces

      logical ifslip


      n     = nx1*ny1*nz1


      ifslip = .false.

      if (ifslip) then

        if (icalld.eq.0) then

          offst = 0.2

          call opzero(ut1,ut2,ut3)
     
          do i=1,nelv
            ymean = vlsum(ym1(1,1,1,i),n)/n
!            xmean = vlsum(xm1(1,1,1,i),n)/n

            if (ymean.gt.1) then
              flowtype(i)=1
            else
              flowtype(i)=2
            endif
            ymean = flowtype(i) + 0.
            call cadd(ut1(1,1,1,i),ymean,n)
          enddo  
     

          call copy(ut2,ut1,n*nelv)
          call dsavg(ut2)

!          call outpost(ut1,ut2,ut3,pr,t,'   ')


!         Mark the faces we wish to change
          nfaces = 2**ndim

          call izero(intels,nelv)
          call izero(intfaces,nelv*nfaces)

          eps = 1.0e-12

          do i=1,nelv
            if (flowtype(i).eq.1) then    
              do j=1,nfaces
                call facexs(facevals,ut2(1,1,1,i),j,0)
                ymean = vlsum(facevals,lx1*lz1)/(lx1*lz1)
                if (abs((ymean-1.5)).lt.eps) then
                  intels(i) = 1
                  intfaces(j,i) = 1
                  write(6,*) i,j
                endif
              enddo
            endif
          enddo 


          call opzero(ut1,ut2,ut3)
          call rzero(ut3,n*nelv)
          do i=1,nelv
            if (intels(i).eq.1) then
              do j=1,nfaces
                if (intfaces(j,i).eq.1) then
                  write(6,*) i,j
                  call facev(ut1,i,j,offst,nx1,ny1,nz1)
                endif
              enddo  
            endif
          enddo      

          call outpost(ut1,ut2,ut3,pr,t,'slp')

          icalld = icalld + 1

        endif      ! ifcalld.eq.0

      else        ! ifslip == .false.

        call opzero(ut1,ut2,ut3)
        call opzero(ut4,ut5,ut6)    

      endif       ! ifslip

!      call exitt 
      
      
      if (istep.gt.0) then
!        call outpost(ut4,ut5,ut6,usrdiv,t,'ut2')
!        call outpost(vx,vy,vz,pr,t,'   ')
!        call exitt
      endif

      return
      end subroutine slp_mark_faces
c-----------------------------------------------------------------------


