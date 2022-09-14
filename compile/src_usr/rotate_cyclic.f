!====================================================================== 
!
!     Author: Prabal Negi     
!     Description: Rotations for Cyclic Boundaries
!
!
!======================================================================       

      subroutine rotate_cyc_all(r1,r2,r3,idir)

!     Rotate in the x-y, y-z or z-x planes.

      implicit none  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'

      real r1(lx1,ly1,lz1,1)
     $   , r2(lx1,ly1,lz1,1)
     $   , r3(lx1,ly1,lz1,1)

      integer idir
      integer e,f,nface
      logical ifxy,ifyz,ifzx
      real length, dotprod
      integer j1,j2,js1,jf1,jskip1,js2,jf2,jskip2
      integer k
      real tol
      real cost,sint,rnor,rtn1


      call rotate_cyc_yz(r1,r2,r3,idir)
      return
 
c     (1) Face n-t transformation

!     In principle this can be made more general so that the periodic
!     face can be arbitrarily aligned in space. One needs a consistent
!     method of rotation on all the faces.
      ifxy = .false.
      ifyz = .false.
      ifzx = .false.

      tol  = 1.0e-12;
      nface = 2*ldim
      do e=1,nelfld(ifield)
      do f=1,nface

         if(cbc(f,e,ifield) .eq. 'P  '.or.cbc(f,e,ifield).eq.'p  ')then

            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,f)
            if (idir.eq.1) then
              k=0
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1
                k=k+1

                if (abs(unx(k,1,f,e)).lt.tol) ifyz = .true.
                if (abs(uny(k,1,f,e)).lt.tol) ifzx = .true.
                if (abs(unz(k,1,f,e)).lt.tol) ifxy = .true.

                if (ifxy) then
                  length = unx(k,1,f,e)**2 + uny(k,1,f,e)**2
                  length = sqrt(length)
                 
                  dotprod = unx(k,1,f,e)*ym1(j1,j2,1,e)
     $                     -uny(k,1,f,e)*xm1(j1,j2,1,e)
                  cost =  unx(k,1,f,e)/length
                  sint =  uny(k,1,f,e)/length
                  rnor = ( r1(j1,j2,1,e)*cost + r2(j1,j2,1,e)*sint )
                  rtn1 = (-r1(j1,j2,1,e)*sint + r2(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r1(j1,j2,1,e) = rnor
                     r2(j1,j2,1,e) = rtn1
                  else
                     r1(j1,j2,1,e) =-rnor
                     r2(j1,j2,1,e) =-rtn1
                  endif
                elseif (ifyz) then  
                  length = uny(k,1,f,e)**2 + unz(k,1,f,e)**2
                  length = sqrt(length)

                  dotprod = uny(k,1,f,e)*zm1(j1,j2,1,e)
     $                     -unz(k,1,f,e)*ym1(j1,j2,1,e)
                  cost =  uny(k,1,f,e)/length
                  sint =  unz(k,1,f,e)/length
                  rnor = ( r2(j1,j2,1,e)*cost + r3(j1,j2,1,e)*sint )
                  rtn1 = (-r2(j1,j2,1,e)*sint + r3(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r2(j1,j2,1,e) = rnor
                     r3(j1,j2,1,e) = rtn1
                  else
                     r2(j1,j2,1,e) =-rnor
                     r3(j1,j2,1,e) =-rtn1
                  endif
                elseif (ifzx) then
                  length = unz(k,1,f,e)**2 + unx(k,1,f,e)**2
                  length = sqrt(length)
                 
                  dotprod = unz(k,1,f,e)*xm1(j1,j2,1,e)
     $                     -unx(k,1,f,e)*zm1(j1,j2,1,e)
                  cost =  unz(k,1,f,e)/length
                  sint =  unx(k,1,f,e)/length
                  rnor = ( r3(j1,j2,1,e)*cost + r1(j1,j2,1,e)*sint )
                  rtn1 = (-r3(j1,j2,1,e)*sint + r1(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r3(j1,j2,1,e) = rnor
                     r1(j1,j2,1,e) = rtn1
                  else
                     r3(j1,j2,1,e) =-rnor
                     r1(j1,j2,1,e) =-rtn1
                  endif
                endif    ! ifxy...
              enddo     ! j2
              enddo     ! j1

            else    ! reverse rotate

              k=0
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1
                k=k+1

                if (abs(unx(k,1,f,e)).lt.tol) ifyz = .true.
                if (abs(uny(k,1,f,e)).lt.tol) ifzx = .true.
                if (abs(unz(k,1,f,e)).lt.tol) ifxy = .true.

                if (ifxy) then
                  length = unx(k,1,f,e)**2 + uny(k,1,f,e)**2
                  length = sqrt(length)
                 
                  dotprod = unx(k,1,f,e)*ym1(j1,j2,1,e)
     $                     -uny(k,1,f,e)*xm1(j1,j2,1,e)
                  cost =  unx(k,1,f,e)/length
                  sint =  uny(k,1,f,e)/length
                  rnor = (r1(j1,j2,1,e)*cost - r2(j1,j2,1,e)*sint )
                  rtn1 = (r1(j1,j2,1,e)*sint + r2(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r1(j1,j2,1,e) = rnor
                     r2(j1,j2,1,e) = rtn1
                  else
                     r1(j1,j2,1,e) =-rnor
                     r2(j1,j2,1,e) =-rtn1
                  endif
                elseif (ifyz) then  
                  length = uny(k,1,f,e)**2 + unz(k,1,f,e)**2
                  length = sqrt(length)

                  dotprod = uny(k,1,f,e)*zm1(j1,j2,1,e)
     $                     -unz(k,1,f,e)*ym1(j1,j2,1,e)
                  cost =  uny(k,1,f,e)/length
                  sint =  unz(k,1,f,e)/length
                  rnor = (r2(j1,j2,1,e)*cost - r3(j1,j2,1,e)*sint )
                  rtn1 = (r2(j1,j2,1,e)*sint + r3(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r2(j1,j2,1,e) = rnor
                     r3(j1,j2,1,e) = rtn1
                  else
                     r2(j1,j2,1,e) =-rnor
                     r3(j1,j2,1,e) =-rtn1
                  endif
                elseif (ifzx) then
                  length = unz(k,1,f,e)**2 + unx(k,1,f,e)**2
                  length = sqrt(length)
                 
                  dotprod = unz(k,1,f,e)*xm1(j1,j2,1,e)
     $                     -unx(k,1,f,e)*zm1(j1,j2,1,e)
                  cost =  unz(k,1,f,e)/length
                  sint =  unx(k,1,f,e)/length
                  rnor = (r3(j1,j2,1,e)*cost - r1(j1,j2,1,e)*sint )
                  rtn1 = (r3(j1,j2,1,e)*sint + r1(j1,j2,1,e)*cost )
                  if (dotprod .ge. 0.0) then 
                     r3(j1,j2,1,e) = rnor
                     r1(j1,j2,1,e) = rtn1
                  else
                     r3(j1,j2,1,e) =-rnor
                     r1(j1,j2,1,e) =-rtn1
                  endif
                endif  
              enddo           ! j2
              enddo           ! j1
            endif             ! idir
          endif               ! if cb.eq.'P  '
      enddo                   ! f
      enddo                   ! e

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_cyc_xy(r1,r2,r3,idir)

!     Rotate in the x-y plane.
!     Angles are based on x,y locations

      implicit none  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'         ! Should have pi

      real r1(lx1,ly1,lz1,1)
     $   , r2(lx1,ly1,lz1,1)
     $   , r3(lx1,ly1,lz1,1)

      integer idir
      integer e,f,nface
      logical ifxy,ifyz,ifzx
      real length, dotprod
      integer j1,j2,js1,jf1,jskip1,js2,jf2,jskip2
      integer k
      real tol
      real theta,phi
      real cost,sint,rnor,rtn1
      real x,y,z
      real ux,uy,uz


!     In principle this can be made more general so that the periodic
!     face can be arbitrarily aligned in space. One needs a consistent
!     method of rotation on all the faces.
      ifxy = .true.
      ifyz = .false.
      ifzx = .false.

      tol  = 1.0e-08;
      nface = 2*ldim
      do e=1,nelfld(ifield)
      do f=1,nface

         if(cbc(f,e,ifield) .eq. 'P  '.or.cbc(f,e,ifield).eq.'p  ')then

            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,f)
            if (idir.eq.1) then
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1

                x=xm1(j1,j2,1,e)      
                y=ym1(j1,j2,1,e)      
                z=zm1(j1,j2,1,e)

                phi     = atan2(-x,y)     ! Angle of the normal
                theta   = atan2(y,x)      ! Angle of the tangent

                ux  = r1(j1,j2,1,e)
                uy  = r2(j1,j2,1,e)

                call rotate2d(rnor,rtn1,ux,uy,phi)

                r1(j1,j2,1,e) = rnor
                r2(j1,j2,1,e) = rtn1

              enddo     ! j2
              enddo     ! j1

            else    ! reverse rotate
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1

                x=xm1(j1,j2,1,e)      
                y=ym1(j1,j2,1,e)      
                z=zm1(j1,j2,1,e)

                phi     = atan2(-x,y)     ! Angle of the normal
                theta   = atan2(y,x)      ! Angle of the tangent

                ux  = r1(j1,j2,1,e)
                uy  = r2(j1,j2,1,e)

                call rotate2d(rnor,rtn1,ux,uy,-phi)

                r1(j1,j2,1,e) = rnor
                r2(j1,j2,1,e) = rtn1

              enddo           ! j2
              enddo           ! j1
            endif             ! idir
          endif               ! if cb.eq.'P  '
      enddo                   ! f
      enddo                   ! e

      return
      end subroutine rotate_cyc_xy

!---------------------------------------------------------------------- 

      subroutine rotate_cyc_yz(r1,r2,r3,idir)

!     Rotate in the y-z plane.
!     Angles are based on y,z locations

      implicit none  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'         ! Should have pi

      real r1(lx1,ly1,lz1,1)
     $   , r2(lx1,ly1,lz1,1)
     $   , r3(lx1,ly1,lz1,1)

      integer idir
      integer e,f,nface
      logical ifxy,ifyz,ifzx
      real length, dotprod
      integer j1,j2,js1,jf1,jskip1,js2,jf2,jskip2
      integer k
      real tol
      real theta,phi
      real cost,sint,rnor,rtn1
      real x,y,z
      real ux,uy,uz


!     In principle this can be made more general so that the periodic
!     face can be arbitrarily aligned in space. One needs a consistent
!     method of rotation on all the faces.
      ifxy = .false.
      ifyz = .true.
      ifzx = .false.

      tol  = 1.0e-08;
      nface = 2*ldim
      do e=1,nelfld(ifield)
      do f=1,nface

         if(cbc(f,e,ifield) .eq. 'P  '.or.cbc(f,e,ifield).eq.'p  ')then

            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,f)
            if (idir.eq.1) then
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1

                x=xm1(j1,j2,1,e)      
                y=ym1(j1,j2,1,e)      
                z=zm1(j1,j2,1,e)

                phi     = atan2(y,z)     ! Angle of the normal
                theta   = atan2(z,y)      ! Angle of the tangent

                uy  = r2(j1,j2,1,e)
                uz  = r3(j1,j2,1,e)

                call rotate2d(rtn1,rnor,uy,uz,theta)

                r2(j1,j2,1,e) = rtn1
                r3(j1,j2,1,e) = rnor

              enddo     ! j2
              enddo     ! j1

            else    ! reverse rotate
              do j2=js2,jf2,jskip2
              do j1=js1,jf1,jskip1

                x=xm1(j1,j2,1,e)      
                y=ym1(j1,j2,1,e)      
                z=zm1(j1,j2,1,e)

                phi     = atan2(y,z)     ! Angle of the normal
                theta   = atan2(z,y)      ! Angle of the tangent

                uy  = r2(j1,j2,1,e)
                uz  = r3(j1,j2,1,e)

                call rotate2d(rtn1,rnor,uy,uz,-theta)

                r2(j1,j2,1,e) = rtn1
                r3(j1,j2,1,e) = rnor

              enddo           ! j2
              enddo           ! j1
            endif             ! idir
          endif               ! if cb.eq.'P  '
      enddo                   ! f
      enddo                   ! e

      return
      end subroutine rotate_cyc_yz

!---------------------------------------------------------------------- 

      subroutine rotate2d(rx,ry,rx0,ry0,th)

      implicit none

      real rx,ry,rx0,ry0
      real th        ! theta, assumed to be in radians

      rx = rx0*cos(th) + ry0*sin(th)
      ry = -rx0*sin(th) + ry0*cos(th)
      

      return
      end subroutine rotate2d
!---------------------------------------------------------------------- 





