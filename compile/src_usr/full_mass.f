!======================================================================
!     Author: Prabal Negi      
!     Description: Evaluations using full mass matrices
!     Routines:   setup_interp_fm     : lx1 -> lxfm setup
!                 intp_fm()           : Interpolate lx1 -> lxfm      
!
!====================================================================== 

      subroutine setup_fm()

      implicit none

!     Set up weights/interpolator from lx1 to lxfm mesh 

      include 'SIZE'
      include 'WZ'

      include 'FULLMASS'

      integer icalld
      save icalld
      data icalld /0/

!     If already initialized      
      if (icalld.gt.0) return

!     Gauss-Legendre Mesh 
      call zwgll (fm_z,fm_wght,lxfm)

!     Interpolator from M1 mesh to lxfm Mesh      
      call igllm (fm_jgl,fm_jglt,zgm1(1,1),fm_z,lx1,lxfm,lx1,lxfm)

      icalld = icalld+1

      return
      end
c-----------------------------------------------------------------------

      subroutine intp_fm(fldf,fld)

      implicit none

      include 'SIZE'
      include 'GEOM'

      include 'FULLMASS'

      real fld  (lx1*ly1*lz1)
      real fldf (lxfm**ldim)

      integer iz,i1,i2


      if (ndim.eq.2) then
        call mxm (fm_jgl,lxfm,fld,lx1,fm_wk1,ly1)
        call mxm (fm_wk1,lxfm,fm_jglt,ly1,fldf,lxfm)
      else        
        call mxm (fm_jgl,lxfm,fld,lx1,fm_wk1,ly1*lz1)
        i1=1
        i2=1
        do iz=1,lz1
          call mxm (fm_wk1(i1),lxfm,fm_jglt,ly1,fm_wk2(i2),lxfm)
          i1=i1+lxfm*ly1
          i2=i2+lxfm*lxfm
        enddo
        call mxm  (fm_wk2,lxfm*lxfm,fm_jglt,lz1,fldf,lxfm)
      endif        


      return
      end subroutine intp_fm
!---------------------------------------------------------------------- 

      function vlsc2_fm(u,v)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      include 'FULLMASS'

      real vlsc2_fm
      real s
      real u(lx1*ly1*lz1,lelt)
      real v(lx1*ly1*lz1,lelt)

      integer e,i,j,k,ii
      integer lxfm3

      real vlsum

      call setup_fm()

      lxfm3 = lxfm**ndim

      s = 0.0
      do e=1,nelv

!       u_lx1 -> u_lxfm          
        call intp_fm(fm_wk3,u(1,e)) 

!       v_lx1 -> v_lxfm
        call intp_fm(fm_wk4,v(1,e))
        call col2(fm_wk3,fm_wk4,lxfm3)

!       J_lx1 -> J_lxfm
        call intp_fm(fm_wk4,jacm1(1,1,1,e))
        call col2(fm_wk3,fm_wk4,lxfm3)

!       W*J*\rho*u
        if (ndim.eq.2) then
          do j=1,lxfm
          do i=1,lxfm
            ii = i + (j-1)*lxfm
            fm_wk3(ii) = fm_wk3(ii)*fm_wght(i)*fm_wght(j)
          enddo
          enddo
        else          
          do k=1,lxfm
          do j=1,lxfm
          do i=1,lxfm
            ii = i + (j-1)*lxfm + (k-1)*lxfm*lxfm
            fm_wk3(ii) = fm_wk3(ii)*fm_wght(i)*fm_wght(j)*fm_wght(k)
          enddo
          enddo
          enddo
        endif  

        s = s + vlsum(fm_wk3,lxfm3) 

      enddo
      
      vlsc2_fm = s


      return
      end function vlsc2_fm
!---------------------------------------------------------------------- 

      function glsc2_fm(u,v)

      implicit none

      include 'SIZE'

      real u(1),v(1)
      real sc,tmp
      real glsc2_fm
      real vlsc2_fm

      sc = vlsc2_fm(u,v)
      call gop(sc,tmp,'+  ',1)

      glsc2_fm = sc

      return 
      end function glsc2_fm 

!---------------------------------------------------------------------- 



