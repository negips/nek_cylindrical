!======================================================================
!     Author: Prabal Negi      
!     Description: Routines for 3D cylindrical solve implementation
!
!====================================================================== 
      subroutine cdabdtp_cyl (ap,wp,h1,h2,h2inv,intype)

!     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
!     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls

      real           ap    (lx2,ly2,lz2,1)
      real           wp    (lx2,ly2,lz2,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      real           h2inv (lx1,ly1,lz1,1)

      real ta1,ta2,ta3,tb1,tb2,tb3
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
     $ ,             tb1 (lx1,ly1,lz1,lelv)
     $ ,             tb2 (lx1,ly1,lz1,lelv)
     $ ,             tb3 (lx1,ly1,lz1,lelv)

      integer intype
      real tolhin,dtbdi

      call opgradt_cyl(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
        tolhin=tolhs
        call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
        if (ifanls) then
          dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
          call opbinv1(tb1,tb2,tb3,ta1,ta2,ta3,dtbdi)
        else
          call opbinv (tb1,tb2,tb3,ta1,ta2,ta3,h2inv)
        endif
      endif
      call opdiv_cyl (ap,tb1,tb2,tb3)

      return
      end
C-----------------------------------------------------------------------

      subroutine opgradt_cyl(outx,outy,outz,inpfld)

!     Compute DTx, DTy, DTz of an input field INPFLD
!     in Cylindrical coordinates

      implicit none

      include 'SIZE'
      include 'GEOM'

      real outx   (lx1,ly1,lz1,1)
      real outy   (lx1,ly1,lz1,1)
      real outz   (lx1,ly1,lz1,1)
      real inpfld (lx2,ly2,lz2,1)

      call cdtp_cyl (outx,inpfld,rxm2,sxm2,txm2,1)
      call cdtp_cyl (outy,inpfld,rym2,sym2,tym2,2)
      if (ldim.eq.3) 
     $   call cdtp_cyl (outz,inpfld,rzm2,szm2,tzm2,3)

      return
      end
!-----------------------------------------------------------------------

      subroutine cdtp_cyl (dtx,x,rm2,sm2,tm2,isd)

!     Compute DT*X (Cylindrical Coordinates)
!     I have assumed all the cross geometric factors are zero.
!     i.e. dr/dy = dr/dz = ds/dx = ds/dz = dt/dx = dt/dy = 0.        
!     I assume the mass matrix already contains a multiplication by R
!     Here I also assume 'R' is the 'y' direction.
!     We can try generalizing some other time. 

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'

      include 'CTIMER'

      real dtx  (lx1*ly1*lz1,lelv)
      real x    (lx2*ly2*lz2,lelv)
      real rm2  (lx2*ly2*lz2,lelv)
      real sm2  (lx2*ly2*lz2,lelv)
      real tm2  (lx2*ly2*lz2,lelv)

      real wx,ta1,ta2,ta3
      common /ctmp1/ wx  (lx1*ly1*lz1)
     $ ,             ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1,lz1)

      REAL           DUAX(LX1)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer isd
      integer e
      integer nxyz1,nxyz2,nyz1,nyz2,nxy1,ly12
      integer n1,n2,i1,i2,iz

#ifdef TIMER
      if (icalld.eq.0) tcdtp=0.0
      icalld=icalld+1
      ncdtp=icalld
      etime1=dnekclock()
#endif

      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2
      nyz2  = ly2*lz2
      nxy1  = lx1*ly1

      n1    = lx1*ly1
      n2    = lx1*ly2

      do e=1,nelv
C       Collocate with weights
        if(ifsplit) then
!         Not implemented
          if (nio.eq.0) then
            write(6,*)
     $        'cdtp_cyl not implemented for Pn-Pn'
            call exitt
          endif   
!          call col3 (wx,bm1(1,1,1,e),x(1,e),nxyz1)
!          call invcol2(wx,jacm1(1,1,1,e),nxyz1)
        else
          call col3 (wx,w3m2,x(1,e),nxyz2)
          if (isd.ne.3) then
            call col2 (wx,ym2(1,1,1,e),nxyz2)
          endif
        endif
C
        if (ldim.eq.2) then

!         Not implemented
          if (nio.eq.0) then
            write(6,*)
     $        'cdtp_cyl not implemented for 2D yet.'
            call exitt
          endif   

!          if (.not.ifdfrm(e) .and. ifalgn(e)) then
!C
!             if (      ifrsxy(e).and.isd.eq.1  .or. 
!     $            .not.ifrsxy(e).and.isd.eq.2) then
!C
!                call col3 (ta1,wx,rm2(1,e),nxyz2)
!                call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
!                call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)
!             else
!                call col3 (ta1,wx,sm2(1,e),nxyz2)
!                call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!                call mxm  (ta2,lx1,dym12,ly2,dtx(1,e),ly1)
!             endif
!          else
!             call col3 (ta1,wx,rm2(1,e),nxyz2)
!             call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
!             call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)
!
!             call col3 (ta1,wx,sm2(1,e),nxyz2)
!             call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!             call mxm  (ta2,lx1,dym12,ly2,ta1,ly1)
!
!             call add2 (dtx(1,e),ta1,nxyz1)
!          endif

        else
          if (ifsplit) then

!             call col3 (ta1,wx,rm2(1,e),nxyz2)
!             call mxm  (dxtm12,lx1,ta1,lx2,dtx(1,e),nyz2)
!             call col3 (ta1,wx,sm2(1,e),nxyz2)
!             i1 = 1
!             i2 = 1
!             do iz=1,lz2
!                call mxm  (ta1(i2),lx1,dym12,ly2,ta2(i1),ly1)
!                i1 = i1 + n1
!                i2 = i2 + n2
!             enddo
!             call add2 (dtx(1,e),ta2,nxyz1)
!             call col3 (ta1,wx,tm2(1,e),nxyz2)
!             call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
!             call add2 (dtx(1,e),ta2,nxyz1)
          else
!           (dv/dr)*(dr/dx_i)*W*p
            call col3 (ta1,wx,rm2(1,e),nxyz2)
            call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
              call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
              i1 = i1 + n1
              i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,izm12,lz2,dtx(1,e),lz1)

!           (dv/ds)*(ds/dx_i)*W*p
            call col3 (ta1,wx,sm2(1,e),nxyz2)
            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
              call mxm  (ta2(i2),lx1,dym12,ly2,ta1(i1),ly1)
              i1 = i1 + n1
              i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,izm12,lz2,ta2,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

!           (dv/dt)*(dt/dx_i)*W*p
            call col3 (ta1,wx,tm2(1,e),nxyz2)
            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
              call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
              i1 = i1 + n1
              i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

!           Additional term in the Radial direction            
!           (v/R)*W*p
            if (isd.eq.2) then
              call invcol3(ta1,wx,ym2(1,1,1,e),nxyz2)
              call col2(ta1,jacm2(1,1,1,e),nxyz2)
              call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
              i1 = 1
              i2 = 1
              do iz=1,lz2
                 call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
                 i1 = i1 + n1
                 i2 = i2 + n2
              enddo
              call mxm  (ta1,nxy1,izm12,lz2,ta2,lz1)
              call add2 (dtx(1,e),ta2,nxyz1)
            endif    ! isd.eq.2

          endif      ! ifsplit
        endif        ! ndim.eq.2 

      enddo

#ifdef TIMER
      tcdtp=tcdtp+(dnekclock()-etime1)
#endif
      return
      end
!---------------------------------------------------------------------- 

      subroutine opdiv_cyl(outfld,inpx,inpy,inpz)

!     Compute OUTFLD = SUMi Di*INPi, 
!     the divergence of the vector field (INPX,INPY,INPZ)


      implicit none

      include 'SIZE'
      include 'GEOM'

      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
      
      integer iflg,ntot2

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call rzero(outfld,ntot2)

      call multd_cyl (work,inpx,rxm2,sxm2,txm2,1,iflg)
      call copy  (outfld,work,ntot2)
      call multd_cyl (work,inpy,rym2,sym2,tym2,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
        call multd_cyl (work,inpz,rzm2,szm2,tzm2,3,iflg)
        call add2  (outfld,work,ntot2)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine multd_cyl (du,u,rm2,sm2,tm2,isd,iflg)

!     Compute D*X
!     X    : input variable, defined on M1
!     DX   : output variable, defined on M2 (note: D is rectangular)   
!     RM2 : RXM2, RYM2 or RZM2
!     SM2 : SXM2, SYM2 or SZM2
!     TM2 : TXM2, TYM2 or TZM2
!     ISD : spatial direction (x=1,y=2,z=3)
!     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'
      include 'CTIMER'

      real           du   (lx2*ly2*lz2,lelv)
      real           u    (lx1*ly1*lz1,lelv)
      real           rm2  (lx2*ly2*lz2,lelv)
      real           sm2  (lx2*ly2*lz2,lelv)
      real           tm2  (lx2*ly2*lz2,lelv)

      integer isd,iflg

      real ta1,ta2,ta3
      common /ctmp1/ ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1*lz1)

      real           duax(lx1)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      integer e,i1,i2,iz
      integer nxy1,nyz1,nxy2,nxyz1,nxyz2,n1,n2


#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      nyz1  = ly1*lz1
      nxy2  = lx2*ly2
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      n1    = lx2*ly1
      n2    = lx2*ly2

      if (ndim.eq.2) then
!       Not implemented
        if (nio.eq.0) then
          write(6,*)
     $      'multd_cyl not implemented for 2D yet.'
          call exitt
        endif
      endif 

      if (ifsplit) then
!       Not implemented
        if (nio.eq.0) then
          write(6,*)
     $      'multd_cyl not implemented for Pn-Pn'
          call exitt
        endif
      endif 

      do e=1,nelv

!       du/dr
        call mxm (dxm12,lx2,u(1,e),lx1,ta1,nyz1)
        i1=1
        i2=1
        do iz=1,lz1
          call mxm (ta1(i1),lx2,iytm12,ly1,ta2(i2),ly2)
          i1=i1+n1
          i2=i2+n2
        enddo
        call mxm  (ta2,nxy2,iztm12,lz1,du(1,e),lz2)
!       dr/dx_i*du/dr        
        call col2 (du(1,e),rm2(1,e),nxyz2)

!       du/ds        
        call mxm  (ixm12,lx2,u(1,e),lx1,ta3,nyz1)
        i1=1
        i2=1
        do iz=1,lz1
          call mxm (ta3(i1),lx2,dytm12,ly1,ta2(i2),ly2)
          i1=i1+n1
          i2=i2+n2
        enddo
        call mxm     (ta2,nxy2,iztm12,lz1,ta1,lz2)
!       ds/dx_i*du/ds        
        call addcol3 (du(1,e),ta1,sm2(1,e),nxyz2)

!       du/dt        
        call mxm  (ixm12,lx2,u(1,e),lx1,ta3,nyz1)
        i1=1
        i2=1
        do iz=1,lz1
          call mxm (ta3(i1),lx2,iytm12,ly1,ta2(i2),ly2)
          i1=i1+n1
          i2=i2+n2
        enddo
        call mxm (ta2,nxy2,dztm12,lz1,ta3,lz2)
!       dt/dx_i*du/dt        
        call addcol3 (du(1,e),ta3,tm2(1,e),nxyz2)

!       Collocate with the weights and Radius on the pressure mesh
        call col2 (du(1,e),w3m2,nxyz2)
        if (isd.ne.3) then
          call col2 (du(1,e),ym2(1,1,1,e),nxyz2)
        endif

!       Add additional Radial term
        if (isd.eq.2) then
!         I12*u        
          call mxm (ixm12,lx2,u(1,e),lx1,ta1,nyz1)
          i1=1
          i2=1
          do iz=1,lz1
            call mxm (ta1(i1),lx2,iytm12,ly1,ta2(i2),ly2)
            i1=i1+n1
            i2=i2+n2
          enddo
          call mxm  (ta2,nxy2,iztm12,lz1,ta3,lz2)

!         W*I12*u          
          call col3 (ta1,w3m2,ta3,nxyz2)
!         J*W*I12*u          
          call col2 (ta1,jacm2(1,1,1,e),nxyz2)
          call add2 (du(1,e),ta1,nxyz2)
        endif
      enddo       ! e=1,nelv

#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      END
c-----------------------------------------------------------------------





