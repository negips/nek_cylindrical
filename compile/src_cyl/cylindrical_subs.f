!======================================================================
!     Author: Prabal Negi      
!     Description: Routines for 3D cylindrical solve implementation
!     Routines:   init_cyl()              : Cylindrical Coordinates initialization
!                 reset_geom_cyl()        : Reset MASS/Geometric factors 
!                 cdabdtp_cyl()           : Pressure Pseudolaplacian
!                 opdiv_cyl()             : Cylindrical Divergence
!                 multd_cyl()             : D*u = q*(D*u)
!                 multd_rho_cyl()         : \rho*D*u = \rho*q*(D*u)
!                 opgradt_cyl()           : Pressure gradient term
!                 cdtp_cyl()              : (D^T)*p = p*(D*v)
!                 convect_cylindircal_rho : Cylindrical convective term with density
!                 convect_cylindircal     : Cylindrical convective term
!                 dealias_rho_uv          : Dealiased \rho*u*v      
!                 dealias_uv              : Dealiased u*v
!                 axhmsf_cyl              : Ax (cylindrical and coupled)
!                 stnrate_cyl             : 1/2*(grad(u) + grad(u)^T)      
!                 stress_cyl              : Sij = 2*mu*Eij
!                 div_stress_cyl          : Grad(v)\dot Sij
!                 ttxyz_cyl               : Grad(v)\dot      
!
!======================================================================
      subroutine init_cyl()
      
      implicit none

      include 'SIZE'
      include 'GEOM'          ! ym1,ym2
      include 'MASS'          ! volvm1
      include 'TSTEP'         ! ifield
      include 'SOLN'          ! vdiff,vtrans
      include 'INPUT'

      include 'CYLINDRICAL'

      integer n,n2

      real glmin,glmax,glsc2
      real r1,r2,omg1,omg2,d

      integer ifld

      real nufld
      common /scrcg/ nufld(lx1,ly1,lz1,lelt)

      real nuavg
      character*132 str

      n = lx1*ly1*lz1*nelv
      call copy(cyl_radius,ym1,n)       ! Radial coordinate (M1)

      n2 = lx2*ly2*lz2*nelv
      call copy(cyl_radius2,ym2,n2)     ! Radial coordinate (M2)

      cyl_rady(1) = glmin(cyl_radius,n)
      cyl_rady(2) = glmax(cyl_radius,n)

!     Negative value means we specify \Omega1*R1
      if (cyl_omega(1).lt.0) then
        cyl_omega(1) = -cyl_omega(1)/cyl_rady(1)
      endif 

!     Negative value means we specify \Omega2*R2
      if (cyl_omega(1).lt.0) then
        cyl_omega(2) = -cyl_omega(2)/cyl_rady(2)
      endif  

!     Laminar Taylor-Couette Parameters
!----------------------------------------       

!     Laminar Taylor Couette solution
      omg1        = cyl_omega(1)
      omg2        = cyl_omega(2)
      r1          = cyl_rady(1)
      r2          = cyl_rady(2)
      d           = r2-r1

      cyl_tc_ab(1) = (omg2*r2*r2 - omg1*r1*r1)/(r2**2 - r1**2)
      cyl_tc_ab(2) = (omg1 - omg2)*(r1**2)*(r2**2)/(r2**2 - r1**2)

!     For the Narrow gap approximation
!     nu = mu/rho 
      call invcol3(nufld,vdiff,vtrans,n)
      nuavg = glsc2(nufld,bm1,n)/volvm1

      cyl_ta_num   = 4.0*(omg1*r1*r1 - omg2*r2*r2)*(omg1*(d**4)) 
     $                  /((r2**2 - r1**2)*(nuavg**2))

      call blank(str,132)
      write(str,'(A9,1x,E12.4E2)') 'Radius 1:', r1
      call mntr_log(cyl_id,cyl_log,str)

      call blank(str,132)
      write(str,'(A9,1x,E12.4E2)') 'Radius 2:', r2
      call mntr_log(cyl_id,cyl_log,str)

      call blank(str,132)
      write(str,'(A16,1x,E12.4E2)') 'Taylor Coutte A:', cyl_tc_ab(1)
      call mntr_log(cyl_id,cyl_log,str)

      call blank(str,132)
      write(str,'(A16,1x,E12.4E2)') 'Taylor Coutte B:', cyl_tc_ab(2)
      call mntr_log(cyl_id,cyl_log,str)

      call blank(str,132)
      write(str,'(A14,1x,E16.8E2)') 'Taylor Number:', cyl_ta_num
      call mntr_log(cyl_id,cyl_log,str)
!---------------------------------------- 

!     Add the factor of Radius into mass matrix and
!     geometric factors
      call reset_geom_cyl()

!     Preconditioner
      param(42)=0       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=1       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=0       ! 0: E based Schwartz (SEM), 1: A based Schwartz (FEM)


      return
      end subroutine init_cyl
!----------------------------------------------------------------------

      subroutine reset_geom_cyl()

      implicit none

      include 'SIZE'
      include 'GEOM'          ! GiM1
      include 'MASS'
      include 'TSTEP'         ! ifield
      include 'INPUT'
      include 'WZ'            ! W3M1

      include 'CYLINDRICAL'

      integer e,i,n,n2,nxyz1
      real r,ri,wght,jaci

      integer ifld

!     For cylindrical solver
!     I probably need to do this in at every time step at igeom=2
!     After the geometry has been regenerated
      n = lx1*ly1*lz1*nelv
      call copy(cyl_radius,ym1,n)       ! Radial coordinate (M1)

      n2 = lx2*ly2*lz2*nelv
      call copy(cyl_radius2,ym2,n2)     ! Radial coordinate (M2)

!     BM2*R
      call col2(bm2,cyl_radius2,n2)
      call invers2(bm2inv,bm2,n2)

      ifld   = ifield
      ifield = 1
!     BM1*R      
      call col2    (bm1,cyl_radius,n)
      call copy    (binvm1,bm1,n)
      call dssum   (binvm1,lx1,ly1,lz1)
      call invcol1 (binvm1,n)

!     Geometric factors for the integrated del-squared operator
      nxyz1 = lx1*ly1*lz1      
      if (ldim.eq.2) then
        do e=1,nelv
          do i=1,nxyz1
            r             = cyl_radius(i,1,1,e)       ! Radius
            jaci          = 1.0/jacm1(i,1,1,e)        ! 1/Jac
            wght          = w3m1(i,1,1)               ! W

!           G11          
            g1m1(i,1,1,e) =  r*jaci*wght*(rxm1(i,1,1,e)*rxm1(i,1,1,e)
     $                        +rym1(i,1,1,e)*rym1(i,1,1,e))

!           G22          
            g2m1(i,1,1,e) =  r*jaci*wght*(sxm1(i,1,1,e)*sxm1(i,1,1,e)
     $                        +sym1(i,1,1,e)*sym1(i,1,1,e))

!           G12          
            g4m1(i,1,1,e) =  r*jaci*wght*(rxm1(i,1,1,e)*sxm1(i,1,1,e)
     $                        +rym1(i,1,1,e)*sxm1(i,1,1,e))
          enddo
        enddo  
      else
        do e=1,nelv
          do i=1,nxyz1
            r             = cyl_radius(i,1,1,e)       ! Radius
            ri            = 1.0/cyl_radius(i,1,1,e)   ! 1/R
            jaci          = 1.0/jacm1(i,1,1,e)        ! 1/Jac
            wght          = w3m1(i,1,1)               ! W

!           G11          
            g1m1(i,1,1,e) =   r*jaci*wght*(rxm1(i,1,1,e)*rxm1(i,1,1,e))
     $                     +  r*jaci*wght*(rym1(i,1,1,e)*rym1(i,1,1,e))
     $                     + ri*jaci*wght*(rzm1(i,1,1,e)*rzm1(i,1,1,e))

!           G22          
            g2m1(i,1,1,e) =   r*jaci*wght*(sxm1(i,1,1,e)*sxm1(i,1,1,e))
     $                     +  r*jaci*wght*(sym1(i,1,1,e)*sym1(i,1,1,e))
     $                     + ri*jaci*wght*(szm1(i,1,1,e)*szm1(i,1,1,e)) 

!           G33 
            g3m1(i,1,1,e) =   r*jaci*wght*(txm1(i,1,1,e)*txm1(i,1,1,e))
     $                     +  r*jaci*wght*(tym1(i,1,1,e)*tym1(i,1,1,e))
     $                     + ri*jaci*wght*(tzm1(i,1,1,e)*tzm1(i,1,1,e)) 

!           G12 
            g4m1(i,1,1,e) =   r*jaci*wght*(rxm1(i,1,1,e)*sxm1(i,1,1,e))
     $                     +  r*jaci*wght*(rym1(i,1,1,e)*sym1(i,1,1,e))
     $                     + ri*jaci*wght*(rzm1(i,1,1,e)*szm1(i,1,1,e))

!           G13
            g5m1(i,1,1,e) =   r*jaci*wght*(rxm1(i,1,1,e)*txm1(i,1,1,e))
     $                     +  r*jaci*wght*(rym1(i,1,1,e)*tym1(i,1,1,e))
     $                     + ri*jaci*wght*(rzm1(i,1,1,e)*tzm1(i,1,1,e))

!           G23
            g6m1(i,1,1,e) =   r*jaci*wght*(sxm1(i,1,1,e)*txm1(i,1,1,e))
     $                     +  r*jaci*wght*(sym1(i,1,1,e)*tym1(i,1,1,e))
     $                     + ri*jaci*wght*(szm1(i,1,1,e)*tzm1(i,1,1,e))

          enddo   ! i=1,nxyz1
        enddo     ! e=1,nelv
      endif

      if (ifheat) then ! temperature mass matrix
        ifield = 2
        n      = lx1*ly1*lz1*nelt
        call copy    (bintm1,bm1,n)
        call dssum   (bintm1,lx1,ly1,lz1)
        call invcol1 (bintm1,n)
      endif

      ifield = ifld

      return 
      end subroutine reset_geom_cyl
!----------------------------------------------------------------------       
      subroutine cdabdtp_cyl (ap,wp,h1,h2,h2inv,intype)

!     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
!     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls
      include 'SOLN'          ! V?MASK

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

      logical iffullmass
      integer maxiter

      iffullmass = .false.

      call opgradt_cyl(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
        tolhin=tolhs
        call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
        if (ifanls) then
          dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
          call opbinv1(tb1,tb2,tb3,ta1,ta2,ta3,dtbdi)
        else
          if (iffullmass) then
            maxiter = 100
            call opcopy(tb1,tb2,tb3,ta1,ta2,ta3)
            call op_gmres(tb1,tb2,tb3,h2,maxiter)
!            call cyl_gmres(tb1,h2,v1mask,maxiter)
!            call cyl_gmres(tb2,h2,v2mask,maxiter)
!            call cyl_gmres(tb3,h2,v3mask,maxiter)
          else
            call opbinv (tb1,tb2,tb3,ta1,ta2,ta3,h2inv)
          endif  
        endif
      endif
      call opdiv_cyl (ap,tb1,tb2,tb3)
!      call opdiv_rho_cyl (ap,tb1,tb2,tb3)

      return
      end
C-----------------------------------------------------------------------

      subroutine opgradt_cyl(outx,outy,outz,inpfld)

!     Compute DTx, DTy, DTz of an input field INPFLD
!     in Cylindrical coordinates

      implicit none

      include 'SIZE'
      include 'GEOM'

      real outx   (lx1,ly1,lz1,lelv)
      real outy   (lx1,ly1,lz1,lelv)
      real outz   (lx1,ly1,lz1,lelv)
      real inpfld (lx2,ly2,lz2,lelv)

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
            call tensor3_op(dtx(1,e),ta1,lx2,ly2,lz2,
     $                      dxtm12,iym12,izm12,lx1,ly1,lz1)

!           (dv/ds)*(ds/dx_i)*W*p
            call col3 (ta1,wx,sm2(1,e),nxyz2)
            call tensor3_op(ta2,ta1,lx2,ly2,lz2,
     $                      ixtm12,dym12,izm12,lx1,ly1,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

!           (dv/dt)*(dt/dx_i)*W*p
            call col3 (ta1,wx,tm2(1,e),nxyz2)
            call tensor3_op(ta2,ta1,lx2,ly2,lz2,
     $                      ixtm12,iym12,dzm12,lx1,ly1,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

!           Additional term in the Radial direction            
!           (v/R)*W*p
            if (isd.eq.2) then
              call invcol3(ta1,wx,ym2(1,1,1,e),nxyz2)
              call col2(ta1,jacm2(1,1,1,e),nxyz2)
              call tensor3_op(ta2,ta1,lx2,ly2,lz2,
     $                        ixtm12,iym12,izm12,lx1,ly1,lz1)
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

      real outfld (lx2,ly2,lz2,lelv)
      real inpx   (lx1,ly1,lz1,lelv)
      real inpy   (lx1,ly1,lz1,lelv)
      real inpz   (lx1,ly1,lz1,lelv)

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
      subroutine opdiv_rho_cyl(outfld,inpx,inpy,inpz)

!     Compute OUTFLD = SUMi Di*INPi, 
!     the divergence of the vector field (INPX,INPY,INPZ)


      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'    ! vtrans

      real outfld (lx2,ly2,lz2,lelv)
      real inpx   (lx1,ly1,lz1,lelv)
      real inpy   (lx1,ly1,lz1,lelv)
      real inpz   (lx1,ly1,lz1,lelv)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
      
      integer iflg,ntot2

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call rzero(outfld,ntot2)

      call multd_rho_cyl (work,inpx,vtrans,rxm2,sxm2,txm2,1,iflg)
      call copy  (outfld,work,ntot2)
      call multd_rho_cyl (work,inpy,vtrans,rym2,sym2,tym2,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
        call multd_rho_cyl (work,inpz,vtrans,rzm2,szm2,tzm2,3,iflg)
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
      end
c-----------------------------------------------------------------------

      subroutine multd_rho_cyl (du,u,rho,rm2,sm2,tm2,isd,iflg)

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
      real           rho  (lx1*ly1*lz1,lelv)    ! Density
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
        call tensor3_op(du(1,e),u(1,e),lx1,ly1,lz1,dxm12,iytm12,
     $                  iztm12,lx2,ly2,lz2)
!       dr/dx_i*du/dr        
        call col2 (du(1,e),rm2(1,e),nxyz2)

!       du/ds
        call tensor3_op(ta1,u(1,e),lx1,ly1,lz1,ixm12,dytm12,
     $                  iztm12,lx2,ly2,lz2)
!       ds/dx_i*du/ds        
        call addcol3 (du(1,e),ta1,sm2(1,e),nxyz2)

!       du/dt
        call tensor3_op(ta3,u(1,e),lx1,ly1,lz1,ixm12,iytm12,
     $                  dztm12,lx2,ly2,lz2)
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
!         ta3 = u (M1 -> M2)
          call tensor3_op(ta3,u(1,e),lx1,ly1,lz1,ixm12,iytm12,
     $                    iztm12,lx2,ly2,lz2)

!         W*I12*u          
          call col3 (ta1,w3m2,ta3,nxyz2)
!         J*W*I12*u          
          call col2 (ta1,jacm2(1,1,1,e),nxyz2)
          call add2 (du(1,e),ta1,nxyz2)
        endif

!       ta3 = \rho (M1 -> M2)
        call tensor3_op(ta3,rho(1,e),lx1,ly1,lz1,ixm12,iytm12,
     $                  iztm12,lx2,ly2,lz2)
        
!       Collocate with density
        call col2(du,ta3,nxyz2) 
       
      enddo       ! e=1,nelv

#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      end subroutine multd_rho_cyl
c-----------------------------------------------------------------------
      subroutine convect_cylindrical_rho(bdu,rho,u,cx,cy,cz)

!     v*\rho*U \dot \grad(U)        
!     Compute dealiased form:  (J^T)*(J*rho)*(Bf)*(J*C) .Grad Ju w/ correct Jacobians
!     Grad is the cylindrical grad operator

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT' 

      real bdu(1),u(1),rho(1),cx(1),cy(1),cz(1)

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e,i
      integer iu,ic,ib,ir
      integer nxyzu,nxyzc,nxyz1,nxyzd,nxyzr

      real radf(ltd)
      real rhof(ltd)

!     Geometric factors Mapping.
!     3D:      
!     rxm1(1,1,1,e) -->  rx(1,1,e)
!     rym1(1,1,1,e) -->  rx(1,2,e)
!     rzm1(1,1,1,e) -->  rx(1,3,e)
!     sxm1(1,1,1,e) -->  rx(1,4,e)
!     sym1(1,1,1,e) -->  rx(1,5,e)
!     szm1(1,1,1,e) -->  rx(1,6,e)
!     txm1(1,1,1,e) -->  rx(1,7,e)
!     tym1(1,1,1,e) -->  rx(1,8,e)
!     tzm1(1,1,1,e) -->  rx(1,9,e)

!     2D:
!     rxm1(1,1,1,e) -->  rx(1,1,e)
!     rym1(1,1,1,e) -->  rx(1,2,e)
!     sxm1(1,1,1,e) -->  rx(1,3,e)
!     sym1(1,1,1,e) -->  rx(1,4,e)

   
   
!     Put the geometric factors on the fine mesh
!     Also includes multiplication by Weights      
      call set_dealias_rx

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
!      if (ifcf) nxyzc = nxyzd

      nxyzr = nxyz1


      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ir = 1    ! pointer to scalar field \rho
      ib = 1    ! pointer to scalar field Bdu


      do e=1,nelv

!       Interpolate convecting field      
        call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0) ! 0 --> forward
        call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0) ! 0 --> forward
        if (if3d) call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Radius to fine mesh
        call intp_rstd(radf,ym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate density to fine mesh
        call intp_rstd(rhof,rho(ir),lx1,lxd,if3d,0)      ! 0 --> forward
       
        if (if3d) then  ! Convert convector F to r-s-t coordinates

          do i=1,nxyzd
            tr(i,1) = radf(i)*rx(i,1,e)*fx(i) + radf(i)*rx(i,2,e)*fy(i)
     $                   + rx(i,3,e)*fz(i)
            tr(i,2) = radf(i)*rx(i,4,e)*fx(i) + radf(i)*rx(i,5,e)*fy(i)
     $                   + rx(i,6,e)*fz(i)
            tr(i,3) = radf(i)*rx(i,7,e)*fx(i) + radf(i)*rx(i,8,e)*fy(i)
     $                   + rx(i,9,e)*fz(i)
          enddo

!         Collocate with density on fine mesh
          call col2(tr(i,1),rhof,nxyzd)
          call col2(tr(i,2),rhof,nxyzd)
          call col2(tr(i,3),rhof,nxyzd)

        else

          do i=1,nxyzd
            tr(i,1) = radf(i)*rx(i,1,e)*fx(i) + radf(i)*rx(i,2,e)*fy(i)
            tr(i,2) = radf(i)*rx(i,3,e)*fx(i) + radf(i)*rx(i,4,e)*fy(i)
          enddo

!         Collocate with density on fine mesh
          call col2(tr(i,1),rhof,nxyzd)
          call col2(tr(i,2),rhof,nxyzd)

        endif

!       Interpolate convected field      
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward
!       Gradients on the Fine Reference mesh.
        call grad_rst(ur,us,ut,uf,lxd,if3d)

        if (if3d) then
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
             uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
          enddo
        else
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
             uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
          enddo
        endif
        call intp_rstd(bdu(ib),uf,lx1,lxd,if3d,1) ! Project back to coarse

        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1
        ir = ir + nxyzr

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dealias_rho_uv(ruv,rho,u,v)

!     Compute dealiased form:  J^T Bf *Jv .Ju w/ correct Jacobians

!     For the cylindrical case, the radius in the denominator gets
!     cancelled by the radius coming from the jacobian.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      real ruv(1),rho(1),u(1),v(1)

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
      integer iu,iv,iuv,ir

      integer nxyz1,nxyzv,nxyzd,nxyzu,nxyzr

      real zd,wd
      common /dealias1/ zd(lxd),wd(lxd)

      integer i,j,l

      real rhof(ltd)

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzv = nxyz1
!      if (ifcf) nxyzc = nxyzd

      nxyzr = nxyz1

      iu  = 1    ! pointer to scalar field u
      iv  = 1    ! pointer to vector field v
      ir  = 1    ! pointer to vector field \rho
      iuv = 1    ! pointer to scalar uv 

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

!       Interpolate v 
        call intp_rstd(fx,v(iv),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate u 
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate \rho 
        call intp_rstd(rhof,rho(ir),lx1,lxd,if3d,0) ! 0 --> forward
       
!       Interpolate Jacobian (Probably only needs to be done once) 
        call intp_rstd(jacm1d,jacm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> forward

        do i=1,nxyzd
          tr(i) = rhof(i)*wd2(i)*jacm1d(i)*uf(i)*fx(i)
        enddo

        call intp_rstd(ruv(iuv),tr,lx1,lxd,if3d,1) ! Project back to coarse

        iv  = iv  + nxyzv
        iu  = iu  + nxyzu
        ir  = ir  + nxyzr
        iuv = iuv + nxyz1

      enddo

      return
      end subroutine dealias_rho_uv
!-----------------------------------------------------------------------

      subroutine convect_cylindrical(bdu,u,cx,cy,cz)

!     Compute dealiased form:  J^T Bf *JC .Grad Ju w/ correct Jacobians
!     Grad is the cylindrical grad operator

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT' 

      real bdu(1),u(1),cx(1),cy(1),cz(1)

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e,i
      integer iu,ic,ib
      integer nxyzu,nxyzc,nxyz1,nxyzd

      real radf(ltd)

!     Geometric factors Mapping.
!     3D:      
!     rxm1(1,1,1,e) -->  rx(1,1,e)
!     rym1(1,1,1,e) -->  rx(1,2,e)
!     rzm1(1,1,1,e) -->  rx(1,3,e)
!     sxm1(1,1,1,e) -->  rx(1,4,e)
!     sym1(1,1,1,e) -->  rx(1,5,e)
!     szm1(1,1,1,e) -->  rx(1,6,e)
!     txm1(1,1,1,e) -->  rx(1,7,e)
!     tym1(1,1,1,e) -->  rx(1,8,e)
!     tzm1(1,1,1,e) -->  rx(1,9,e)

!     2D:
!     rxm1(1,1,1,e) -->  rx(1,1,e)
!     rym1(1,1,1,e) -->  rx(1,2,e)
!     sxm1(1,1,1,e) -->  rx(1,3,e)
!     sym1(1,1,1,e) -->  rx(1,4,e)

   
   
!     Put the geometrix factors on the fine mesh
!     Also includes multiplication by Weights      
      call set_dealias_rx

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
!      if (ifcf) nxyzc = nxyzd

      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu


      do e=1,nelv

!       Interpolate convecting field      
        call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0) ! 0 --> forward
        call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0) ! 0 --> forward
        if (if3d) call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Radius to fine mesh
        call intp_rstd(radf,ym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> forward

        if (if3d) then  ! Convert convector F to r-s-t coordinates

          do i=1,nxyzd
            tr(i,1) = radf(i)*rx(i,1,e)*fx(i) + radf(i)*rx(i,2,e)*fy(i)
     $                   + rx(i,3,e)*fz(i)
            tr(i,2) = radf(i)*rx(i,4,e)*fx(i) + radf(i)*rx(i,5,e)*fy(i)
     $                   + rx(i,6,e)*fz(i)
            tr(i,3) = radf(i)*rx(i,7,e)*fx(i) + radf(i)*rx(i,8,e)*fy(i)
     $                   + rx(i,9,e)*fz(i)
          enddo

        else

          do i=1,nxyzd
            tr(i,1) = radf(i)*rx(i,1,e)*fx(i) + radf(i)*rx(i,2,e)*fy(i)
            tr(i,2) = radf(i)*rx(i,3,e)*fx(i) + radf(i)*rx(i,4,e)*fy(i)
          enddo

        endif

!       Interpolate convected field      
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward
!       Gradients on the Fine Reference mesh.
        call grad_rst(ur,us,ut,uf,lxd,if3d)

        if (if3d) then
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
             uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
          enddo
        else
          do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
             uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
          enddo
        endif
        call intp_rstd(bdu(ib),uf,lx1,lxd,if3d,1) ! Project back to coarse

        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dealias_uv(uv,u,v)

!     Compute dealiased form:  J^T Bf *Jv .Ju w/ correct Jacobians

!     For the cylindrical case, the radius in the denominator gets
!     cancelled by the radius coming from the jacobian.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      real uv(1),u(1),v(1)

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
      integer iu,iv,iuv

      integer nxyz1,nxyzv,nxyzd,nxyzu

      real zd,wd
      common /dealias1/ zd(lxd),wd(lxd)

      integer i,j,l

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzv = nxyz1
!      if (ifcf) nxyzc = nxyzd

      iu  = 1    ! pointer to scalar field u
      iv  = 1    ! pointer to vector field v
      iuv = 1    ! pointer to scalar uv 

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

!       Interpolate v 
        call intp_rstd(fx,v(iv),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate u 
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Jacobian (Probably only needs to be done once) 
        call intp_rstd(jacm1d,jacm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> forward

        do i=1,nxyzd
          tr(i) = wd2(i)*jacm1d(i)*uf(i)*fx(i)
        enddo
        call intp_rstd(uv(iuv),tr,lx1,lxd,if3d,1) ! Project back to coarse

        iv  = iv  + nxyzv
        iu  = iu  + nxyzu
        iuv = iuv + nxyz1

      enddo

      return
      end subroutine dealias_uv
!-----------------------------------------------------------------------

      subroutine axhmsf_cyl(Au1,Au2,Au3,u1,u2,u3,h1,h2)

!     Fluid (MATMOD .GE. 0) :  Hij Uj = Aij*Uj + H2*B*Ui 

      implicit none

      include 'SIZE'
      include 'TSTEP'   ! nelfld
      include 'MASS'

      real u1(1),u2(1),u3(1)
      real Au1(1),Au2(1),Au3(1)
      real h1(1),h2(1)

      integer matmod,nel,ntot1

      logical ifdfrm,iffast,ifh2,ifsolv
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv


      matmod = 0                ! Newtonian Fluids
      nel    = nelfld(ifield)   ! ifield should be already set

      ntot1 = lx1*ly1*lz1*nel


!     Common blocks are used in succession.

!     Aij*Uj      
      call stnrate_cyl(u1,u2,u3,nel,matmod)
      call stress_cyl (h1,h2,nel,matmod)
      call div_stress_cyl(Au1,Au2,Au3,nel)   ! aijuj

!     Add Helmholtz contributions
!     + H2*B*U      
      call addcol4 (Au1,bm1,h2,u1,ntot1)
      call addcol4 (Au2,bm1,h2,u2,ntot1)
      call addcol4 (Au3,bm1,h2,u3,ntot1)

      return
      end subroutine axhmsf_cyl             
!-----------------------------------------------------------------------
      subroutine stnrate_cyl (u1,u2,u3,nel,matmod)

C     Compute strainrates
C     CAUTION : Stresses and strainrates share the same scratch commons
      
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'

      real exx,exy,exz,eyx
      common /ctmp1/ exx(lx1*ly1*lz1*lelt)
     $             , exy(lx1*ly1*lz1*lelt)
     $             , exz(lx1*ly1*lz1*lelt)
     $             , eyx(lx1*ly1*lz1*lelt)
      real eyy,eyz
      common /ctmp0/ eyy(lx1*ly1*lz1*lelt)
     $             , eyz(lx1*ly1*lz1*lelt)

      real ezx,ezy,ezz
      common /axcyl_ez/ ezx(lx1*ly1*lz1*lelt)
     $                , ezy(lx1*ly1*lz1*lelt)
     $                , ezz(lx1*ly1*lz1*lelt)

      real radi,wk
      common /axcyl_wk/ radi(lx1*ly1*lz1*lelt)  ! radius inverse
     $                , wk(lx1*ly1*lz1*lelt)    ! work array

      real u1(lx1,ly1,lz1,1)
      real u2(lx1,ly1,lz1,1)
      real u3(lx1,ly1,lz1,1)

      integer nt1,nel,matmod

      nt1 = lx1*ly1*lz1*nel

!     1/R      
      call invers2(radi,ym1,nt1)

      call rzero3 (exx,exy,exz,nt1)
      call rzero3 (eyx,eyy,eyz,nt1)
      call rzero3 (ezx,ezy,ezz,nt1)

!     uxyz does not zero out variables.
!     values are just added on

!     X-Direction      
      call uxyz  (u1,exx,exy,exz,nel)
!     exz = 1/R*du/\theta
      call col2(exz,radi,nt1)      

      call invcol2 (exx,jacm1,nt1)
      call invcol2 (exy,jacm1,nt1)
      call invcol2 (exz,jacm1,nt1)

!     R-Direction      
      call uxyz  (u2,eyx,eyy,eyz,nel)
!     eyz = 1/R*dv/\theta
      call col2(eyz,radi,nt1)

      call invcol2 (eyx,jacm1,nt1)
      call invcol2 (eyy,jacm1,nt1)
      call invcol2 (eyz,jacm1,nt1)


!     \theta-Direction      
      if (ldim.eq.3) then
        call uxyz   (u3,ezx,ezy,ezz,nel)

        call invcol2 (ezx,jacm1,nt1)
        call invcol2 (ezy,jacm1,nt1)
        call invcol2 (ezz,jacm1,nt1)

!       ezz = 1/R*dw/\theta
        call col2(ezz,radi,nt1)

!       ezz = 1/R*dw/\theta + v/R
        call col3(wk,u2,radi,nt1)
        call add2(ezz,wk,nt1)

!       ezy = dw/dR - w/R
        call col3(wk,u3,radi,nt1)
        call sub2(ezy,wk,nt1)
      endif   

!     1/2*(\grad u + (\grad u)^T)
!     1/2*(exx + exx) = exx

!     1/2*(exy + eyx)
      call copy(wk,exy,nt1)
      call add2(exy,eyx,nt1)
      call add2(eyx,wk,nt1)
      call cmult(exy,0.5,nt1)
      call cmult(eyx,0.5,nt1)

!     1/2*(exz + exz)
      call copy(wk,exz,nt1)
      call add2(exz,ezx,nt1)
      call add2(ezx,wk,nt1)
      call cmult(exz,0.5,nt1)
      call cmult(ezx,0.5,nt1)

!     1/2*(eyx + exy)
!     Already calculated
      
!     1/2*(eyy + eyy) = eyy

!     1/2*(eyz + ezy)      
      call copy(wk,eyz,nt1)
      call add2(eyz,ezy,nt1)
      call add2(ezy,wk,nt1)
      call cmult(eyz,0.5,nt1)
      call cmult(ezy,0.5,nt1)

!     1/2*(ezx + exz)
!     Already calculated

!     1/2*(ezy + eyz)
!     Already calculated

!     1/2*(ezz + ezz) = ezz      


      return
      end subroutine stnrate_cyl
!-----------------------------------------------------------------------

      subroutine stress_cyl (h1,h2,nel,matmod)

C     MATMOD.GE.0        Fluid material models
C     MATMOD.LT.0        Solid material models
C
C     CAUTION : Stresses and strainrates share the same scratch commons

!     Tij = 2*\mu*Eij

      implicit none

      include 'SIZE'

      real exx,exy,exz,eyx
      common /ctmp1/ exx(lx1*ly1*lz1*lelt)
     $             , exy(lx1*ly1*lz1*lelt)
     $             , exz(lx1*ly1*lz1*lelt)
     $             , eyx(lx1*ly1*lz1*lelt)
      real eyy,eyz
      common /ctmp0/ eyy(lx1*ly1*lz1*lelt)
     $             , eyz(lx1*ly1*lz1*lelt)

      real ezx,ezy,ezz
      common /Axcyl_ez/ ezx(lx1*ly1*lz1*lelt)
     $                , ezy(lx1*ly1*lz1*lelt)
     $                , ezz(lx1*ly1*lz1*lelt)

      real radi,wk
      common /Axcyl_wk/ radi(lx1*ly1*lz1*lelt)  ! radius inverse
     $                , wk(lx1*ly1*lz1*lelt)    ! work array

      real t11,t22,t33,hii
      common /scrsf/ t11(lx1,ly1,lz1,lelt)
     $             , t22(lx1,ly1,lz1,lelt)
     $             , t33(lx1,ly1,lz1,lelt)
     $             , hii(lx1,ly1,lz1,lelt)

      real h1(lx1,ly1,lz1,1),h2(lx1,ly1,lz1,1)

      integer nt1,nel,matmod
      real const

      nt1 = lx1*ly1*lz1*nel

      if (matmod.eq.0) then

!        newtonian fluids
         const = 2.0
         call cmult2 (hii,h1,const,nt1)

         call col2   (exx,hii,nt1)
         call col2   (exy,hii,nt1)
         call col2   (exz,hii,nt1)

         call col2   (eyx,hii,nt1)
         call col2   (eyy,hii,nt1)
         call col2   (eyz,hii,nt1)

         call col2   (ezx,hii,nt1)
         call col2   (ezy,hii,nt1)
         call col2   (ezz,hii,nt1)
        

      elseif (matmod.eq.-1) then
!       elastic solids
        if (nio.eq.0) then
          write(6,*) 'Cylindrical solver not implemented for solids'
        endif
        call exitt
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine div_stress_cyl (au1,au2,au3,nel)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      real exx,exy,exz,eyx
      common /ctmp1/ exx(lx1*ly1*lz1*lelt)
     $             , exy(lx1*ly1*lz1*lelt)
     $             , exz(lx1*ly1*lz1*lelt)
     $             , eyx(lx1*ly1*lz1*lelt)
      real eyy,eyz
      common /ctmp0/ eyy(lx1*ly1*lz1*lelt)
     $             , eyz(lx1*ly1*lz1*lelt)

      real ezx,ezy,ezz
      common /Axcyl_ez/ ezx(lx1*ly1*lz1*lelt)
     $                , ezy(lx1*ly1*lz1*lelt)
     $                , ezz(lx1*ly1*lz1*lelt)

      real radi,wk
      common /Axcyl_wk/ radi(lx1*ly1*lz1*lelt)  ! radius inverse
     $                , wk(lx1*ly1*lz1*lelt)    ! work array

      real au1(lx1,ly1,lz1,1),
     $     au2(lx1,ly1,lz1,1),
     $     au3(lx1,ly1,lz1,1)

      integer nel
      integer nt1

      nt1 = lx1*ly1*lz1*nel

!     X-Direction      
      call ttxyz_cyl(au1,exx,exy,exz,nel)
!     No additional terms in x

!     R-Direction      
      call ttxyz_cyl(au2,eyx,eyy,eyz,nel)

!     ezz = 1/R*dw/\theta + v/R     
!     (v_hat/R)*(ezz)*W*Jac*R
!     Assuming BM1 contains the factor of R
!     (and also Jac. of course)      
      call col3(wk,ezz,radi,nt1)    ! wk = ezz/R
      call col2(wk,bm1,nt1)         ! wk = (ezz/R)*W*Jac*R

!     + (v_hat/R)*(ezz)*W*Jac*R
      call add2(au2,wk,nt1)

!     \theta-Direction      
      call ttxyz_cyl(au3,ezx,ezy,ezz,nel)

!     eyz = (1/R*dv/d\theta + dw/dR - w/R)
!     (v_hat/R)*(eyz)*W*Jac*R
!     Assuming BM1 contains the factor of R
!     (and also Jac. of course)      
      call col3(wk,eyz,radi,nt1)    ! wk = eyz/R
      call col2(wk,bm1,nt1)         ! wk = (eyz/R)*W*Jac*R

!     - (v_hat/R)*(eyz)*W*Jac*R
      call add2(au3,wk,nt1)

      return
      end
!-----------------------------------------------------------------------

      subroutine ttxyz_cyl (ff,tx,ty,tz,nel)

      implicit none

      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'WZ'

      real tx(lx1,ly1,lz1,1),
     $     ty(lx1,ly1,lz1,1),
     $     tz(lx1,ly1,lz1,1),
     $     ff(lx1*ly1*lz1,1)

      real fr,fs,ft
      common /scrsf/ fr(lx1*ly1*lz1,lelt)
     $             , fs(lx1*ly1*lz1,lelt)
     $             , ft(lx1*ly1*lz1,lelt)

      real ys(lx1)

      integer nel,isd
      integer iel,ix,iy,iz
      integer nxyz1,nt1

      real radi,wk
      common /Axcyl_wk/ radi(lx1*ly1*lz1*lelt)  ! radius inverse
     $                , wk(lx1*ly1*lz1*lelt)    ! work array

      nxyz1 = lx1*ly1*lz1
      nt1   = nxyz1*nel

!     1/R      
      call invers2(radi,ym1,nt1)

      call col3    (fr,rxm1,tx,nt1)
      call addcol3 (fr,rym1,ty,nt1)
      call col3    (fs,sxm1,tx,nt1)
      call addcol3 (fs,sym1,ty,nt1)

      if (ldim.eq.3) then
!       wk = 1/R*dr/dz
        call col3(wk,rzm1,tz,nt1)
        call col2(wk,radi,nt1)
        call addcol3 (fr,wk,tz,nt1)

!       wk = 1/R*ds/dz
        call col3(wk,szm1,tz,nt1)
        call col2(wk,radi,nt1)
        call addcol3 (fs,wk,tz,nt1)

        call col3    (ft,txm1,tx,nt1)
        call addcol3 (ft,tym1,ty,nt1)
!       wk = 1/R*dt/dz
        call col3(wk,tzm1,tz,nt1)
        call col2(wk,radi,nt1)
        call addcol3 (ft,wk,tz,nt1)
      endif

!     Multiply by Radius            
      call col2(fr(1,iel),ym1,nt1)
      call col2(fs(1,iel),ym1,nt1)
      if (ldim.eq.3) call col2(ft(1,iel),ym1,nt1)

!     Multiply by Weights          
      do iel=1,nel
        call col2(fr(1,iel),w3m1,nxyz1)
        call col2(fs(1,iel),w3m1,nxyz1)
        call col2(ft(1,iel),w3m1,nxyz1)
      enddo

!     grad(v).(tx,ty,tz)       
      do iel=1,nel
        call ttrst (ff(1,iel),fr(1,iel),fs(1,iel),
     $              ft(1,iel),wk)
      enddo

!     Additional terms due to gradients of unit vectors are done
!     outside this routine since they require cross terms 


      return
      end
c-----------------------------------------------------------------------





















