!======================================================================
!     Routines for introducing third component in a 2d simulation
!     Author: Prabal S. Negi
!     Right now we assume the third component is homogeneous.
!     Later, a fourier dependence can be added for the linearized solve.
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine initp_f3d()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include 'F3D'

      integer icalld
      save icalld
      data icalld /0/

      integer nxyz,ntot1
      integer p41,p42,p43,p44
      logical ifaxis_init           ! Did we initialize with ifaxis
      

      if (.not.iff3d) return

      if (nio.eq.0) then
        write(6,'(A18)') 'F3D: Initializing'
      endif

      ifaxis_init = ifaxis

      if (mod(npert,2).eq.1) then
        if (nio.eq.0) then
          write(6,'(A7,1x,I2)') 'NPERT =', npert
          write(6,*) 'Need both real and imaginary parts'
          write(6,*) 'Ensure NPERT = 2 for 3D perturbation solve'
        endif  
        call exitt
      endif

      if (ifcyl_f3d.and..not.ifaxis) then
        if (nio.eq.0) then
          write(6,*) 'IFAXIS:', IFAXIS
          write(6,*) 'IFCYL_F3D:', IFCYL_F3D
          write(6,*) 
     $      'IFAXIS initialization required for Cylindrical solver'
        endif
        call exitt
      endif  

      if (.not.ifcyl_f3d.and.ifaxis) then
        if (nio.eq.0) then
          write(6,*) 'IFAXIS:', IFAXIS
          write(6,*) 'IFCYL_F3D:', IFCYL_F3D
          write(6,*) 
     $      'IFAXIS initialized without Cylindrical solver'
        endif
        call exitt
      endif  


!     Only Homogeneous nonlinear solver implemented.
      if (npert.eq.0) k_f3d = 0.0

      nxyz  = lx1*ly1*lz1
      ntot1 = nxyz*nelv
     
      call init_pertfld_f3d() 

!     Need to initialize some variables
!     V3MASK
      call copy(v3mask,v1mask,ntot1)

!     Velocities can be initialized from useric. 

      if (nio.eq.0) then
        write(6,'(A19,1x,G10.3)') 'IFF3D: Wavenumber=',k_f3d
        write(6,'(A7,1x,L1)') 'IFF3D:', iff3d
        write(6,'(A11,1x,L1)') 'IFCYL_F3D:',ifcyl_f3d

!       Set ifaxis to false 
        write(6,'(A23,1x,L1)') 'IFAXIS initialization:',ifaxis_init
        if (ifcyl_f3d) ifaxis = .false.
        write(6,'(A8,1x,L1)') 'IFAXIS:',ifaxis

        write(6,'(A24)') 'F3D: Initializion done.'

!       Write out Preconditioner settings:
        p41 = param(41)
        p42 = param(42)
        p43 = param(43)
        p44 = param(44)

        write(6,'(A24)') 'Preconditioner Settings'
        if (p41.eq.0) then
          write(6,'(A11,1x,I1,A44)') 'Param(41)=',p41,
     $      ': Additive Spectral-Element Multigrid (SEMG)'
        elseif(p41.eq.1) then
          write(6,'(A11,1x,I1,A42)') 'Param(41)=',p41,
     $      ': Hybrid Spectral-Element Multigrid (SEMG)'
        else
          write(6,'(A11,1x,I1)') 'Param(41)=',p41
        endif  

        if (p42.eq.0) then
          write(6,'(A11,1x,I1,A33)') 'Param(42)=',p42,
     $      ': GMRES with nonsymmetric weights'
        elseif(p42.eq.1) then
          write(6,'(A11,1x,I1,A21)') 'Param(42)=',p42,
     $      ': PCG without weights'
        else
          write(6,'(A11,1x,I1)') 'Param(42)=',p42
        endif  

        if (p43.eq.0) then
          write(6,'(A11,1x,I1,A28)') 'Param(43)=',p43,
     $      ': Additive multilevel scheme'
        elseif(p43.eq.1) then
          write(6,'(A11,1x,I1,A27)') 'Param(43)=',p43,
     $      ': Original two level scheme'
        else
          write(6,'(A11,1x,I1)') 'Param(43)=',p43
        endif  

        if (p44.eq.0) then
          write(6,'(A11,1x,I1,A47)') 'Param(44)=',p44,
     $      ': Top level Schwartz based on restrictions of E'
        elseif(p44.eq.1) then
          write(6,'(A11,1x,I1,A47)') 'Param(44)=',p44,
     $      ': Top level Schwartz based on restrictions of A'
        else
          write(6,'(A11,1x,I1)') 'Param(44)=',p44
        endif  

      endif

      return
      end subroutine initp_f3d
!----------------------------------------------------------------------

      subroutine fluidp_f3d (igeom)

!     Driver for perturbation velocity

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      include 'F3D'

      include 'TEST'          ! test arrays

      integer igeom
      integer jp0,jp2
      integer i,j

!      if (ifcyl_f3d) then
        call fluidp_cyl(igeom)
        return
!      endif  

      jp = 0
      do i = 1,npert,2
        jp0 = i

!       Solve Momentum            
        do j = 0,1
          jp = jp0 + j    
          if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1      format(i9,1pe14.7,' Perturbation Solve (Momentum):',i5)

          call perturbv_mom_f3d(igeom)

        enddo

!       Solve Pressure            
        jp = jp0 + j 
        if (nio.eq.0.and.igeom.eq.2) write(6,2) istep,time,jp
   2    format(i9,1pe14.7,' Perturbation Solve (Pressure):',i5)

        call incomprp_f3d(igeom) 
      
!       Pressure/Velocity Correction            
        do j = 1,npert,2
          jp = jp0 + j    
          call velpr_update_f3d(igeom)
        enddo
      enddo ! i=1,jp2


      jp=0   ! set jp to zero, for baseline flow

      return
      end subroutine fluidp_f3d
!-----------------------------------------------------------------------
      subroutine fluidp_cyl (igeom)

!     Driver for perturbation velocity
!     In Cylindrical Coordinates

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'

      include 'F3D'

      include 'TEST'

      integer igeom
      integer jp0,jp2
      integer i,j

      integer ntot

      jp2 = npert/2
      jp = 0
      do i = 0,jp2-1

        jp0 = 2*i

!       Build RHS for momentum equations 
        do j = 1,2
          jp = jp0 + j    
          if (nio.eq.0.and.igeom.eq.2) write(6,1) istep,time,jp
   1      format(i9,1pe14.7,' Build RHS (Momentum):',i5)

          call buildrhs_cyl(igeom)
        enddo

!       Solve momentum equations        
        if (igeom.gt.1) then
          jp = jp0+1
          call solvemom_cyl(igeom)
        endif

!       Solve Pressure
        if (igeom.gt.1) then
          
          jp = jp0+1 
          call incomprp_cyl(igeom) 

!         Pressure/Velocity Correction            
          jp = jp0+1
          call velpr_update_f3d(igeom)

        endif     ! ifgeom

      enddo ! i=1,npert,2


      jp=0   ! set jp to zero, for baseline flow

      return
      end subroutine fluidp_cyl
!-----------------------------------------------------------------------
      subroutine perturbv_mom_f3d (igeom)

      implicit none

!     Solve the convection-diffusion equation for the perturbation field, 
!     with projection onto a div-free space.


      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      real resv1,resv2,resv3
      real dv1,dv2,dv3
      real h1,h2

      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      integer intype,igeom
      integer ntot1


      ifield = 1

      ntot1 = lx1*ly1*lz1*nelv

      if (igeom.eq.1) then

!        Old geometry, old velocity

         call makefp_f3d
         call lagfieldp_f3d

!        Add third component of convective term   
!         call advab_w_f3d   

      else
c
c        New geometry, new velocity
c
         intype = -1
         call sethlm_f3dp(h1,h2,intype)
         call cresvipp_f3d(resv1,resv2,resv3,h1,h2)

         if3d = .true.   
         call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         if3d = .false.

         call add2(vxp(1,jp),dv1,ntot1)
         call add2(vyp(1,jp),dv2,ntot1)
         call add2(vzp(1,jp),dv3,ntot1)

      endif

      return
      end subroutine perturbv_mom_f3d
!-----------------------------------------------------------------------

      subroutine buildrhs_cyl (igeom)

      implicit none

!     Solve the convection-diffusion equation for the perturbation field, 
!     with projection onto a div-free space.


      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      integer intype,igeom
      integer ntot1


      ifield = 1

      if (igeom.eq.1) then

!        Old geometry, old velocity

         call makefp_f3d
         call lagfieldp_f3d

      endif

      return
      end subroutine buildrhs_cyl

!---------------------------------------------------------------------- 
      subroutine solvemom_cyl(igeom)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      include 'TEST'      

      real dv1,dv2,dv3
      real h1,h2

      real resv1r,resv2r,resv3r
      real dv1r,dv2r,dv3r
      common /scrns1/  resv1r(lx1,ly1,lz1,lelv)
     $ ,               resv2r(lx1,ly1,lz1,lelv)
     $ ,               resv3r(lx1,ly1,lz1,lelv)
     $ ,               dv1r  (lx1,ly1,lz1,lelv)
     $ ,               dv2r  (lx1,ly1,lz1,lelv)
     $ ,               dv3r  (lx1,ly1,lz1,lelv)

      real resv1i,resv2i,resv3i
      real dv1i,dv2i,dv3i
      common /scrns2/ resv1i(lx1,ly1,lz1,lelv)
     $ ,              resv2i(lx1,ly1,lz1,lelv)
     $ ,              resv3i(lx1,ly1,lz1,lelv)
     $ ,              dv1i  (lx1,ly1,lz1,lelv)
     $ ,              dv2i  (lx1,ly1,lz1,lelv)
     $ ,              dv3i  (lx1,ly1,lz1,lelv)
    
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      integer intype,igeom
      integer ntot1,ntot2

      integer jpr,jpi


      jpr = jp
      jpi = jp + 1

      if (igeom.gt.1) then
        if (mod(jp,2).ne.1) then
          if (nid.eq.0) write(6,*) 'something went wrong for jp loop'
          if (nid.eq.0) write(6,*) 'jp = ', jp
          if (nid.eq.0) write(6,*) 'Exiting in solvemom_cyl'
          call exitt
        endif          

        if (nio.eq.0.and.igeom.eq.2) write(6,2) istep,time,jpr,jpi
   2    format(i9,1pe14.7,
     $  ' Cylindrical Perturbation Solve (Momentum) jp:',i3,',',i3)

        ifield = 1

        ntot1 = lx1*ly1*lz1*nelv

        intype = -1
        call sethlm_f3dp(h1,h2,intype)

        call cresvipp_cyl(resv1r,resv2r,resv3r,
     $                    resv1i,resv2i,resv3i,h1,h2)

!       Solve        
        call ophinv_cyl(dv1r,dv2r,dv3r,dv1i,dv2i,dv3i,
     $                  resv1r,resv2r,resv3r,resv1i,resv2i,resv3i,
     $                  h1,h2,tolhv,nmxv)

!       Update Velocity (real)
        call add2_3(vxp(1,jpr),vyp(1,jpr),vzp(1,jpr),
     $              dv1r,dv2r,dv3r,ntot1)

!       Update Velocity (imaginary)
        call add2_3(vxp(1,jpi),vyp(1,jpi),vzp(1,jpi),
     $              dv1i,dv2i,dv3i,ntot1)


      endif

      return
      end subroutine solvemom_cyl
!-----------------------------------------------------------------------

      subroutine init_pertfld_f3d

      implicit none

      include 'SIZE'
      include 'SOLN'    ! jp
      include 'TSTEP'   ! ifield

      integer i,n2
      integer ifld

      n2 = lx2*ly2*lz2*lelv

      ifld = ifield
      ifield = 1
      do i = 1,npert
        jp = i
        call nekuic

!       zero out pressure            
        call rzero(prp(1,jp),n2)

        call dsavg(vxp(1,jp))
        call dsavg(vyp(1,jp))
        call dsavg(vzp(1,jp))
      enddo
      jp = 0
      ifield = ifld

      return
      end subroutine init_pertfld_f3d
!----------------------------------------------------------------------
      subroutine makefp_f3d

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      include 'F3D'

      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

!     Build user defined forcing
      if3d = .true.
      call makeufp_f3d
      if3d = .false.
      if (filterType.eq.2) call make_hpf

      if (ifcyl_f3d) then
        call advabp_cyl_f3d
      else
        call advabp_f3d
      endif
!      if (ifnav.and.(.not.ifchar).and.(ifadj)) call advabp_adjoint_f3dp
      if (iftran) call makextp_f3d
      call makebdfp_f3d



      return
      end subroutine makefp_f3d

!----------------------------------------------------------------------

      subroutine makeufp_f3d

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      include 'F3D'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,iel,i,j,k,ielg
      integer ijke


      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfxp(1,jp),ntot1)
      call rzero(bfyp(1,jp),ntot1)
      call rzero(bfzp(1,jp),ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            ijke = i+lx1*((j-1)+ly1*((k-1) + lz1*(iel-1)))
            bfxp(ijke,jp) = ffx
            bfyp(ijke,jp) = ffy
            bfzp(ijke,jp) = ffz
 100  continue

!     Not sure why we multiply by density
!      call col2(bfzp(1,jp),vtrans(1,1,1,1,ifield),nx1*ny1*nz1*nelv)

      call col2  (bfxp(1,jp),bm1,ntot1)
      call col2  (bfyp(1,jp),bm1,ntot1)
      call col2  (bfzp(1,jp),bm1,ntot1)
      time = time+dt

      return
      end subroutine makeufp_f3d

!-----------------------------------------------------------------------
      subroutine advabp_f3d

!     Eulerian scheme, add convection term to forcing function
!     at current time step.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include 'F3D'

      real ta1,ta2,ta3
      real tb1,tb2,tb3
      common /scrns/ ta1 (lx1*ly1*lz1*lelv)
     $ ,             ta2 (lx1*ly1*lz1*lelv)
     $ ,             ta3 (lx1*ly1*lz1*lelv)
     $ ,             tb1 (lx1*ly1*lz1*lelv)
     $ ,             tb2 (lx1*ly1*lz1*lelv)
     $ ,             tb3 (lx1*ly1*lz1*lelv)

      integer i,ntot1
      real tmp


      ntot1 = lx1*ly1*lz1*nelv

      if (if3d.or.iff3d) then
         call copy  (tb1,vx,ntot1)                   ! Save velocity
         call copy  (tb2,vy,ntot1)                   ! Save velocity
         call copy  (tb3,vz,ntot1)                   ! Save velocity

!        U <-- dU
         call copy  (vx,vxp(1,jp),ntot1)                   ! Save velocity
         call copy  (vy,vyp(1,jp),ntot1)                   ! Save velocity
         call copy  (vz,vzp(1,jp),ntot1)                   ! Save velocity

         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call convop  (ta3,tb3)

!        Restore velocity
         call copy  (vx,tb1,ntot1)
         call copy  (vy,tb2,ntot1)
         call copy  (vz,tb3,ntot1)

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))
         call convop  (ta3,vzp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
         enddo

!        Add z convection for all components
!        Assuming dU/dz, dV/dz, dW/dz = 0 of course
         if (mod(jp,2).eq.1) then
            call convect_w_f3d(ta1,vxp(1,jp+1),vz)           
            call convect_w_f3d(ta2,vyp(1,jp+1),vz)           
            call convect_w_f3d(ta3,vzp(1,jp+1),vz)           

            do i=1,ntot1
               tmp = -k_f3d*vtrans(i,1,1,1,ifield)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
               bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
            enddo

         else
            call convect_w_f3d(ta1,vxp(1,jp-1),vz)           
            call convect_w_f3d(ta2,vyp(1,jp-1),vz)           
            call convect_w_f3d(ta3,vzp(1,jp-1),vz)           

            do i=1,ntot1
               tmp = k_f3d*vtrans(i,1,1,1,ifield)
               bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
               bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
               bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
            enddo

         endif

      else  ! 2D without fourier third component

         call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
         call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
         call convop  (ta1,tb1)                                ! du.grad U
         call convop  (ta2,tb2)
         call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

         call convop  (ta1,vxp(1,jp))       !  U.grad dU
         call convop  (ta2,vyp(1,jp))

         do i=1,ntot1
            tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
         enddo

      endif


      return
      end subroutine advabp_f3d
!--------------------------------------------------------------------

      subroutine advabp_cyl_f3d

!     Eulerian scheme, add convection term to forcing function
!     at current time step.
!     Done for Cylindrical coordinates
!     Coordinate Transformation: y --> R
!                                x --> x
!                                z --> theta (Done with Fourier)        

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'GEOM'

      include 'F3D'

      real ta1,ta2,ta3
      real tb1,tb2,tb3
      common /scrns/ ta1 (lx1*ly1*lz1*lelv)
     $ ,             ta2 (lx1*ly1*lz1*lelv)
     $ ,             ta3 (lx1*ly1*lz1*lelv)
     $ ,             tb1 (lx1*ly1*lz1*lelv)
     $ ,             tb2 (lx1*ly1*lz1*lelv)
     $ ,             tb3 (lx1*ly1*lz1*lelv)

      real rinv,rhom
      common /scrch/ rinv(lx1*ly1*lz1*lelt),
     $               rhom(lx1*ly1*lz1*lelt)    

      integer i,ntot1
      real tmp,tmp2
      integer ji


      real tmp5(lx1,ly1,lz1,lelt)
      real tmp6(lx1,ly1,lz1,lelt)
      real tmp7(lx1,ly1,lz1,lelt)
      real tmp8(lx2,ly2,lz2,lelt)

      common /testtmp2/ tmp5,tmp6,tmp7,tmp8
      real glmax


      ji = mod(jp,2)

      ntot1 = lx1*ly1*lz1*nelv

      if (if3d.or.iff3d) then
        call copy  (tb1,vx,ntot1)                   ! Save velocity
        call copy  (tb2,vy,ntot1)                   ! Save velocity
        call copy  (tb3,vz,ntot1)                   ! Save velocity

!       U <-- dU
        call copy  (vx,vxp(1,jp),ntot1)                   ! Save velocity
        call copy  (vy,vyp(1,jp),ntot1)                   ! Save velocity
        call copy  (vz,vzp(1,jp),ntot1)                   ! Save velocity

        call convop  (ta1,tb1)                                ! du.grad U
        call convop  (ta2,tb2)
        call convop  (ta3,tb3)

!       Restore velocity
        call copy  (vx,tb1,ntot1)
        call copy  (vy,tb2,ntot1)
        call copy  (vz,tb3,ntot1)

        do i=1,ntot1
          tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
          bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
          bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
        enddo

        call convop  (ta1,vxp(1,jp))       !  U.grad dU
        call convop  (ta2,vyp(1,jp))
        call convop  (ta3,vzp(1,jp))

        do i=1,ntot1
          tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
          bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
          bfzp(i,jp) = bfzp(i,jp)-tmp*ta3(i)
        enddo

!       Add additional convection terms
!       Assuming dU/dz, dV/dz, dW/dz = 0

!       rinv = 1/R
        if (ifcyl_f3d) then 
          call invers2(rinv,ym1,ntot1)
        else
          call rone(rinv,ntot1) 
        endif  

!       rhom = \rho*BM1        
!        call col3(rhom,bm1,vtrans(1,1,1,1,ifield),ntot1)
        call copy(rhom,vtrans(1,1,1,1,ifield),ntot1)

        if (ji.eq.1) then
!         Real part

!         x component
          call convect_w_f3d(ta1,vxp(1,jp+1),vz)  ! U\theta*ux'_i

          do i=1,ntot1
            tmp = -k_f3d*rinv(i)*rhom(i)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          enddo

!         R component
          call convect_w_f3d(ta1,vyp(1,jp+1),vz)  ! U\theta*uR'_i
          call convect_w_f3d(ta2,vzp(1,jp),vz)    ! U\theta*u\theta'_r
          do i=1,ntot1
            tmp  = -k_f3d*rinv(i)*rhom(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta1(i)
            if (ifcyl_f3d) then
              tmp2       = -2.*rinv(i)*rhom(i)
              bfyp(i,jp) = bfyp(i,jp)-tmp2*ta2(i)
            endif
          enddo

!         \theta component
          call convect_w_f3d(ta1,vzp(1,jp+1),vz)  ! U\theta*u\theta'_i
          call convect_w_f3d(ta2,vyp(1,jp),vz)    ! U\theta*uR'_r
          call convect_w_f3d(ta3,vzp(1,jp),vy)    ! UR*u\theta'_r

          do i=1,ntot1
            tmp  = -k_f3d*rinv(i)*rhom(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta1(i)
            if (ifcyl_f3d) then
              tmp2   =  rinv(i)*rhom(i)
              bfzp(i,jp) = bfzp(i,jp)-tmp2*ta2(i)-tmp2*ta3(i)
            endif
          enddo
          
        else
!         Imaginary part

!         x component
          call convect_w_f3d(ta1,vxp(1,jp-1),vz)  ! U\theta*ux'_r

          do i=1,ntot1
            tmp = k_f3d*rinv(i)*rhom(i)
            bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          enddo

!         R component
          call convect_w_f3d(ta1,vyp(1,jp-1),vz)  ! U\theta*uR'_r
          call convect_w_f3d(ta2,vzp(1,jp),vz)    ! U\theta*u\theta'_i

          do i=1,ntot1
            tmp  =  k_f3d*rinv(i)*rhom(i)
            bfyp(i,jp) = bfyp(i,jp)-tmp*ta1(i)
            if (ifcyl_f3d) then
              tmp2   = -2.*rinv(i)*rhom(i)
              bfyp(i,jp) = bfyp(i,jp)-tmp2*ta2(i)
            endif      
          enddo

!         \theta component
          call convect_w_f3d(ta1,vzp(1,jp-1),vz)  ! U\theta*u\theta'_r
          call convect_w_f3d(ta2,vyp(1,jp),vz)    ! U\theta*uR'_i
          call convect_w_f3d(ta3,vzp(1,jp),vy)    ! UR*u\theta'_i

          do i=1,ntot1
            tmp  =  k_f3d*rinv(i)*rhom(i)
            bfzp(i,jp) = bfzp(i,jp)-tmp*ta1(i)
            if (ifcyl_f3d) then
              tmp2   =  rinv(i)*rhom(i)
              bfzp(i,jp) = bfzp(i,jp)-tmp2*ta2(i)-tmp2*ta3(i)
            endif
          enddo

        endif         ! ji.eq.1  

      else  ! 2D without Fourier third component

        call opcopy  (tb1,tb2,tb3,vx,vy,vz)                   ! Save velocity
        call opcopy  (vx,vy,vz,vxp(1,jp),vyp(1,jp),vzp(1,jp)) ! U <-- dU
        call convop  (ta1,tb1)                                ! du.grad U
        call convop  (ta2,tb2)
        call opcopy  (vx,vy,vz,tb1,tb2,tb3)  ! Restore velocity

        do i=1,ntot1
          tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
          bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
        enddo

        call convop  (ta1,vxp(1,jp))       !  U.grad dU
        call convop  (ta2,vyp(1,jp))

        do i=1,ntot1
          tmp = bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)
          bfxp(i,jp) = bfxp(i,jp)-tmp*ta1(i)
          bfyp(i,jp) = bfyp(i,jp)-tmp*ta2(i)
        enddo

      endif       ! if3d .or. iff3d

      return
      end subroutine advabp_cyl_f3d
!--------------------------------------------------------------------

      subroutine makextp_f3d

!     Add extrapolation terms to perturbation source terms

!     (nek5 equivalent for velocity is "makeabf")

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include 'F3D'

      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      ntot1 = lx1*ly1*lz1*nelv

      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,exx1p(1,jp),exx2p(1,jp),ab1,ab2,ntot1)
      call add3s2 (ta2,exy1p(1,jp),exy2p(1,jp),ab1,ab2,ntot1)
      call copy   (exx2p(1,jp),exx1p(1,jp),ntot1)
      call copy   (exy2p(1,jp),exy1p(1,jp),ntot1)
      call copy   (exx1p(1,jp),bfxp (1,jp),ntot1)
      call copy   (exy1p(1,jp),bfyp (1,jp),ntot1)
      call add2s1 (bfxp(1,jp),ta1,ab0,ntot1)
      call add2s1 (bfyp(1,jp),ta2,ab0,ntot1)
      if (if3d.or.iff3d) then
         call add3s2 (ta3,exz1p(1,jp),exz2p(1,jp),ab1,ab2,ntot1)
         call copy   (exz2p(1,jp),exz1p(1,jp),ntot1)
         call copy   (exz1p(1,jp),bfzp (1,jp),ntot1)
         call add2s1 (bfzp(1,jp),ta3,ab0,ntot1)
      endif
c
      return
      end subroutine makextp_f3d
!-----------------------------------------------------------------------

      subroutine makebdfp_f3d

!     Add contributions to perturbation source from lagged BD terms.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      include 'F3D'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)


      integer ilag,ntot1
      real const



      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
      call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      call opcolv3c (tb1,tb2,tb3
     $              ,vxp(1,jp),vyp(1,jp),vzp(1,jp),bm1,bd(2))

      if (iff3d) then
        call col3(tb3,vzp(1,jp),bm1,ntot1)
        call cmult(tb3,bd(2),ntot1)
      endif

!     Add contribution from lag terms
      do ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))
         else
            call opcolv3c(ta1,ta2,ta3,vxlagp(1,ilag-1,jp),
     $                                vylagp(1,ilag-1,jp),
     $                                vzlagp(1,ilag-1,jp),
     $                                bm1                   ,bd(ilag+1))
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
      enddo
      call opadd2col (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),tb1,tb2,tb3,h2)

!     Add contribution of lag terms to vzp
      if (iff3d) then
        do ilag=2,nbd
           if (ifgeom) then
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1lag(1,1,1,1,ilag-1),
     $                                ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           else
              call col3(ta3,vzlagp(1,ilag-1,jp),bm1,ntot1) 
              call cmult(ta3,bd(ilag+1),ntot1)           
           endif
           call add2(tb3,ta3,ntot1)
        enddo
        call add2col2(bfzp(1,jp),tb3,h2,ntot1)
      endif       ! iff3d


      return
      end subroutine makebdfp_f3d
!-----------------------------------------------------------------------

      subroutine lagfieldp_f3d

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      include 'F3D'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do ilag=nbdinp-1,2,-1
         call opcopy
     $     (vxlagp(1,ilag  ,jp),vylagp(1,ilag  ,jp),vzlagp(1,ilag  ,jp)
     $     ,vxlagp(1,ilag-1,jp),vylagp(1,ilag-1,jp),vzlagp(1,ilag-1,jp))
      enddo
      call opcopy(vxlagp(1,1,jp),vylagp(1,1,jp),vzlagp(1,1,jp)
     $           ,vxp   (1,jp)  ,vyp   (1,jp)  ,vzp   (1,jp) )

!     iff3d
      if (iff3d) then
        do ilag=nbdinp-1,2,-1
          call copy(vzlagp(1,ilag,jp),vzlagp(1,ilag-1,jp),ntot1)
        enddo
        call copy(vzlagp(1,1,jp),vzp(1,jp),ntot1)
      endif


      return
      end subroutine lagfieldp_f3d
!----------------------------------------------------------------------
      subroutine ophx_f3d (out1,out2,out3,inp1,inp2,inp3,h1,h2)

!     OUT = (H1*A+H2*B) * INP  

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include 'F3D'


      real out1 (lx1,ly1,lz1,1)
      real out2 (lx1,ly1,lz1,1)
      real out3 (lx1,ly1,lz1,1)
      real inp1 (lx1,ly1,lz1,1)
      real inp2 (lx1,ly1,lz1,1)
      real inp3 (lx1,ly1,lz1,1)
      real h1   (lx1,ly1,lz1,1)
      real h2   (lx1,ly1,lz1,1)

      integer imesh,matmod
      

      imesh = 1

      if (ifstrs) then
         matmod = 0
         call axhmsf (out1,out2,out3,inp1,inp2,inp3,h1,h2,matmod)
      else

!        the numbers are only needed for axis-symmetric formulation
!        need to come back to this later. 
         call axhelm (out1,inp1,h1,h2,imesh,1)
         call axhelm (out2,inp2,h1,h2,imesh,2)
         call axhelm (out3,inp3,h1,h2,imesh,3)
      endif

      return
      end subroutine ophx_f3d
!-----------------------------------------------------------------------
      subroutine cresvipp_f3d (resv1,resv2,resv3,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'    ! v?mask
      include 'MASS'    ! bm1

      include 'F3D'

      real           resv1 (lx1,ly1,lz1,lelv)
      real           resv2 (lx1,ly1,lz1,lelv)
      real           resv3 (lx1,ly1,lz1,lelv)
      real           h1    (lx1,ly1,lz1,lelv)
      real           h2    (lx1,ly1,lz1,lelv)

      real w1,w2,w3
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      real const
      

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if3d = .true.
      call bcdirvc (vxp(1,jp),vyp(1,jp),vzp(1,jp),
     $              v1mask,v2mask,v3mask)
      if3d = .false.

!     prabal. We don't care about traction conditions for now.
!     Maybe need to look at it if added stiffness terms are needed
!     Or if surface tension is needed
!      call bcneutr

      if (mod(jp,2).eq.1) then
!       Real part              

!       Need to Extrapolate both pressures at the same time    
        call extrapprp (prextr_f3d(1,1))
        call opgradt (resv1,resv2,resv3,prextr_f3d(1,1))


        jp = jp + 1
        call extrapprp (prextr_f3d(1,2))
        jp = jp - 1

        call map21_all_f3d(resv3,prextr_f3d(1,2))
        const = k_f3d
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      else
!       Imaginary part

        call extrapprp (prextr_f3d(1,2))
        call opgradt (resv1,resv2,resv3,prextr_f3d(1,2))

        jp = jp - 1
        call extrapprp (prextr_f3d(1,1))
        jp = jp + 1

        call map21_all_f3d(resv3,prextr_f3d(1,1))
        const = -k_f3d
        call cmult(resv3,const,ntot1)
        call col2(resv3,bm1,ntot1)
      endif    

      call add2(resv1,bfxp(1,jp),ntot1)
      call add2(resv2,bfyp(1,jp),ntot1)
      call add2(resv3,bfzp(1,jp),ntot1)

!     prabal
      call ophx_f3d(w1,w2,w3,vxp(1,jp),vyp(1,jp),
     $              vzp(1,jp),h1,h2)
      call sub2(resv1,w1,ntot1)
      call sub2(resv2,w2,ntot1)
      call sub2(resv3,w3,ntot1)

      return
      end subroutine cresvipp_f3d
!-----------------------------------------------------------------------

      subroutine cresvipp_cyl (resv1r,resv2r,resv3r,
     $                         resv1i,resv2i,resv3i,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      implicit none

      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'SOLN'    ! v?mask
      include 'MASS'    ! bm1

      include 'F3D'

      include 'TEST'

      real           resv1r(lx1,ly1,lz1,lelv)
      real           resv2r(lx1,ly1,lz1,lelv)
      real           resv3r(lx1,ly1,lz1,lelv)

      real           resv1i(lx1,ly1,lz1,lelv)
      real           resv2i(lx1,ly1,lz1,lelv)
      real           resv3i(lx1,ly1,lz1,lelv)

      real           h1(lx1,ly1,lz1,lelv)
      real           h2(lx1,ly1,lz1,lelv)

      real w1r,w2r,w3r
      common /scruz/ w1r(lx1,ly1,lz1,lelv)
     $ ,             w2r(lx1,ly1,lz1,lelv)
     $ ,             w3r(lx1,ly1,lz1,lelv)

      real w1i,w2i,w3i
      common /scruz2/ w1i(lx1,ly1,lz1,lelv)
     $ ,              w2i(lx1,ly1,lz1,lelv)
     $ ,              w3i(lx1,ly1,lz1,lelv)


      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      real const

      integer jpr,jpi           ! jp for real and imaginary parts


      jpr = jp
      jpi = jp + 1

      if (mod(jp,2).ne.1) then
        if (nid.eq.0) write(6,*) 'Something went wrong for jp loop'
        if (nid.eq.0) write(6,*) 'Exiting in cresvipp_cyl'
        call exitt
      endif          

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if3d = .true.
!     Real part      
      call bcdirvc_cyl(vxp(1,jpr),vyp(1,jpr),vzp(1,jpr),
     $              v1mask,v2mask,v3mask)

!     Imaginary part
      call bcdirvc_cyl(vxp(1,jpi),vyp(1,jpi),vzp(1,jpi),
     $              v1mask,v2mask,v3mask)
      if3d = .false.

!     prabal. We don't care about traction conditions for now.
!     Maybe need to look at it if added stiffness terms are needed
!     Or if surface tension is needed
!      call bcneutr

!     Real part              
!     Need to Extrapolate both pressures at the same time
      jp = jpr 
      call extrapprp (prextr_f3d(1,1))
      jp = jpi
      call extrapprp (prextr_f3d(1,2))
      jp = jpr

!     Note: theta gradient goes to imaginary part      
      call opgradt_f3d(resv1r,resv2r,resv3i,prextr_f3d(1,1))
!     Note: theta gradient goes to real part 
      call opgradt_f3d(resv1i,resv2i,resv3r,prextr_f3d(1,2))
      call chsign(resv3r,ntot1)     ! i*i = -1

!     Real      
      call add2_3(resv1r,resv2r,resv3r,
     $            bfxp(1,jpr),bfyp(1,jpr),bfzp(1,jpr),ntot1)
!     Imaginary      
      call add2_3(resv1i,resv2i,resv3i,
     $            bfxp(1,jpi),bfyp(1,jpi),bfzp(1,jpi),ntot1)


!     Ax
      call axhmsf_cyl(w1r,w2r,w3r,w1i,w2i,w3i,
     $                vxp(1,jpr),vyp(1,jpr),vzp(1,jpr),
     $                vxp(1,jpi),vyp(1,jpi),vzp(1,jpi),
     $                h1,h2)

      call sub2(resv1r,w1r,ntot1)
      call sub2(resv2r,w2r,ntot1)
      call sub2(resv3r,w3r,ntot1)

      call sub2(resv1i,w1i,ntot1)
      call sub2(resv2i,w2i,ntot1)
      call sub2(resv3i,w3i,ntot1)

      return
      end subroutine cresvipp_cyl
!-----------------------------------------------------------------------

      subroutine sethlm_f3dp(h1,h2,intype)

      implicit none

c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTYPE =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN' 
      include 'TSTEP'   ! ifield

      include 'F3D'

      real h1(1),h2(1)

      integer intype
      integer nel,ntot1

      real dtbd

      real k2

!     For the cylindrical case we use the stress formulation.
!     All relevant changes to the operators are handled during the
!     matrix vector products

      nel   = nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel

      k2    = k_f3d**2

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         if (intype.eq.0) then
            call rzero (h2,ntot1)
         else
            call cmult2(h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)

!            if (.not.ifcyl_f3d) then
!!             Add second derivative of the 3rd direction to the operator
!              call add2s2(h2,vdiff(1,1,1,1,ifield),k2,ntot1)
!            endif  
         endif
      else
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         call rzero (h2,ntot1)
      endif


      return
      end subroutine sethlm_f3dp

!----------------------------------------------------------------------
      subroutine incomprp_f3d (igeom)
c
c     Project U onto the closest incompressible field
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c

      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'CTIMER'

      include 'F3D'

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      integer ntot1,ntot2,intype,istart

      real dtbd,const

      integer jpi

      integer igeom

      if (igeom.eq.1) return

      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld

      jpi = 1
      if (mod(jp,2).eq.0) jpi = 2

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      dtbd   = bd(1)/dt

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (prcorr_f3d(1,jpi),vxp(1,jp),vyp(1,jp),vzp(1,jp))

      if (jpi.eq.1) then
        call map12_all_f3d(dp2,vzp(1,jpi+1))
        const = -k_f3d
      else
        call map12_all_f3d(dp2,vzp(1,jpi-1))
        const =  k_f3d
      endif

      call cmult(dp2,const,ntot2)

!      call add2s2(prcorr_f3d(1,jpi),dp2,const,ntot2)

      call chsign  (prcorr_f3d(1,jpi),ntot2)
      call ortho   (prcorr_f3d(1,jpi))


      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      intype =  2             ! Changing integration type here.
                              ! Need to modify cdabdtp accordingly
                              ! Also need to modify uzprec

      call esolver (prcorr_f3d(1,jpi),h1,h2,h2inv,intype)


      return
      end subroutine incomprp_f3d
!------------------------------------------------------------------------

      subroutine incomprp_cyl (igeom)
c
c     Project U onto the closest incompressible field
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c

      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'CTIMER'
      include 'GEOM'          ! YM1,YM2
      include 'MASS'

      include 'F3D'

      include 'TEST'

      real h1,h2,h2inv
      common /scrvh/ h1   (lx1,ly1,lz1,lelv)
     $ ,             h2   (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv(lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      real dummy
      common /scrcg/ dummy(lx1*ly1*lz1*lelt) 


      integer ntot1,ntot2,intype,istart

      real bddt,bddti,const

      integer jpr,jpi

      integer igeom

      real dnorm
      real divv,bdivv
      common /scruz/  divv (lx2,ly2,lz2,lelv)
     $ ,              bdivv(lx2,ly2,lz2,lelv)

      real wk1(lx1*ly1*lz1*lelv),wk2(lx1*ly1*lz1*lelv)
      real wk3(lx1*ly1*lz1*lelv)
      common /scrsf2/ wk1,wk2,wk3

      real glsc2

      if (igeom.eq.1) return

      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld

      jpr = jp
      jpi = jp+1

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      bddt   = bd(1)/dt
      bddti  = dt/bd(1)

      call rzero (h1,ntot1)
!      call copy  (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2(h2,vtrans(1,1,1,1,ifield),bddt,ntot1)
      call invers2(h2inv,h2,ntot1)

!     Note: OPDIV already contains the mass matrix multiplication

!     Note, we take imaginary part of vzp
      call cmult2(wk1,vzp(1,jpi),-1.0,ntot1)
      call opdiv_f3d(prcorr_f3d(1,1),vxp(1,jpr),vyp(1,jpr),wk1)
      call chsign(prcorr_f3d(1,1),ntot2)
      call ortho (prcorr_f3d(1,1))

!!     Imaginary part      
!     Note, we take real part of vzp
      call cmult2(wk1,vzp(1,jpr),1.0,ntot1)
      call opdiv_f3d(prcorr_f3d(1,2),vxp(1,jpi),vyp(1,jpi),wk1)
      call chsign(prcorr_f3d(1,2),ntot2)
      call ortho (prcorr_f3d(1,2))

      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      intype =  1             ! Changing integration type here.
                              ! Need to modify cdabdtp accordingly
                              ! Also need to modify uzprec

      if (nio.eq.0.and.igeom.eq.2) write(6,3) istep,time,jpr
      call esolver_f3d(prcorr_f3d(1,1),h1,h2,h2inv,intype)
 
      if (nio.eq.0.and.igeom.eq.2) write(6,3) istep,time,jpi
      call esolver_f3d(prcorr_f3d(1,2),h1,h2,h2inv,intype)

   3  format(i9,1pe14.7,' Perturbation Solve (Pressure):',i5)

      return
      end subroutine incomprp_cyl
!------------------------------------------------------------------------
      subroutine velpr_update_f3d (igeom)

!     Update Pressure and velocities based on pressure correction


      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'MASS'
      include 'CTIMER'
      include 'GEOM'          ! YM1

      include 'F3D'

      include 'TEST'

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns1/ w1    (lx1,ly1,lz1,lelv)
     $ ,              w2    (lx1,ly1,lz1,lelv)
     $ ,              w3    (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
     $ ,              dp    (lx2,ly2,lz2,lelv)

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      integer ntot1,ntot2

      real bddt,bddti,const

      integer jpr,jpi

      integer igeom

      if (igeom.eq.1) return

      jpr = jp
      jpi = jp + 1

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      bddt   = bd(1)/dt
      bddti  = dt/bd(1)

!     In principle these should still be preserved from the last routine.
!     But I calculate it again for clarity
      call rzero   (h1,ntot1)
!      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult2  (h2,vtrans(1,1,1,1,ifield),bddt,ntot1)
      call invers2 (h2inv,h2,ntot1)

!!    Update Pressure
      jp = jpr
      call lagpresp
      call add3(prp(1,jpr),prextr_f3d(1,1),prcorr_f3d(1,1),ntot2)

      jp = jpi
      call lagpresp
      call add3(prp(1,jpi),prextr_f3d(1,2),prcorr_f3d(1,2),ntot2)

      jp = jpr

      call opgradt_f3d(w1,w2,w3,prcorr_f3d(1,1))

      if3d = .true.
      call opbinv_f3d(dv1,dv2,dv3,w1,w2,w3,h2inv)
      if3d = .false.

      call add2(vxp(1,jpr),dv1,ntot1)
      call add2(vyp(1,jpr),dv2,ntot1)
!      call add2(vzp(1,jpr),dv3,ntot1)
      call add2(vzp(1,jpi),dv3,ntot1)           ! Imaginary part gets
                                                ! updated here

      call opgradt_f3d(w1,w2,w3,prcorr_f3d(1,2))

      if3d = .true.
      call opbinv_f3d(dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      if3d = .false.

      call add2(vxp(1,jpi),dv1,ntot1)
      call add2(vyp(1,jpi),dv2,ntot1)
!      call add2(vzp(1,jpi),dv3,ntot1)
      call sub2(vzp(1,jpr),dv3,ntot1)     ! Real part gets updated here
                                          ! Note the change from add to
                                          ! sub

                                          
      return
      end subroutine velpr_update_f3d
!------------------------------------------------------------------------
      subroutine map12_all_f3d(pm2,pm1)

      implicit none

      include 'SIZE'

      real pm1(lx1,ly1,lz1,lelv)
      real pm2(lx2,ly2,lz2,lelv)
      integer e

      do e=1,nelv
         call map12 (pm2(1,1,1,e),pm1(1,1,1,e),e)
      enddo
   
      return
      end subroutine map12_all_f3d
!-----------------------------------------------------------------------

      subroutine map21_all_f3d (y,x)

!     Map X from mesh M2 to mesh M1 (Y)
      
      implicit none

      INCLUDE 'SIZE'

      real x(lx2,ly2,lz2,lelv)
      real y(lx1,ly1,lz1,lelv)

      integer e

      do e=1,nelv
        call map21t(y(1,1,1,e),x(1,1,1,e),e)
      enddo

      return
      end subroutine map21_all_f3d

!-----------------------------------------------------------------------

      subroutine opbinv_f3d (out1,out2,out3,inp1,inp2,inp3,h2inv)

!     Compute OUT = (H2*B)-1 * INP   (explicit)


      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'

      include 'F3D'     ! iff3d

      real out1  (1)
      real out2  (1)
      real out3  (1)
      real inp1  (1)
      real inp2  (1)
      real inp3  (1)
      real h2inv (1)

      real tmp
      integer i,isbcnt,ntot

      include 'OPCTR'
C
#ifdef TIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'opbinv'
      endif
#endif


      ntot=lx1*ly1*lz1*nelv

!      call opmask  (inp1,inp2,inp3)
!      call opdssum (inp1,inp2,inp3)

!     mask      
      call col2_3(inp1,inp2,inp3,v1mask,v2mask,v3mask,ntot)
!     dssum
      call dssum3(inp1,inp2,inp3)


#ifdef TIMER
      isbcnt = ntot*(1+ldim)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

      call invcol3 (out1,bm1,h2inv,ntot)  ! this is expensive and should
      call dssum   (out1,lx1,ly1,lz1)     ! be changed (pff, 3/18/09)
      if (iff3d) then
         do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
            out3(i)=inp3(i)*tmp
         enddo
      else
         do i=1,ntot
            tmp = 1./out1(i)
            out1(i)=inp1(i)*tmp
            out2(i)=inp2(i)*tmp
         enddo
      endif

      return
      end subroutine opbinv_f3d
c-----------------------------------------------------------------------

      subroutine convect_w_f3d(cku,u,Cz)

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
      end subroutine convect_w_f3d
!-----------------------------------------------------------------------
      subroutine cdabdtp_f3d(ap,wp,h1,h2,h2inv,intype)

!     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
!     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p
!     INTYPE= 2  Compute the matrix-vector product    D(B/DT)(-1)DT*p
!                  with fourier in 3rd component 

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'MASS'
      include 'F3D'
      include 'GEOM'          ! YM2

      real           ap    (lx2,ly2,lz2,lelv)
      real           wp    (lx2,ly2,lz2,lelv)
      real           h1    (lx1,ly1,lz1,lelv)
      real           h2    (lx1,ly1,lz1,lelv)
      real           h2inv (lx1,ly1,lz1,lelv)

      real ta1,ta2,ta3,tb1,tb2,tb3
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
     $ ,             tb1 (lx1,ly1,lz1,lelv)
     $ ,             tb2 (lx1,ly1,lz1,lelv)
     $ ,             tb3 (lx1,ly1,lz1,lelv)

      real ttmp2(lx2,ly2,lz2,lelv)         ! lazy work. Should use a scratch array

      integer ntot1,ntot2,intype

      real const
      logical ifaxis_old

      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv

!!     (D^T)P
!     (pdv/dx; pdv/dR + pv/R; -kpv/R)*BM1
       call opgradt_f3d(ta1,ta2,ta3,wp)

!!    ((B*beta/dt)^-1)*(D^T)P
      if3d = .true.  ! Also do this for the third component
      call opbinv_f3d (tb1,tb2,tb3,ta1,ta2,ta3,h2inv)
      if3d = .false.

      call chsign(tb3,ntot1)  ! since we need i*i = -1

      call opdiv_f3d(ap,tb1,tb2,tb3)

      return
      end subroutine cdabdtp_f3d

!-----------------------------------------------------------------------

      subroutine bcdirvc_cyl(v1,v2,v3,mask1,mask2,mask3)

c     apply dirichlet boundary conditions to surface of vector (v1,v2,v3).
c     use ifield as a guide to which boundary conditions are to be applied.

      implicit none

      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
      INCLUDE 'CTIMER'

      include 'F3D'

      real tmp1,tmp2,tmp3
      common /scruz/ tmp1(lx1,ly1,lz1,lelv)
     $             , tmp2(lx1,ly1,lz1,lelv)
     $             , tmp3(lx1,ly1,lz1,lelv)

      real tmq1,tmq2,tmq3      
      common /scrmg/ tmq1(lx1,ly1,lz1,lelv)
     $             , tmq2(lx1,ly1,lz1,lelv)
     $             , tmq3(lx1,ly1,lz1,lelv)
c
      real v1(lx1,ly1,lz1,lelv),v2(lx1,ly1,lz1,lelv)
     $    ,v3(lx1,ly1,lz1,lelv)
      real mask1(lx1,ly1,lz1,lelv),mask2(lx1,ly1,lz1,lelv)
     $    ,mask3(lx1,ly1,lz1,lelv)
c
      common  /nekcb/ cb
      character cb*3
      character*1 cb1(3)
      equivalence (cb1,cb)

      character cbb*3
c
      logical ifonbc

      integer nfaces,nxyz,nel,ntot,ie,iface,isweep
      real bc1,bc2,bc3

      real tmpl1(lx1,ly1,lz1),tmpl2(lx1,ly1,lz1),tmpl3(lx1,ly1,lz1)
c
      ifonbc = .false.
c
      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dnekclock()

      nfaces=2*ldim
      nxyz  =lx1*ly1*lz1
      nel   =nelfld(ifield)
      ntot  =nxyz*nel

      call rzero(tmp1,ntot)
      call rzero(tmp2,ntot)
      call rzero(tmp3,ntot)

!     Velocity boundary conditions

      do 2100 isweep=1,2
         do 2000 ie=1,nel
         do 2000 iface=1,nfaces
            cb  = cbc(iface,ie,ifield)
            bc1 = bc(1,iface,ie,ifield)
            bc2 = bc(2,iface,ie,ifield)
            bc3 = bc(3,iface,ie,ifield)

            if (cb.eq.'V  ' .or. cb.eq.'VL ') then
!               prabal. WS, WSL not supported yet
!     $          cb.eq.'WS ' .or. cb.eq.'WSL') then
               call facev (tmp1,ie,iface,bc1,lx1,ly1,lz1)
               call facev (tmp2,ie,iface,bc2,lx1,ly1,lz1)
               call facev (tmp3,ie,iface,bc3,lx1,ly1,lz1)

!              False for V, VL  
!               if ( ifqinp(iface,ie) )
!     $         call globrot (tmp1(1,1,1,ie),tmp2(1,1,1,ie),
!     $                       tmp3(1,1,1,ie),ie,iface)
            endif

            if (cb.eq.'v  ' .or. cb.eq.'vl ' .or.
     $          cb.eq.'mv ' .or. cb.eq.'mvn' ) then
!               prabal. Below conditions not supported yet.
!     $          cb.eq.'ws ' .or. cb.eq.'wsl' .or.
!     $          cb1(1).eq.'d'.or.cb1(2).eq.'d'.or.cb1(3).eq.'d') then

                call faceiv (cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,lx1,ly1,lz1)

                if ( ifqinp(iface,ie) )
     $          call globrot (tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                        tmp3(1,1,1,ie),ie,iface)
            endif

!           In case of Symmetry, the amimuthal velocity needs to
!           be specified. tmpl1 and tmpl2 are placeholder arrays to be
!           discarded. The wall-normal component is zero (set via
!           application of masks) and the tangential component is
!           evaluated.
            if (cb.eq.'SYM') then
                call chcopy(cbb,cb,3)
                cb = '  d'
                call faceiv(cb,tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $            tmp3(1,1,1,ie),ie,iface,lx1,ly1,lz1)
                call chcopy(cb,cbb,3)
            endif

            IF (CB.EQ.'ON ' .OR. CB.EQ.'on ') then   ! 5/21/01 pff
                ifonbc =.true.
                call faceiv ('v  ',tmp1(1,1,1,ie),tmp2(1,1,1,ie),
     $                       tmp3(1,1,1,ie),ie,iface,lx1,ly1,lz1)
            ENDIF

 2000    CONTINUE
         do 2010 ie=1,nel
         do 2010 iface=1,nfaces
            if (cbc(iface,ie,ifield).eq.'W  ') then
               call facev (tmp1,ie,iface,0.0,lx1,ly1,lz1)
               call facev (tmp2,ie,iface,0.0,lx1,ly1,lz1)
               if (if3d) call facev (tmp3,ie,iface,0.0,lx1,ly1,lz1)
            endif
 2010    continue
C
C        Take care of Neumann-Dirichlet shared edges...
C
         if (isweep.eq.1) then
!            call opdsop(tmp1,tmp2,tmp3,'MXA')
           call dsop(tmp1,'MXA',lx1,ly1,lz1)
           call dsop(tmp2,'MXA',lx1,ly1,lz1)
           call dsop(tmp3,'MXA',lx1,ly1,lz1)
         else
!            call opdsop(tmp1,tmp2,tmp3,'MNA')
           call dsop(tmp1,'MNA',lx1,ly1,lz1)
           call dsop(tmp2,'MNA',lx1,ly1,lz1)
           call dsop(tmp3,'MNA',lx1,ly1,lz1)
         endif
 2100 CONTINUE
c
c     copy temporary array to velocity arrays.
c
      if ( .not.ifstrs ) then
         call col2(v1,mask1,ntot)
         call col2(v2,mask2,ntot)
         if (iff3d) call col2(v3,mask3,ntot)
         if (ifonbc) then
            call antimsk1(tmp1,mask1,ntot)
            call antimsk1(tmp2,mask2,ntot)
            if (iff3d) call antimsk1(tmp3,mask3,ntot)
         endif
      else
!         call rmask (v1,v2,v3,nelv)
!        rmask needs to be figured out properly           
         call col2_3(v1,v2,v3,v1mask,v2mask,v3mask,ntot)
      endif

      call add2_3(v1,v2,v3,tmp1,tmp2,tmp3,ntot)      
!      call add2(v1,tmp1,ntot)
!      call add2(v2,tmp2,ntot)
!      call add2(v3,tmp3,ntot)

!      if (ifneknekc) call fix_surface_flux

      tusbc=tusbc+(dnekclock()-etime1)

      return
      end subroutine bcdirvc_cyl
!-----------------------------------------------------------------------
      subroutine opgradt_f3d(outx,outy,outz,inpfld)

!     Compute DTx, DTy, DTz of an input field INPFLD
!     INPFLD is on Pressure grid        

      implicit none  

      include 'SIZE'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'F3D'

      real outx   (1)
      real outy   (1)
      real outz   (1)
      real inpfld (1)

      real const

      integer ntot1,ntot2

      real divv,bdivv
      common /scruz2/  divv (lx2,ly2,lz2,lelv)
     $ ,               bdivv(lx2,ly2,lz2,lelv)

      logical ifaxis_old

!      ifaxis = .true.
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

!     The ifaxis operators need to be correctly set inside cdtp
      ifaxis_old = ifaxis
      if (ifcyl_f3d) ifaxis = .true.

      call cdtp (outx,inpfld,rxm2,sxm2,txm2,1)
      call cdtp (outy,inpfld,rym2,sym2,tym2,2)


!     Note: This is done in cdtp if ifaxis==T      
!     Outz is used as a work array      
!     BM1*p/R
      if (ifcyl_f3d.and..not.ifaxis) then
!        ntot1 = lx1*ly1*lz1*nelv
!        call map21_all_f3d(outz,inpfld) 
!        call col2(outz,bm1,ntot1)
!        call invcol2(outz,ym1,ntot1)
!        call add2(outy,outz,ntot1)

        call copy(divv,inpfld,ntot2)
        call invcol2(divv,ym2,ntot2)
        call map21_weak(outz,divv)
        call add2(outy,outz,ntot1)
      endif  


!     In principle the third component comes from the imaginary part
!     Or from the real part if the other two components are imaginary      
!     I assume this is handled outside the routine      
!     BM1*p*dv/dtheta = BM1*k*p*v
!      call map21_all_f3d(outz,inpfld) 
!      call col2(outz,bm1,ntot1)            ! opgradt includes a mass matrix
!      const = -k_f3d
!      call cmult(outz,const,ntot1)
!      if (ifcyl_f3d) call invcol2(outz,ym1,ntot1)      ! 1/R 

      call copy(divv,inpfld,ntot2)
      const = -k_f3d
      call cmult(divv,const,ntot2)
      if (ifcyl_f3d) call invcol2(divv,ym2,ntot2)     ! 1/R
      call map21_weak(outz,divv)

      ifaxis = ifaxis_old

      return
      end subroutine opgradt_f3d
!-----------------------------------------------------------------------
      
      subroutine opdiv_f3d(outfld,inx,iny,inz)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MASS'
      include 'F3D'

      real outfld (1)   ! Pressure Mesh
      real inx    (1)   ! Vel. Mesh   
      real iny    (1)   ! Vel. Mesh
      real inz    (1)   ! Vel. Mesh

      real dummy
      common /scrcg2/ dummy(lx1*ly1*lz1*lelt)    ! In principle I only need
                                                ! the pressure sized mesh
      integer ntot2
      logical ifaxis_old      

      ntot2 = lx2*ly2*lz2*nelv

      ifaxis_old = ifaxis

!     2D Divergence
!     The ifaxis operators need to be correctly set 
!     inside multd (called in opdiv)
      if (ifcyl_f3d) ifaxis = .true.
     
      call opdiv  (outfld,inx,iny,inz)

!     1/R*B*(dp/dR)
      if (ifcyl_f3d.and..not.ifaxis) then
!       We calculate this term ourself            
        call map12_all_f3d(dummy,iny)
        call invcol2(dummy,ym2,ntot2)
!        call chsign(dummy,ntot2)
        call Xaddcol3(outfld,dummy,bm2,ntot2)
      endif        
    
!     Map third component to pressure grid 
      call map12_all_f3d(dummy,inz)
      call cmult(dummy,k_f3d,ntot2)
      if (ifcyl_f3d) call invcol2(dummy,ym2,ntot2)     ! 1/R
      call Xaddcol3(outfld,dummy,bm2,ntot2)

      ifaxis = ifaxis_old

      return
      end subroutine opdiv_f3d        
!-----------------------------------------------------------------------

      subroutine chkdiv_f3d

!     Check Divergence for the complex variable case 
      
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

      real divv,bdivv
      common /scruz/ divv (lx2,ly2,lz2,lelv)
     $ ,             bdivv(lx2,ly2,lz2,lelv)

      real dummy1,dummy2
      common /screv/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt)    

      integer jpr,jpi
      integer ntot1,ntot2

      real glsc2
      real dnorm

      jpr = 1
      jpi = 2

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

!     Real
      call cmult2(dummy1,vzp(1,jpi),-1.0,ntot1)
      call opdiv_f3d(bdivv,vxp(1,jpr),vyp(1,jpr),dummy1)
      call col3 (divv,bdivv,bm2inv,ntot2)
      dnorm = sqrt(glsc2(divv,bdivv,ntot2)/volvm2) 
      if (nio.eq.0) write (6,*) istep,' Real: Dnorm', dnorm

!     Imaginary
      call cmult2(dummy1,vzp(1,jpr),1.0,ntot1)
      call opdiv_f3d(bdivv,vxp(1,jpi),vyp(1,jpi),dummy1)
      call col3 (divv,bdivv,bm2inv,ntot2)
      dnorm = sqrt(glsc2(divv,bdivv,ntot2)/volvm2) 
      if (nio.eq.0) write (6,*) istep,' Imaginary: Dnorm', dnorm

      return
      end subroutine chkdiv_f3d
!---------------------------------------------------------------------- 

      subroutine map21_weak(y,x)

!     Map X from mesh M2 to mesh M1 (Y)
      
      implicit none

      include 'SIZE'
!      include 'GEOM'
      include 'MASS'
      include 'IXYZ'

      include 'F3D'

      real x(lx2,ly2,lz2,lelv)
      real y(lx1,ly1,lz1,lelv)
      real wk1,wk2
      common /screv/ wk1(lx1,ly1,lz1,lelt),
     $               wk2(lx1,ly1,lz1,lelt)    

      integer e,nxyz1,nxyz2

      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      do e=1,nelv
       if (ifcyl_f3d) call ifaxisop_f3d(e)
        call col3    (y(1,1,1,e),x(1,1,1,e),bm2(1,1,1,e),nxyz2)
!       I assume 1/R is already factored in x            
!        call invcol2 (y(1,1,1,e),ym2(1,1,1,e),nxyz2)
        call mxm     (ixtm12,lx1,y(1,1,1,e),lx2,wk1(1,1,1,e),ly2)
        call mxm     (wk1(1,1,1,e),lx1,iym12,ly2,y(1,1,1,e),ly1)
      enddo

      return
      end subroutine map21_weak

!-----------------------------------------------------------------------
      subroutine esolver_f3d(res,h1,h2,h2inv,intype)

!     Choose E-solver
!--------------------------------------------------------------------

      implicit none
        
      include 'SIZE'
      include 'ESOLV'
      include 'INPUT'
      include 'F3D'

      real res   (lx2,ly2,lz2,lelv)
      real h1    (lx1,ly1,lz1,lelv)
      real h2    (lx1,ly1,lz1,lelv)
      real h2inv (lx1,ly1,lz1,lelv)

!      common /scruz/ wk1(lx2*ly2*lz2*lelv)
!     $             , wk2(lx2*ly2*lz2*lelv)
!     $             , wk3(lx2*ly2*lz2*lelv)

      include 'CTIMER'
      real kwave2
      real icg
      integer intype

      if (icalld.eq.0) teslv=0.0

      if (abs(k_f3d).lt.1.0e-12) then
        call ortho(res) !Ensure that residual is orthogonal to null space
      endif  

      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()

      if (.not. ifsplit) then
         if (param(42).eq.1) then
            call uzawa_f3d (res,h1,h2,h2inv,intype,icg)
         else
            call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
         endif
      else
         write(6,*) 'ERROR: E-solver does not exist PnPn'
         call exitt
      endif

      teslv=teslv+(dnekclock()-etime1)

      return
      end subroutine esolver_f3d
!-----------------------------------------------------------------------
      subroutine uzawa_f3d(rcg,h1,h2,h2inv,intype,iter)
!-----------------------------------------------------------------------
!
!     Solve the pressure equation by (nested) preconditioned 
!     conjugate gradient iteration.
!     INTYPE =  0  (steady)
!     INTYPE =  1  (explicit)
!     INTYPE = -1  (implicit)
!
!-----------------------------------------------------------------------

      implicit none 

      include 'SIZE'
!      include 'MASS'
      include 'INPUT'         ! PARAM
      include 'TSTEP'         ! ISTEP
      include 'F3D'
!      include 'TOTAL'

      real divex
      common  /ctolpr/ divex
      
      logical          ifprint
      common  /cprint/ ifprint

      real             rcg  (lx2,ly2,lz2,lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      real wp,xcg,pcg,rpcg
      common /scruz/   wp   (lx2,ly2,lz2,lelv)
     $ ,               xcg  (lx2,ly2,lz2,lelv)
     $ ,               pcg  (lx2,ly2,lz2,lelv) 
     $ ,               rpcg (lx2,ly2,lz2,lelv)
 
      real*8 etime1,dnekclock
      integer*8 ntotg,nxyz2

      integer iter,iconv,ntot1,ntot2,intype,nelgv
      real rrp1,rrp2,beta,div0,tolpss,rnorm,ratio,pap,alpha
      real pcgmx,wp_mx,h1_mx,h2_mx,rnrm1,rrpx

      real glsc2,glamax
      

      etime1 = dnekclock()
      DIVEX = 0.
      ITER  = 0
c
      CALL CHKTCG2 (TOLPS,RCG,ICONV)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   TOLPS = abs(param(21))
C
c      IF (ICONV.EQ.1) THEN
c         IF (NID.EQ.0) WRITE(6,9999) ITER,DIVEX,TOLPS
c         return
c      ENDIF

      nxyz2 = lx2*ly2*lz2
      ntot2 = nxyz2*nelv
      ntotg = nxyz2*nelgv

      CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)

      RRP1 = GLSC2 (RPCG,RCG,NTOT2)
      CALL COPY    (PCG,RPCG,NTOT2)
      CALL RZERO   (XCG,NTOT2)
      if (rrp1.eq.0) return
      BETA = 0.
      div0=0.
C
      tolpss = tolps
      DO 1000 ITER=1,4000 !NMXP
C
C        CALL CONVPR  (RCG,tolpss,ICONV,RNORM)
         call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)

         if (iter.eq.1)      div0   = rnorm
         if (param(21).lt.0) tolpss = abs(param(21))*div0

         ratio = rnorm/div0
         IF (IFPRINT.AND.NIO.EQ.0) 
     $   WRITE (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66    format(i5,1p4e12.5,i8,' Divergence')
c
         IF (ICONV.EQ.1.and.iter.gt.1) GOTO 9000
c        IF (ICONV.EQ.1.and.(iter.gt.1.or.istep.le.2)) GOTO 9000
c        IF (ICONV.EQ.1) GOTO 9000
c        if (ratio.le.1.e-5) goto 9000


         IF (ITER .NE. 1) THEN
            BETA = RRP1/RRP2
            CALL ADD2S1 (PCG,RPCG,BETA,NTOT2)
         ENDIF

         call cdabdtp_f3d  (wp,pcg,h1,h2,h2inv,intype)
         PAP   = GLSC2 (PCG,WP,NTOT2)
         IF (PAP.NE.0.) THEN
            ALPHA = RRP1/PAP
         ELSE
            pcgmx = glamax(pcg,ntot2)
            wp_mx = glamax(wp ,ntot2)
            ntot1 = lx1*ly1*lz1*nelv
            h1_mx = glamax(h1 ,ntot1)
            h2_mx = glamax(h2 ,ntot1)
            if (nid.eq.0) write(6,*) 'ERROR: pap=0 in uzawa.'
     $      ,iter,pcgmx,wp_mx,h1_mx,h2_mx
            call exitt
         ENDIF
         CALL ADD2S2 (XCG,PCG,ALPHA,NTOT2)
         CALL ADD2S2 (RCG,WP,-ALPHA,NTOT2)

         if (iter.eq.-1) then
            call convprn (iconv,rnrm1,rrpx,rcg,rpcg,tolpss)
            if (iconv.eq.1) then
               rnorm = rnrm1
               ratio = rnrm1/div0
               if (nio.eq.0) 
     $         write (6,66) iter,tolpss,rnrm1,div0,ratio,istep
               goto 9000
            endif
         endif

!        prabal            
!         call ortho(rcg)
         if (abs(k_f3d).lt.1.0e-12) then
           call ortho(rcg)
         endif  

         RRP2 = RRP1
         CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
c        RRP1 = GLSC2 (RPCG,RCG,NTOT2)
!     prabal
!      call copy(rpcg,rcg,ntot2)



 1000 CONTINUE
      if (nid.eq.0) WRITE (6,3001) ITER,RNORM,tolpss
c     if (istep.gt.20) CALL EMERXIT
 3001 FORMAT(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
 9000 CONTINUE

      divex = rnorm
      iter  = iter-1

      if (iter.gt.0) call copy (rcg,xcg,ntot2)
!      call ortho(rcg)
      if (abs(k_f3d).lt.1.0e-12) then
        call ortho(rcg)
      endif  

      etime1 = dnekclock()-etime1
      IF (NIO.EQ.0) WRITE(6,9999) ISTEP, '  U-Press std. ',
     &                            ITER,DIVEX,div0,tolpss,etime1
 9999 FORMAT(I11,a,I7,1p4E13.4)
19999 FORMAT(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)
C
C
      return
      end subroutine uzawa_f3d
!-----------------------------------------------------------------------

      subroutine ifaxisop_f3d(e)

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'INPUT'
      include 'GEOM'    ! ifrzer

      include 'F3D'

      integer nxyz2,ly12
      integer e               ! element no

      nxyz2 = lx2*ly2*lz2

      ly12   = ly1*ly2
      if (ifrzer(e)) then
         call copy (iym12,iam12,ly12)
         call copy (dym12,dam12,ly12)
         call copy (w3m2,w2am2,nxyz2)
      else
         call copy (iym12,icm12,ly12)
         call copy (dym12,dcm12,ly12)
         call copy (w3m2,w2cm2,nxyz2)
      endif


      return
      end subroutine ifaxisop_f3d        
!---------------------------------------------------------------------- 








