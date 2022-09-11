!> @file arn_usr.f
!! @ingroup arn_arp
!! @brief Set of user defined subroutines for Arpack
!! @author Prabal Negi
!! @date Sept 10, 2022
!
!====================================================================== 

      subroutine arn_nektoworkda()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      integer i

#ifdef ARPACK_DIRECT
      ! no temperature here
      ! A*x = lambda*x
      ! velocity
      if (arna_ifcomplex) then
        i = ipntarp(2)
        call copytocomplex(workda(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(workda(i),v1mask,tst_nv)
        i = i + tst_nv
        call copytocomplex(workda(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(workda(i),v2mask,tst_nv)
        i = i + tst_nv
        if (iff3d) then
          call copytocomplex(workda(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(workda(i),v3mask,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copytocomplex(workda(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif  
      else 
        i = ipntarp(2) 
        call col3(workda(i),VXP,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(workda(i),VYP,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then 
          call col3(workda(i),VZP,V3MASK,tst_nv)
          i = i + tst_nv
        endif
!       pressure
        if (arna_ifpr) then
          call copy(workda(i),prp,tst_np)        
          i = i + tst_np
        endif  
      endif  
#else
      ! velocity
      ! A*x = lambda*M*x
      if (arna_ifcomplex) then
        i = ipntarp(2)
        call copytocomplex(workda(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(workda(i),v1mask,tst_nv)
        i = i + tst_nv
        call copytocomplex(workda(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(workda(i),v2mask,tst_nv)
        i = i + tst_nv
        if (iff3d) then
          call copytocomplex(workda(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(workda(i),v3mask,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copytocomplex(workda(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif
        if (ifheat) then
          call copytocomplex(workda(i),tp(1,1,1),tp(1,1,2),tst_nt)
          call col2_cr(workda(i),tmask,tst_nt)
          i = i + tst_nt
        endif
      else
        i = ipntarp(2) 
        call col3(workda(i),VXP,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(workda(i),VYP,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then
          call col3(workda(i),VZP,V3MASK,tst_nv)
          i = i + tst_nv
        endif
!       pressure
        if (arna_ifpr) then
          call copy(workda(i),prp,tst_np)        
          i = i + tst_np
        endif  
!       temperature
        if (IFHEAT) then
          call col3(workda(i),TP,TMASK,tst_nt)
          i = i + tst_nt
        endif  
        ! this may be not necessary, but ARPACK manual is not clear about it
        !call col3(workda(ipntarp(1)),VXP,BM1,tst_nv)
        !call col3(workda(ipntarp(1)+tst_nv),VYP,BM1,tst_nv)
        !if (IF3D) call col3(workda(ipntarp(1)+2*tst_nv),VZP,BM1,tst_nv)
      endif       ! arna_ifcomplex 
#endif

      return
      end subroutine arn_nektoworkda
!---------------------------------------------------------------------- 

      subroutine arn_innerprod()

      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'MASS'            ! BM1
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      integer i,ix,iy

!     From Documentation:         
!     IDO =  2: compute  Y = M * X  where
!               IPNTR(1) is the pointer into WORKD for X,
!               IPNTR(2) is the pointer into WORKD for Y.

      if (arna_ifcomplex) then
        ix = ipntarp(1)
        iy = ipntarp(2) 
!       velocities 
        call ccopy(workda(iy),workda(ix),arna_ns)
        i = iy
        call col2_cr(workda(i),bm1,tst_nv)
        call col2_cr(workda(i),v1mask,tst_nv)
        i = i+tst_nv
        call col2_cr(workda(i),bm1,tst_nv)
        call col2_cr(workda(i),v2mask,tst_nv)
        i = i+tst_nv
        if (iff3d) then
          call col2_cr(workda(i),bm1,tst_nv)
          call col2_cr(workda(i),v3mask,tst_nv)
          i = i+tst_nv
        endif
!       pressure                  
        if (arna_ifpr) then
          call czero(workda(i),tst_np)
          i = i+tst_np
        endif
!       temperature
        if (ifheat) then
          call col2_cr(workda(i),bm1,tst_nt)
          call col2_cr(workda(i),tmask,tst_nt)
!         rescale coefficient of temperature
          call cmult_cr(workda(i),1.0,tst_nt)
          i = i+tst_nt
        endif
      else
        ix = ipntarp(1)
        iy = ipntarp(2) 
        i  = iy
!       velocity
        call col3(workda(i),BM1,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(workda(i),BM1,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then
          call col3(workda(i),BM1,V3MASK,tst_nv)
          i = i + tst_nv
        endif
!       Pressure             
        if (arna_ifpr) then
          call rzero(workda(i),tst_np)
!          call copy(workda(i),bm2,tst_np)       ! prabal
          i = i + tst_np
        endif     

!       Temperature
        if(IFHEAT) then
           call col3(workda(i),BM1,TMASK,tst_nt)
           i = i + tst_nt
!          Temperature coefficients
           call cht_weight_fun (workda(ipntarp(2)),
     $          workda(ipntarp(2)+tst_nv),
     $          workda(ipntarp(2)+2*tst_nv),
     $          workda(ipntarp(2)+NDIM*tst_nv),1.0)
        endif

        call col2(workda(iy),workda(ix),arna_ns)
      endif        ! arna_ifcomplex

      return
      end subroutine arn_innerprod
!---------------------------------------------------------------------- 

      subroutine arn_workdatonek()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      integer i

!     move renormed data back to nekton
      if (arna_ifcomplex) then
        i = ipntarp(1)
        call copytoreal(vxp(1,1),vxp(1,2),workda(i),tst_nv)
        i = i+tst_nv
        call copytoreal(vyp(1,1),vyp(1,2),workda(i),tst_nv)
        i = i+tst_nv
        if (iff3d) then 
          call copytoreal(vzp(1,1),vzp(1,2),workda(i),tst_nv)
          i = i+tst_nv
        endif
        if (arna_ifpr) then 
          call copytoreal(prp(1,1),prp(1,2),workda(i),tst_np)
          i = i+tst_np
        endif
        if (ifheat) then
          call copytoreal(tp(1,1,1),tp(1,1,2),workda(i),tst_nt)
          i = i+tst_nt
        endif  
      else  
        ! velocity
        i = ipntarp(1)
        call copy(VXP,workda(i),tst_nv)
        i = i + tst_nv
        call copy(VYP,workda(i),tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then 
          call copy(VZP,workda(i),tst_nv)
          i = i + tst_nv
        endif
!       Pressure 
        if (arna_ifpr) then
          call copy(prp,workda(i),tst_np)
          i = i + tst_np
        endif  
!       Temperature
        if (IFHEAT) then
          call copy(TP,workda(i),tst_nt)
          i = i + tst_nt
        endif  

      endif          ! arna_ifcomplex

!     make sure the velocity and temperature fields are continuous at
!     element faces and edges
      call tst_dssum        ! Should not be there.


      return
      end subroutine arn_workdatonek
!---------------------------------------------------------------------- 

      subroutine arn_nektoresida()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'ARN_ARPD'
      include 'TSTEPPERD'

      include 'F3D'

      integer i

      ! if no restart fill RESIDA with initial conditions
      ! V?MASK removes points at the wall and inflow
#ifdef ARPACK_DIRECT
      ! A*x = lambda*x
      ! velocity
      if (arna_ifcomplex) then
        i = 1
        call copytocomplex(resida(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(resida(i),v1mask,tst_nv) 
        i = i + tst_nv

        call copytocomplex(resida(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(resida(i),v2mask,tst_nv) 
        i = i + tst_nv
        
        if (iff3d) then
          call copytocomplex(resida(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(resida(i),v3mask,tst_nv)
          i = i + tst_nv
        endif  
         
        if (arna_ifpr) then
          call copytocomplex(resida(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif  
      else           ! real arithmetic
        i = 1
        call col3(resida(i),VXP,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(resida(i),VYP,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then
          call col3(resida(i),VZP,V3MASK,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then 
          call copy(resida(i),prp,tst_np)
          i = i + tst_np
        endif  
      endif          ! arna_ifcomplex             

      ! no temperature here
#else
      ! A*x = lambda*M*x
      if (arna_ifcomplex) then   
        i = 1
        call copytocomplex(resida(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(resida(i),v1mask,tst_nv) 
        i = i + tst_nv

        call copytocomplex(resida(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(resida(i),v2mask,tst_nv) 
        i = i + tst_nv
        
        if (iff3d) then
          call copytocomplex(resida(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(resida(i),v3mask,tst_nv)
          i = i + tst_nv
        endif  
         
        if (arna_ifpr) then
          call copytocomplex(resida(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif

        if (ifheat) then
          call copytocomplex(resida(i),tp(1,1,1),tp(1,1,2),tst_nt)
          call col2_cr(resida(i),tmask,tst_nt)
          i = i + tst_nt
        endif 

      else           ! real arithmetic
        ! velocity
        i = 1
        call col3(resida(i),VXP,v1mask,tst_nv)
        i = i + tst_nv
        call col3(resida(i),VYP,v2mask,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then
          call col3(resida(i),VZP,v3mask,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copy(resida(i),PRP,tst_np)
          i = i + tst_np
        endif  
        ! temperature
        if (IFHEAT) then
          call col3(resida(i),TP,tmask,tst_nt)
          i = i + tst_nt
        endif  
      endif          ! arna_ifcomplex

#endif

      return
      end subroutine arn_nektoresida
!---------------------------------------------------------------------- 

