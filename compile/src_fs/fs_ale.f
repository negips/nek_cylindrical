!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for free surface mesh movement.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine gen_mapping_mvb

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'

      include 'FS_ALE'

      integer i,j,k,n,n1

      integer igl_running_sum       ! function
      integer iglsum

      real xmean,ymean,vol
      character cb*3
      real temp1(lx1,ly1,lz1,lelv)
      integer kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface,ie,ieg
      integer ix,iy,iz,nfaces
      real vlsum
      integer ifield, nxyz
      real temp0(lx1,ly1,lz1)
      character str*132

      integer gl_fs_nel

      logical outflds
      real dmy1,dmy2,dmy3
      common /scrsf/ dmy1(lx1,ly1,lz1,lelt),
     $               dmy2(lx1,ly1,lz1,lelt),
     $               dmy3(lx1,ly1,lz1,lelt)

      real sint,sarea
      real AW_int             ! Air water interface
      real tol
      logical ifcrnr
      logical cornerpoint     ! function


      AW_int  = 1.23
      tol     = 1.0e-6
      outflds = .true.

      call mntr_log(fs_id,fs_log,'Generating Global Mapping.')

      n = lx1*ly1*lz1*nelv

!     Keep original mesh location      
      call copy(xm0_fs,xm1,n)
      call copy(ym0_fs,ym1,n)
      if (if3d) call copy(zm0_fs,zm1,n)


!     Create Global numbering and GS handle
!      call fs_map()
      call fs_map1()

      call rzero(fs_mask,n)

!     Multiplicity (Inverse) (vmult)
      call rzero(fs_vmult,n)
      nfaces = 2*ndim
      ifield = 1
      nxyz   = lx1*ly1*lz1
      nx     = lx1
      ny     = ly1
      nz     = lz1
      fs_nel = 0
      do ie=1,nelv
        do iface=1,nfaces
          cb  = cbc(iface,ie,ifield)
          ieg = lglel(ie)
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)

          fs_cbc(iface,ie) = '   '

          call surface_int(sint,sarea,xm1,ie,iface)
          sint = sint/sarea
!         Change these conditions for actual case            
          if (cb.eq.'O  ') then
!          if (cb.eq.'E  '.and.(abs(sint-AW_int).lt.tol)) then
             fs_nel = fs_nel + 1
             fs_cbc(iface,ie) = 'INT'
             fs_elno(fs_nel)  = ie
             fs_iface(fs_nel) = iface
             do iz=kz1,kz2
             do iy=ky1,ky2
             do ix=kx1,kx2
               fs_vmult(ix,iy,iz,ie) = 1.0
               fs_mask(ix,iy,iz,ie)  = 1.0
               ifcrnr = cornerpoint(ix,iy,iz)
               if (ifcrnr) fs_crmask(ix,iy,iz,ie) = 1.0
             enddo
             enddo
             enddo
          endif                     
        enddo
      enddo      

      gl_fs_nel = iglsum(fs_nel,1)
      call blank(str,132)
      write(str,'(A23,1x,I5)') 'No. Interface Elements:', gl_fs_nel
      call mntr_log(fs_id,fs_log,str)

      call fgslib_gs_op(fs_gs_handle,fs_vmult,1,1,0)  ! 1 ==> +
      call invcol1(fs_vmult,n)

!      call fgslib_gs_op(fs_gs_handle,temp1,1,1,0)  ! 1 ==> +
!      call col2(temp1,fs_vmult,n)

      call fs_symo_save       ! save SYM/O corner points

      call fs_gen_damping

      if (outflds) then
        do i=1,n
          dmy1(i,1,1,1) = fs_gl_num(i)+0.0
        enddo
        ifto = .true.
        call outpost(dmy1,fs_mask,fs_vmult,pr,fs_damp,'ale')
!        call outpost(dmy1,fs_mask,fs_vmult,pr,fs_vmult,'ale')
      endif  

      call mntr_log(fs_id,fs_log,'Global Mapping Done.')

      return
      end subroutine gen_mapping_mvb
!----------------------------------------------------------------------
      subroutine fs_map

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'

      include 'FS_ALE'

!     Temporary arrays
      real wk1,wk2,wk3
      common /scrns/ wk1(lx1*ly1*lz1*lelt*2),
     $               wk2(lx1*ly1*lz1*lelt*2),
     $               wk3(lx1*ly1*lz1*lelt*3)

      real ta1,ta2,ta3
      common /scrsf/ ta1 (lx1*ly1*lz1*lelv)
     $ ,             ta2 (lx1*ly1*lz1*lelv)
     $ ,             ta3 (lx1*ly1*lz1*lelv)

      real tb1,tb2,tb3,tb4
      common /scruz/ tb1 (lx1*ly1*lz1*lelv)
     $ ,             tb2 (lx1*ly1*lz1*lelv)
     $ ,             tb3 (lx1*ly1*lz1*lelv)
     $ ,             tb4 (lx1*ly1*lz1*lelv) 

      real x1local_u(lx1*lz1*lelv),x2local_u(lx1*lz1*lelv)
      integer local_num(lx1*ly1*lz1*lelv) ! local numbering
      integer local_num_u(lx1*lz1*lelv)   ! unique local numbering
      integer sortind1(lx1*ly1*lz1*lelv)
      integer sortind2(lx1*ly1*lz1*lelv)
      integer i,j,k,n,n1
      integer jloop

      integer nunq,nunq2                  ! Unique local points
      integer n2,n3
      integer glnos
      integer own
      integer nown,unown                  ! unique owned/not-owned points
     
      real tol

      integer igl_running_sum       ! function
      integer iglsum

      integer key(lx1*ly1*lz1*lelv)
      integer ind(lx1*ly1*lz1*lelv)
      integer ind2(lx1*ly1*lz1*lelv)
      integer ind3(lx1*ly1*lz1*lelv)
      integer ninseg(lx1*ly1*lz1*lelv)
      logical ifseg(lx1*ly1*lz1*lelv)
      integer nkey

      integer iseg,nseg,nsort

      integer glnum(lx1*lz1*lelv)
      integer glnum_wk(lx1*lz1*lelv*2)
      integer owner(lx1*lz1*lelv)
      integer glnum2(lx1*lz1*lelv)
      integer owner_wk(lx1*lz1*lelv*2)

      integer cnsort
      integer csteps,ic,ip,ipass,il
      integer msg_id1,msg_id2,msg_id3,msg_id4,msg_id5
      integer srcid,dstid,len
      integer ntot1

      integer irecv
      integer ivlmin,ivlmax

      ntot1 = lx1*ly1*lz1*nelv

      jloop = 1
      if (if3d) jloop = 2

      call rzero3(ta1,ta2,ta3,ntot1)

!     Segments (initialization)        
      do i=1,ntot1
        ifseg(i) = .false.
      enddo  

      ifseg(1)  = .true.
      ninseg(1) = ntot1
      nseg      = 1
      tol       = 1.0e-14

      call copy(ta1,ym1,ntot1)
      call copy(ta2,zm1,ntot1)
      do i=1,ntot1
        ind(i) = i
      enddo  

!     Local sorting
      do ipass=1,2
        do j=1,jloop              ! For 2D we only do one dimension
          nkey = 1
          i    = 1
          do iseg = 1,nseg
            if (j.eq.1) then
              call tuple_sort(ta1(i),1,ninseg(iseg),1,nkey,ind2,ta3)
              call iswap_ip(ind(i),ind2,ninseg(iseg))
              call swap_ip(ta2(i),ind2,ninseg(iseg))
            elseif (j.eq.2) then
              call tuple_sort(ta2(i),1,ninseg(iseg),1,nkey,ind2,ta3)
              call iswap_ip(ind(i),ind2,ninseg(iseg))
              call swap_ip(ta1(i),ind2,ninseg(iseg))
            endif
            i = i + ninseg(iseg) 
          enddo

          do i=2,ntot1
            if (j.eq.1) then
              if ((abs(ta1(i)-ta1(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif  
            elseif (j.eq.2) then
              if ((abs(ta2(i)-ta2(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif    
            endif
          enddo

!         Count up number of different segments
          nseg = 0
          do i=1,ntot1
             if (ifseg(i)) then
                nseg = nseg+1
                ninseg(nseg) = 1
                x1local_u(nseg) = ta1(i)
                x2local_u(nseg) = ta2(i)
                local_num_u(nseg) = nseg  
             else
                ninseg(nseg) = ninseg(nseg) + 1
             endif
             local_num(i) = nseg
          enddo
        enddo     ! j
      enddo       ! ipass

      nunq = nseg  

      n2 = igl_running_sum(nunq)
      n1 = n2 - nunq
      n3 = iglsum(nunq,1)

!!    Cycle 1: Globally Unique points and Ownership
!--------------------------------------------------       
!     Set ownership to current processor      
      call ifill(owner,nid,nunq) 
      call izero(owner_wk,lx1*lz1*lelv*2) 
      call icopy(owner_wk,owner,nunq)

!     First global communication steps to resolve global ownership.
!     This only serves to mark globally unique points.
!     Second global communication is then done for global numbering.
!     Third global communication is done to number the un-owned points.
      

!     Store local unique values in wk1,wk2      
      call copy(wk1(1),x1local_u,nunq)
      call copy(wk2(1),x2local_u,nunq)

      call copy(tb1(1),x1local_u,nunq)
      call copy(tb2(1),x2local_u,nunq)

      csteps = int(log(np+0.0)/log(2.0))
      if (np.gt.2**csteps) csteps = csteps+1

      ifseg(1)  = .true.
      ninseg(1) = nunq
      nseg      = 1
      tol       = 1.0e-14

!     work array size 
      do ic=1,csteps

        il = 2**(ic-1)
        srcid = nid - il
        dstid = nid + il
        if (srcid.lt.0)  srcid = srcid + np
        if (dstid.ge.np) dstid = dstid - np

!       Put local coordinates into work array 
        call copy(wk1(1),tb1,nunq)
        call copy(wk2(1),tb2,nunq)
        call icopy(owner_wk,owner,nunq)

!       set buffer for the number of elements to receive
        len = isize
        msg_id1 = irecv(0,cnsort,len)

!       send local size of array
        call csend(0,nunq,len,dstid,0)

!       wait for length data
        call msgwait(msg_id1)

!       receive x1 data
        len = wdsize*cnsort
        msg_id2 = irecv(1,wk1(nunq+1),len)

!       receive x2 data
        len = wdsize*cnsort
        msg_id3 = irecv(2,wk2(nunq+1),len)

!       receive ownership data
        len = isize*cnsort
        msg_id4 = irecv(3,owner_wk(nunq+1),len)

!       Send data
        len = wdsize*nunq
        call csend(1,tb1,len,dstid,0)
        call csend(2,tb2,len,dstid,0)
        len = isize*nunq
        call csend(3,owner,len,dstid,0)

!       finish communication            
        call msgwait(msg_id2)
        call msgwait(msg_id3)
        call msgwait(msg_id4)

        ifseg(1)  = .true.
        nsort     = nunq + cnsort
        ninseg(1) = nsort
        nseg      = 1
        tol       = 1.0e-14
        ind2(1)   = 1
        do i=2,nsort
          ifseg(i) = .false.
          ind2(i)  = i              ! original array possitions
        enddo  

!       Perform local sorting
        do ip=1,2       ! pass
          do j=1,jloop      ! dimensions of sorting data
            nkey = 1
            i    = 1
            do iseg = 1,nseg
              if (j.eq.1) then
                call tuple_sort(wk1(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call iswap_ip(owner_wk(i),ind,ninseg(iseg))
                call swap_ip(wk2(i),ind,ninseg(iseg))
                call iswap_ip(ind2(i),ind,ninseg(iseg))
              elseif (j.eq.2) then
                call tuple_sort(wk2(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call iswap_ip(owner_wk(i),ind,ninseg(iseg))
                call swap_ip(wk1(i),ind,ninseg(iseg))
                call iswap_ip(ind2(i),ind,ninseg(iseg))
              endif
              i = i + ninseg(iseg) 
            enddo

            do i=2,nsort
              if (j.eq.1) then
                if ((abs(wk1(i)-wk1(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif
              elseif (j.eq.2) then
                if ((abs(wk2(i)-wk2(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif 
              endif
            enddo

!           Count up number of different segments
            nseg = 0
            do i=1,nsort
              if (ifseg(i)) then
                 nseg = nseg+1
                 ninseg(nseg) = 1
              else
                 ninseg(nseg) = ninseg(nseg) + 1
              endif
            enddo
          enddo     ! j
        enddo       ! ip

!       Reassign ownership
        i = 1
        do iseg=1,nseg
          own   = ivlmin(owner_wk(i),ninseg(iseg))
          call ifill(owner_wk(i),own,ninseg(iseg))
          i = i + ninseg(iseg)
        enddo 

!       Reassign local ownership            
        do i=1,nsort
          j = ind2(i)
          if (j.le.nunq) then       ! local point
            owner(j) = owner_wk(i)
          endif
        enddo  

      enddo         ! ic

!!     debugging
!      call nekgsync()
!      if (nid.eq.0) then
!        do i=1,nunq
!          write(6,*), nid,owner(i),tb1(i),tb2(i),nseg
!        enddo
!      endif
!      call nekgsync()
!
!      if (nid.eq.1) then
!        do i=1,nunq
!          write(6,*), nid,owner(i),tb1(i),tb2(i),nseg
!        enddo
!      endif
!      call nekgsync()


!!    Cycle 2: Global Numbering of Owned points.
!-------------------------------------------------- 
      j = 0
      do i=1,nunq
        if (owner(i).eq.nid) then
          j = j+1
          tb1(j)=x1local_u(i)
          tb2(j)=x2local_u(i)
        endif
      enddo
      nunq2 = j                     ! Unique owned elements

      call rzero(wk1,ntot1*2)
      call rzero(wk2,ntot1*2)

      call copy(wk1,tb1,nunq2)
      call copy(wk2,tb2,nunq2)
      call izero(ind3,ntot1)

      do i=1,nunq2
        glnum(i) = i
      enddo
      call izero(glnum_wk,nunq2)

      ifseg(1)  = .true.
      ninseg(1) = nunq2
      nseg      = 1
      tol       = 1.0e-14

      do ic=1,csteps

        il = 2**(ic-1)
        srcid = nid - il
        dstid = nid + il
        if (srcid.lt.0)  srcid = srcid + np
        if (dstid.ge.np) dstid = dstid - np

!       set buffer for the number of elements to receive
        len = isize
        msg_id1 = irecv(0,cnsort,len)

!       send local size of array
        call csend(0,nunq2,len,dstid,0)

!       wait for length data
        call msgwait(msg_id1)

!       Put local coordinates into work array 
        call copy(wk1(1),tb1,nunq2)
        call copy(wk2(1),tb2,nunq2)

!       receive x1 data
        len = wdsize*cnsort
        msg_id2 = irecv(1,wk1(nunq2+1),len)

!       receive x2 data
        len = wdsize*cnsort
        msg_id3 = irecv(2,wk2(nunq2+1),len)

!       Send data
        len = wdsize*nunq2
        call csend(1,tb1(1),len,dstid,0)
        call csend(2,tb2(1),len,dstid,0)

!       finish communication            
        call msgwait(msg_id2)
        call msgwait(msg_id3)

        ifseg(1)  = .true.
        nsort     = nunq2 + cnsort
        ninseg(1) = nsort
        nseg      = 1
        tol       = 1.0e-14

        ind3(1)   = 1
        do i=2,nsort
          ifseg(i) = .false.
          ind3(i)  = i              ! original array possitions
        enddo  

!       Perform local sorting
        do ip=1,2       ! pass
          do j=1,jloop  ! dimensions of sorting data
            nkey = 1
            i    = 1
            do iseg = 1,nseg
              if (j.eq.1) then
                call tuple_sort(wk1(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call swap_ip(wk2(i),ind,ninseg(iseg))
                call iswap_ip(ind3(i),ind,ninseg(iseg))
              elseif (j.eq.2) then
                call tuple_sort(wk2(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call swap_ip(wk1(i),ind,ninseg(iseg))
                call iswap_ip(ind3(i),ind,ninseg(iseg))
              endif
              i = i + ninseg(iseg) 
            enddo

            do i=2,nsort
              if (j.eq.1) then
                if ((abs(wk1(i)-wk1(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif
              elseif (j.eq.2) then
                if ((abs(wk2(i)-wk2(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif 
              endif
            enddo

!           Count up number of different segments
!           We should have as many segments as points.            
            nseg = 0
            do i=1,nsort
              if (ifseg(i)) then
                 nseg = nseg+1
                 ninseg(nseg) = 1
              else
                 ninseg(nseg) = ninseg(nseg) + 1
              endif
            enddo
          enddo     ! j
        enddo       ! ip

!       Adjust global numbering
        do i=1,nsort
          j = ind3(i)   ! position in local array
          if (j.le.nunq2) then       ! local point
            glnum_wk(j) = glnum_wk(j) + (i-j)
            glnum(j)    = glnum(j) + (i-j)
          endif
        enddo  

      enddo         ! ic

!!     debugging
!      call nekgsync()
!      if (nid.eq.0) then
!        do i=1,nunq2
!          write(6,*), nid,glnum(i),glnum_wk(i),tb1(i),tb2(i)
!        enddo
!      endif
!      call nekgsync()
!
!      if (nid.eq.1) then
!        do i=1,nunq2
!          write(6,*), nid,glnum(i),glnum_wk(i),tb1(i),tb2(i)
!        enddo
!      endif
!      call nekgsync()
!      call sleep(1)

!     Cycle 3: Numbering of unowned points
!--------------------------------------------------       

      j = 0
      k = 0
      do i=1,nunq
        if (owner(i).eq.nid) then
          j = j+1
          tb1(j)=x1local_u(i)
          tb2(j)=x2local_u(i)
        else
          k = k+1
          ta1(k)=x1local_u(i)
          ta2(k)=x2local_u(i)
        endif
      enddo
      nunq2  = j                     ! Unique owned elements
      unown  = k                     ! Unique un-owned elements
      

      call rzero(wk1,ntot1*2)
      call rzero(wk2,ntot1*2)
      call ifill(glnum2,-1.0,lz1*lz1*lelv)

      call izero(ind3,ntot1)

      ifseg(1)  = .true.
      ninseg(1) = nunq2
      nseg      = 1
      tol       = 1.0e-14

      do ic=1,csteps

        il = 2**(ic-1)
        srcid = nid - il
        dstid = nid + il
        if (srcid.lt.0)  srcid = srcid + np
        if (dstid.ge.np) dstid = dstid - np

!       Put unowned coordinates into work array 
        call copy(wk1(1),ta1,unown)
        call copy(wk2(1),ta2,unown)
        call icopy(glnum_wk(1),glnum2,unown)


!       set buffer for the number of elements to receive
        len = isize
        msg_id1 = irecv(0,cnsort,len)

!       send local size of array
        call csend(0,nunq2,len,dstid,0)

!       wait for length data
        call msgwait(msg_id1)

!       receive x1 data
        len = wdsize*cnsort
        msg_id2 = irecv(1,wk1(unown+1),len)

!       receive x2 data
        len = wdsize*cnsort
        msg_id3 = irecv(2,wk2(unown+1),len)

!       receive global numbering data
        len = isize*cnsort
        msg_id4 = irecv(3,glnum_wk(unown+1),len)

        
!       Send the owned data
        len = wdsize*nunq2
        call csend(1,tb1,len,dstid,0)
        call csend(2,tb2,len,dstid,0)
        len = isize*nunq2
        call csend(3,glnum,len,dstid,0)

!       finish communication            
        call msgwait(msg_id2)
        call msgwait(msg_id3)
        call msgwait(msg_id4)

        ifseg(1)  = .true.
        nsort     = unown + cnsort
        ninseg(1) = nsort
        nseg      = 1
        tol       = 1.0e-14

        ind3(1)   = 1
        do i=2,nsort
          ifseg(i) = .false.
          ind3(i)  = i              ! original array possitions
        enddo  

!       Perform local sorting
        do ip=1,2       ! pass
          do j=1,jloop  ! dimensions of sorting data
            nkey = 1
            i    = 1
            do iseg = 1,nseg
              if (j.eq.1) then
                call tuple_sort(wk1(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call swap_ip(wk2(i),ind,ninseg(iseg))
                call iswap_ip(ind3(i),ind,ninseg(iseg))
                call iswap_ip(glnum_wk(i),ind,ninseg(iseg))
              elseif (j.eq.2) then
                call tuple_sort(wk2(i),1,ninseg(iseg),1,nkey,ind,ta3)
                call swap_ip(wk1(i),ind,ninseg(iseg))
                call iswap_ip(ind3(i),ind,ninseg(iseg))
                call iswap_ip(glnum_wk(i),ind,ninseg(iseg))
              endif
              i = i + ninseg(iseg) 
            enddo

            do i=2,nsort
              if (j.eq.1) then
                if ((abs(wk1(i)-wk1(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif
              elseif (j.eq.2) then
                if ((abs(wk2(i)-wk2(i-1)).gt.tol)) then
                  ifseg(i) = .true.
                endif 
              endif
            enddo

!           Count up number of different segments
            nseg = 0
            do i=1,nsort
              if (ifseg(i)) then
                 nseg = nseg+1
                 ninseg(nseg) = 1
              else
                 ninseg(nseg) = ninseg(nseg) + 1
              endif
            enddo
          enddo     ! j
        enddo       ! ip


!       Reassign Global Numbering
        i = 1
        do iseg=1,nseg
          glnos = ivlmax(glnum_wk(i),ninseg(iseg))
          call ifill(glnum_wk(i),glnos,ninseg(iseg))
          i = i + ninseg(iseg)
        enddo 

!       Reassign Numbering of unown points            
        do i=1,nsort
          j = ind3(i)
          if (j.le.unown) then       ! local point
            glnum2(j) = glnum_wk(i)
          endif
        enddo  

      enddo         ! ic

!!     debugging
!      call nekgsync()
!      if (nid.eq.0) then
!        do i=1,unown
!          write(6,*), nid,glnum2(i),ta1(i),ta2(i)
!        enddo
!      endif
!      call nekgsync()
!
!      if (nid.eq.1) then
!        do i=1,unown
!          write(6,*), nid,glnum2(i),ta1(i),ta2(i)
!        enddo
!      endif
!      call nekgsync()
!      call sleep(1)


!     Local sort again with numbered points.
!-------------------------------------------------- 

      call copy(ta1,ym1,ntot1)
      call copy(ta2,zm1,ntot1)
      call ifill(glnum_wk,-1,ntot1)
      j = 0
      k = 0
      do i=1,nunq
        ta1(ntot1+i)=x1local_u(i)
        ta2(ntot1+i)=x2local_u(i)
        if (owner(i).eq.nid) then
          j = j+1
          glnum_wk(ntot1+i) = glnum(j)
        else
          k = k+1
          glnum_wk(ntot1+i) = glnum2(k)
        endif
      enddo

      nsort = ntot1+nunq
      do i=1,nsort
        ind(i) = i
      enddo  

      do i=1,nsort
        ifseg(i) = .false.
      enddo  

      ifseg(1)  = .true.
      ninseg(1) = nsort
      nseg      = 1
      tol       = 1.0e-14

!     Local sorting
      do ipass=1,2
        do j=1,jloop              ! For 2D we only do one dimension
          nkey = 1
          i    = 1
          do iseg = 1,nseg
            if (j.eq.1) then
              call tuple_sort(ta1(i),1,ninseg(iseg),1,nkey,ind2,ta3)
              call swap_ip(ta2(i),ind2,ninseg(iseg))
              call iswap_ip(ind(i),ind2,ninseg(iseg))
              call iswap_ip(glnum_wk(i),ind2,ninseg(iseg))
            elseif (j.eq.2) then
              call tuple_sort(ta2(i),1,ninseg(iseg),1,nkey,ind2,ta3)
              call swap_ip(ta1(i),ind2,ninseg(iseg))
              call iswap_ip(ind(i),ind2,ninseg(iseg))
              call iswap_ip(glnum_wk(i),ind2,ninseg(iseg))
            endif
            i = i + ninseg(iseg) 
          enddo

          do i=2,nsort
            if (j.eq.1) then
              if ((abs(ta1(i)-ta1(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif  
            elseif (j.eq.2) then
              if ((abs(ta2(i)-ta2(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif    
            endif
          enddo

!         Count up number of different segments
          nseg = 0
          do i=1,nsort
             if (ifseg(i)) then
                nseg = nseg+1
                ninseg(nseg) = 1
             else
                ninseg(nseg) = ninseg(nseg) + 1
             endif
             local_num(i) = nseg
          enddo
        enddo     ! j
      enddo       ! ipass

!     Reassign global numbering
      i = 1
      do iseg=1,nseg
        glnos = ivlmax(glnum_wk(i),ninseg(iseg))
        call ifill(glnum_wk(i),glnos,ninseg(iseg))
        i = i + ninseg(iseg)
      enddo 

      call i8zero(fs_gl_num,lx1*lz1*lelv)
!     Reassign global numbering 
      do i=1,nsort
        j = ind(i)
        if (j.le.ntot1) then       ! local point
          fs_gl_num(j) = glnum_wk(i)
!          ta3(j) = glnum_wk(i) + 0.0
        endif
      enddo  

!!     Gather-Scatter Setup      
      call fs_setupds(fs_gs_handle,ntot1,fs_gl_num)


      return
      end subroutine fs_map
!---------------------------------------------------------------------- 

      subroutine fs_gen_damping

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'GEOM'

      include 'FS_ALE'

      real dummy
      common /scrcg/ dummy(lx1,ly1,lz1,lelt)

      real x,x0,mu
      integer i,n

      call mntr_log(fs_id,fs_log,'Generating Damping Function')


      n = lx1*ly1*lz1*nelv

!     Find interface position      
      call copy(dummy,xm1,n)
      call col2(dummy,fs_mask,n)
      call fgslib_gs_op(fs_gs_handle,dummy,1,1,0)     ! 1 ==> +
      call col2(dummy,fs_vmult,n)


!     Create damping function      
      mu = 0.20
      do i=1,n
        x0               = dummy(i,1,1,1)
        x                = xm1(i,1,1,1)
        if (abs(x-x0).lt.fs_ofst) then
          fs_damp(i,1,1,1) = 1.0
        else  
          fs_damp(i,1,1,1) = exp(-((x-x0-fs_ofst)/mu)**2)
        endif  
      enddo

      return
      end subroutine fs_gen_damping        

!---------------------------------------------------------------------- 
      subroutine fs_setupds(gs_handle,ntot,glo_num)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
!      include 'NONCON'

      integer gs_handle
      integer*8 glo_num(1)

      integer mid,mp,nekcomm,nekgroup,nekreal      
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ntot
      real t0,t1
      real dnekclock          ! function

      t0 = dnekclock()

c     Initialize gather-scatter code
!      ntot      = lx1*ly1*lz1*nelv
      call fgslib_gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

c     call gs_chkr(glo_num)

      t1 = dnekclock() - t0
      if (nio.eq.0) then
         write(6,1) t1,gs_handle
    1    format('FS setupds time',1pe11.4,' seconds ',i3)
      endif
c
      return
      end
!-----------------------------------------------------------------------

      subroutine fs_mvmesh()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'MVGEOM'
      include 'INPUT'
      include 'FS_ALE'

      integer i,ntot1

      if (.not.ifmvbd) return
      if (.not.ifusermv) return

      param(59) = 1.       ! ifdef = .true.

      ntot1 = lx1*ly1*lz1*nelv

!     Zero out everything except the free surface      
      call opcopy(wx,wy,wz,vx,vy,vz)

      if (fs_ifgrid) then
!       Directly get smooth normals/tangentials by global
!       approximation of the interface.      
        call fs_smoothmv(wx,wy,wz)
      else
!       Move mesh in the normal direction.
!       With possible tangential correction.
!       And projection to smooth field
        call fs_smooth_meshmv(wx,wy,wz)
      endif  

!     Project interface velocity to the interior of the domain.
      call fs_int_project(wx,wy,wz)

!     Spatially damp out the mesh velocity      
      call col2(wx,fs_damp,ntot1)
      call col2(wy,fs_damp,ntot1)
      if (ndim.eq.3) then
        call col2(wz,fs_damp,ntot1)
      else  
!       We don't move mesh along Z
        call rzero(wz,ntot1)
      endif


!     Use Gordon Hall for mesh motion
      if (fs_ifgh) call fs_gh_correction(wx,wy,wz) 

      return
      end subroutine fs_mvmesh        

!-----------------------------------------------------------------------
!
!      subroutine fs_mvmeshn(ux,uy,uz)
!!     Only in 2D for now
!!     Project velocities on to Normal directions
!
!      implicit none
!
!      include 'SIZE'
!      include 'GEOM'
!      include 'INPUT'
!      include 'MASS'
!
!      include 'FS_ALE'
!
!      real ux(lx1,ly1,lz1,lelv)
!      real uy(lx1,ly1,lz1,lelv)
!      real uz(lx1,ly1,lz1,lelv)
!
!      integer i,n,nface
!      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
!      integer j1,j2,j3,nxyz
!      integer ifld
!
!      real rnor,rtn1
!
!      character cb*3
!
!      real dummy1,dummy2,dummy3
!      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
!     $               dummy2(lx1,ly1,lz1,lelt),
!     $               dummy3(lx1,ly1,lz1,lelt)
!
!      integer nsave
!      integer ies(2),ixs(2),iys(2)
!      save ies,ixs,iys,nsave
!      real wallvx(lx1*lelt),wallvy(lx1*lelt)
!      real tol
!      integer icalld
!      save icalld
!      data icalld /0/
!
!      ifld  = 1
!      nxyz  = lx1*ly1*lz1
!      nface = 2*ndim
!
!      do i = 1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wallvx(i) = ux(j1,j2,j3,e)
!        wallvy(i) = 0.0
!      enddo  
!        
!
!      do 200 e=1,nelv
!      do 200 ifc=1,nface
!        cb  = fs_cbc(ifc,e)
!        if (cb.eq.'INT') then
!          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
!          i = 0
!          do 220 j2=js2,jf2,jskip2
!          do 220 j1=js1,jf1,jskip1
!             i = i + 1
!!            normal component         
!             rnor = ( ux(j1,j2,1,e)*unx(i,1,ifc,e) +
!     $                uy(j1,j2,1,e)*uny(i,1,ifc,e) )
!!            remove tangential component
!             ux(j1,j2,1,e) = rnor*unx(i,1,ifc,e)
!             uy(j1,j2,1,e) = rnor*uny(i,1,ifc,e)
!
!
!  220      continue
!        endif                 
!  200 continue
!
!      call dsavg(ux)
!      call dsavg(uy)
!
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!
!        ux(j1,j2,j3,e) = wallvx(i)
!        uy(j1,j2,j3,e) = wallvy(i)
!      enddo   
!
!
!      return
!      end subroutine fs_mvmeshn        
!---------------------------------------------------------------------- 

      subroutine fs_projt(ux,uy,uz,tdir)
!     Project velocities on to tangential directions

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'

      include 'FS_ALE'

      real ux(lx1,ly1,lz1,lelv)
      real uy(lx1,ly1,lz1,lelv)
      real uz(lx1,ly1,lz1,lelv)

      integer i,n,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3,nxyz
      integer ifld

      real rnor,rtn1,rtn2
      integer tdir

      character cb*3

      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)

      integer nsave
      integer ies(lx1),ixs(lx1),iys(lx1)
      save ies,ixs,iys,nsave
      real wallvx(lx1),wallvy(lx1)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      nxyz  = lx1*ly1*lz1
      nface = 2*ndim

      do i = 1,nsave
        wallvx(i) = 0.0
        wallvy(i) = 0.0
      enddo  
        

      do 200 e=1,nelv
        do 200 ifc=1,nface
          cb  = fs_cbc(ifc,e)
          if (cb.eq.'INT') then
            call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
            i = 0
            do 220 j2=js2,jf2,jskip2
            do 220 j1=js1,jf1,jskip1
               i = i + 1
!              normal component         
               rnor = ( ux(j1,j2,1,e)*unx(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*uny(i,1,ifc,e) )
!              tangential compnent            
               rtn1 = ( ux(j1,j2,1,e)*t1x(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*t1y(i,1,ifc,e) )
               rtn2 = ( ux(j1,j2,1,e)*t2x(i,1,ifc,e) +
     $                  uy(j1,j2,1,e)*t2y(i,1,ifc,e) )
!              project to tangential direction
               if (tdir.eq.1) then
                 ux(j1,j2,1,e) = rtn1*t1x(i,1,ifc,e)
                 uy(j1,j2,1,e) = rtn1*t1y(i,1,ifc,e)
                 if (ndim.eq.3) uz(j1,j2,1,e) = rtn1*t1z(i,1,ifc,e)
               elseif (tdir.eq.2) then
                 ux(j1,j2,1,e) = rtn2*t2x(i,1,ifc,e)
                 uy(j1,j2,1,e) = rtn2*t2y(i,1,ifc,e)
                 if (ndim.eq.3) uz(j1,j2,1,e) = rtn2*t2z(i,1,ifc,e)
               else
                 if (nio.eq.0) write(6,*) 'unknown tdir:', tdir
                 call exitt  
               endif
  220        continue
          endif                 
  200   continue

      call dsavg(ux)
      call dsavg(uy)

      do i=1,nsave
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        ux(j1,j2,j3,e) = wallvx(i)
        uy(j1,j2,j3,e) = wallvy(i)
      enddo   

      return
      end subroutine fs_projt 
!---------------------------------------------------------------------- 

      subroutine fs_int_project(wx,wy,wz) 

!     project surface values to the interior of the domain        

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv)
      real wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

      call col2(wx,fs_mask,ntot1)
      call col2(wy,fs_mask,ntot1)
      if (ndim.eq.3) call col2(wz,fs_mask,ntot1)

!     Extend the velocity to the interior of the domain
!     using the custom gather-scatter operator
      call fgslib_gs_op(fs_gs_handle,wx,1,1,0)  ! 1 ==> +
      call fgslib_gs_op(fs_gs_handle,wy,1,1,0)  ! 1 ==> +
      if (ndim.eq.3) call fgslib_gs_op(fs_gs_handle,wz,1,1,0)  ! 1 ==> +

!     Take care of multiplicity 
      call col2(wx,fs_vmult,ntot1)
      call col2(wy,fs_vmult,ntot1)
      if (ndim.eq.3) call col2(wz,fs_vmult,ntot1)

      return
      end subroutine fs_int_project
!----------------------------------------------------------------------

      subroutine fs_symo_save

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'FS_ALE'

      integer nface,nsave
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3
      integer ifld

      character cb*3
      real tol
      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)


!     Store corner point locations for SYM/O conditions
      ifld  = 1 
      tol   = 1.0e-14
      nsave = 0
      nface = 2*ndim
      fs_nsymo = 0
      call rzero(dummy1,lx1*ly1*lz1*nelv)
      do e=1,nelv
      do ifc=1,nface
        cb  = cbc(ifc,e,ifld)
!       Change these conditions for actual case            
        if (cb.eq.'SYM'.or.cb.eq.'O  ') then
          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
          do j2=js2,jf2,jskip2
          do j1=js1,jf1,jskip1
            dummy1(j1,j2,1,e) = dummy1(j1,j2,1,e)+1.0
            if (abs(dummy1(j1,j2,1,e)-2.0).lt.tol) then
              nsave         = nsave + 1
              fs_ie(nsave)  = e
              fs_ix(nsave)  = j1
              fs_iy(nsave)  = j2
              fs_iz(nsave)  = 1
              fs_nsymo      = nsave
            endif
          enddo ! j2
          enddo ! j1
        endif                 
      enddo
      enddo

      return
      end subroutine fs_symo_save
!----------------------------------------------------------------------
      subroutine fs_gh_correction(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'INPUT'   ! cbc
      include 'TSTEP'

      include 'FS_ALE'

!     Mesh velocities      
      real wx(lx1,ly1,lz1,lelv)
      real wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,f,nfaces,ntot1,ifld
      real dti
      integer i,j,k

      character cb*3

!     mask for interior points
      real maski
      common /scrcg/ maski(lx1,ly1,lz1,lelt)

      integer nadd


      call mntr_log(fs_id,fs_log,'Gordon Hall Correction')

      ntot1  = lx1*ly1*lz1*nelv
      nfaces = 2*ndim
      ifld   = 1

!     create a mask for interior points.
      call rzero(maski,ntot1)

      do e=1,nelv     
      do f=1,nfaces   
         cb = cbc(f,e,ifld)
!        fill up edges with ones            
         call facev(maski,e,f,1.0,nx1,ny1,nz1)
!        Put ones on the interface 
!         if (cb.eq.'O  ') call facev(maski,e,f,1.0,nx1,ny1,nz1)
!        Don't move Periodic points         
         if (cb.eq.'P  ') call facev(maski,e,f,0.0,nx1,ny1,nz1)
!        Keeping wall fixed            
         if (cb.eq.'W  ') call facev(maski,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo

!     Only keep the edge values            
      call col2(wx,maski,ntot1)
      call col2(wy,maski,ntot1)
      if (ndim.eq.3) call col2(wz,maski,ntot1)

!     At first order, the displacement is just wx*dt, etc
!     We can do a better approximation if this works.
!     (Although better approximation is probably not needed)      
      call opcmult(wx,wy,wz,dt)           ! This is now displacement
      
!     Project displacement into the interior of the domain.
      nadd = 1    ! 1: corners only, 2: Corners & Edges, 
                  ! 3: Corners, Edges & Faces
      call fs_gordonhall(wx,wy,wz,nadd)

      dti = 1.0/dt
      call opcmult(wx,wy,wz,dti)           ! This is now velocities

      return
      end subroutine fs_gh_correction        
!----------------------------------------------------------------------       
      subroutine fs_gordonhall(mvx,mvy,mvz,nadd) 

!     fix up geometry irregularities
!     routine has been modified from fix_geom
!     to move internal points of the domain.

      implicit none

      include 'SIZE'
      include 'INPUT'   ! cbc
      include 'GEOM'
      include 'TSTEP'   ! ifield
      include 'WZ'

      integer lt
      parameter (lt = lx1*ly1*lz1)
      real xb,yb,zb,tmsk,tmlt,w1,w2
      common /scrns/ xb(lt,lelt),yb(lt,lelt),zb(lt,lelt),
     $               tmsk(lt,lelt),tmlt(lt,lelt),w1(lt),w2(lt)

      real xtmp,ytmp,ztmp
      common /scrsf/ xtmp(lx1,ly1,lz1,lelt),
     $               ytmp(lx1,ly1,lz1,lelt),
     $               ztmp(lx1,ly1,lz1,lelt)


      real mvx(lt,lelt),mvy(lt,lelt),mvz(lt,lelt) ! displaced boundary points
      
      integer i,e,kpass,n,f,nxyz,nfaces
      character*3 cb

      real s,xm,ym,zm,xx,yy,zz
      real glamax
      real glmax

      integer nadd      ! no of passes for which we add displacements.

      n      = nx1*ny1*nz1*nelt
      nxyz   = nx1*ny1*nz1
      nfaces = 2*ndim
      ifield = 1                   ! velocity field
!      if (ifheat) ifield = 2       ! temperature field

!     Use the temporary arrays for x,y,z      
      call opcopy(xtmp,ytmp,ztmp,xm1,ym1,zm1) 

      call rone  (tmlt,n)
      call dssum (tmlt,nx1,ny1,nz1)  ! denominator

      call rone  (tmsk,n)
      do e=1,nelfld(ifield)      ! fill mask where bc is periodic
      do f=1,nfaces              ! so we don't translate periodic bcs (z only)
         cb =cbc(f,e,ifield)
         if (cb.eq.'P  ') call facev (tmsk,e,f,0.0,nx1,ny1,nz1)
!        Keeping wall fixed            
         if (cb.eq.'W  ') call facev (tmsk,e,f,0.0,nx1,ny1,nz1)
      enddo
      enddo

      do kpass = 1,ndim+1     ! extra pass is just to test convergence

         call copy(xb,xtmp,n)
         call copy(yb,ytmp,n)
         if (ndim.eq.3) call copy(zb,ztmp,n)

         call dssum(xb,nx1,ny1,nz1)
         call dssum(yb,nx1,ny1,nz1)
         if (ndim.eq.3) call dssum(zb,nx1,ny1,nz1)

         xm = 0.
         ym = 0.
         zm = 0.

         do e=1,nelt
            do i=1,nxyz                       ! compute averages of geometry
               s     = 1./tmlt(i,e)
               xb(i,e) = s*xb(i,e)
               yb(i,e) = s*yb(i,e)
               if (ndim.eq.3) zb(i,e) = s*zb(i,e)

!              if there are tears between elemnts   
               xb(i,e) = xb(i,e) - xtmp(i,1,1,e)   ! local displacements
               yb(i,e) = yb(i,e) - ytmp(i,1,1,e)
               if (ndim.eq.3) zb(i,e) = zb(i,e) - ztmp(i,1,1,e)

               xb(i,e) = xb(i,e)*tmsk(i,e)
               yb(i,e) = yb(i,e)*tmsk(i,e)
               if (ndim.eq.3) zb(i,e) = zb(i,e)*tmsk(i,e)

!              add in the boundary displacements in the first pass   
               if (kpass.le.nadd) then
                 xb(i,e)=xb(i,e)+mvx(i,e)
                 yb(i,e)=yb(i,e)+mvy(i,e)
                 if (ndim.eq.3) zb(i,e)=zb(i,e)+mvz(i,e)
               endif   

               xm = max(xm,abs(xb(i,e)))
               ym = max(ym,abs(yb(i,e)))
               if (ndim.eq.3) zm = max(zm,abs(zb(i,e)))
            enddo

            if (kpass.le.ndim) then
               call gh_face_extend(xb(1,e),zgm1,nx1,kpass,w1,w2)
               call gh_face_extend(yb(1,e),zgm1,nx1,kpass,w1,w2)
               if (ndim.eq.3)
     $           call gh_face_extend(zb(1,e),zgm1,nx1,kpass,w1,w2)
            endif ! kpass.le.ndim

         enddo

         if (kpass.le.ndim) then
           call add2(xtmp,xb,n)
           call add2(ytmp,yb,n)
           if (ndim.eq.3) call add2(ztmp,zb,n)

           call sub2(mvx,xb,n)
           call sub2(mvy,yb,n)
           if (ndim.eq.3) call sub2(mvz,zb,n)
         endif
        
         xx = glamax(xb,n)
         yy = glamax(yb,n)
         if (ndim.eq.3) zz = glamax(zb,n)

         xm = glmax(xm,1)
         ym = glmax(ym,1)
         if (ndim.eq.3) zm = glmax(zm,1)

!         if (nio.eq.0) write(6,1) xm,ym,zm,xx,yy,zz,kpass
!    1    format(1p6e12.4,' fs GH corr',i2)

      enddo

      call opcopy(mvx,mvy,mvz,xtmp,ytmp,ztmp)
      call opsub2(mvx,mvy,mvz,xm1,ym1,zm1)
!      call opcopy(xm1,ym1,zm1,xtmp,ytmp,ztmp)
!      param(59) = 1.       ! ifdef = .true.
      
      return
      end subroutine fs_gordonhall
c-----------------------------------------------------------------------

      subroutine fs_mvmesh_linear()

!     Linear Approximation within elements for mesh motion.        

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'MVGEOM'
      include 'INPUT'
      include 'FS_ALE'

      integer i,ntot1

      if (.not.ifmvbd) return
      if (.not.ifusermv) return

      param(59) = 1.       ! ifdef = .true.

      ntot1 = lx1*ly1*lz1*nelv

!     Zero out everything except the free surface      
!      call opcopy(wx,wy,wz,vx,vy,vz)
      call col2(wx,fs_crmask,ntot1)
      call col2(wy,fs_crmask,ntot1)
      if (if3d) call col2(wz,fs_crmask,ntot1)

!     Project interface velocity to the interior of the domain.
      call fs_int_project(wx,wy,wz)

!     Spatially damp out the mesh velocity      
      call col2(wx,fs_damp,ntot1)
      call col2(wy,fs_damp,ntot1)
      if (ndim.eq.3) then
        call col2(wz,fs_damp,ntot1)
      else  
!       We don't move mesh along Z
        call rzero(wz,ntot1)
      endif

!     Use Gordon Hall for mesh motion
      if (fs_ifgh) call fs_gh_correction(wx,wy,wz)

!      call outpost(wx,wy,wz,pr,t,'msh') 

      return
      end subroutine fs_mvmesh_linear 

!-----------------------------------------------------------------------


      subroutine fs_get_intpos(pos)

!     Get interface position        


      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'FS_ALE'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real pos(lt)
      real pos2(lt)

      integer i,n

      real mu,x

      real wght,wghte,dmy1,dmy2
      common /scrmg/  wght(lt)            ! Local Weights
     $              , wghte(lt)           ! Sum of weights (extended)
     $              , dmy1(lt)
     $              , dmy2(lt)

      n = lx1*ly1*lz1*nelv

!     Get Weights
      mu = 0.01
      do i=1,n
        x             = abs(t(i,1,1,1,1))       ! absolute distance from interface 
        wght(i)       = exp(-(x/mu)**2)         ! local weights
        pos(i)        = xm1(i,1,1,1)            ! which coordinate?
      enddo

      call copy(wghte,wght,n)

!     Add weights along the mapping and extend to each point      
      call fgslib_gs_op(fs_gs_handle,wghte,1,1,0)  ! 1 ==> +
      call col2(wghte,fs_vmult,n)

      call col2(pos,wght,n)                    ! convolution
!     Extend convolution
      call fgslib_gs_op(fs_gs_handle,pos,1,1,0)  ! 1 ==> +
      call col2(pos,fs_vmult,n)

!     Get interface position
!     Position = convolution/weight
!     Each Point knows the interface postion at that (y,z)
      call invcol2(pos,wghte,n)

      return
      end subroutine fs_get_intpos 
!-----------------------------------------------------------------------

      function cornerpoint(ix,iy,iz)

      implicit none

      include 'SIZE'

      logical cornerpoint
      
      integer ix,iy,iz

      cornerpoint = .false.

!      if (ndim.eq.3) then
      
      if (ix.eq.1  .and.iy.eq.1  .and.iz.eq.1  ) cornerpoint = .true.
      if (ix.eq.1  .and.iy.eq.1  .and.iz.eq.lz1) cornerpoint = .true.
      if (ix.eq.1  .and.iy.eq.ly1.and.iz.eq.1  ) cornerpoint = .true.
      if (ix.eq.1  .and.iy.eq.ly1.and.iz.eq.lz1) cornerpoint = .true.
      if (ix.eq.lx1.and.iy.eq.1  .and.iz.eq.1  ) cornerpoint = .true.
      if (ix.eq.lx1.and.iy.eq.1  .and.iz.eq.lz1) cornerpoint = .true.
      if (ix.eq.lx1.and.iy.eq.ly1.and.iz.eq.1  ) cornerpoint = .true.
      if (ix.eq.lx1.and.iy.eq.ly1.and.iz.eq.lz1) cornerpoint = .true.

      return
      end
!---------------------------------------------------------------------- 

      subroutine fs_solve_springmass(pos,wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real dmy1,dmy2,dmy3,dmy4
      common /scrmg/  dmy1(lt)
     $              , dmy2(lt)
     $              , dmy3(lt)
     $              , dmy4(lt)
     
      real pos(lt)
      real wx(lt),wy(lt),wz(lt)

      real gridvx(lt),gridvx_lag(lt,2)
      real rhs_vx(lt),rhs_vxlag(lt,2)

      common /fs_gridsolv/ gridvx,gridvx_lag,rhs_vx,rhs_vxlag

      real gridx(lt),gridx_lag(lt,2)
      real rhs_x(lt),rhs_xlag(lt,2)

      common /fs_gridsolx/ gridx,gridx_lag,rhs_x,rhs_xlag


      real fs_dt,fs_dtlag(10),fs_bd(10),fs_extk(10)
      real dti

      common /fs_gridtime/fs_dt,fs_dtlag,fs_bd,fs_extk 

      integer i,n
      integer ilag,irst,nbdmsh

      real damp,stiff


      n = lx1*ly1*lz1*nelv

      stiff = uparam(8)
      damp  = uparam(9)

      if (istep.eq.1) then

        call cfill(gridx,1.23,n)

        call rzero(gridvx,n)
!        call cfill(gridvx,1.0,n)

        call rzero(gridx_lag,lt*2)
        call rzero(gridvx_lag,lt*2)

        call rzero(rhs_x,n)
        call rzero(rhs_xlag,lt*2)

        call rzero(rhs_vx,n)
        call rzero(rhs_vxlag,lt*2)

        call rzero(fs_dtlag,10)

      endif

      irst = param(46)

      do i=10,2,-1
        fs_dtlag(i) = fs_dtlag(i-1)
      enddo
      call setdt_simple
      fs_dtlag(1) = dt
      if (istep.eq.1 .and. irst.le.0) fs_dtlag(2) = dt

      call copy(dtlag,fs_dtlag,10)        ! only for now
      fs_dt    = dt

      time     = time + dt

!     BD coefficients      
      call setordbd
      if (irst.gt.0) nbd = nbdinp
      call rzero(fs_bd,10)
      call setbd(fs_bd,fs_dtlag,nbd)
     

!     Extk Coefficients      
      if (param(27).lt.0) then
         nab = nbdinp
      else
         nab = 3
      endif
      if (istep.lt.nab.and.irst.le.0) nab = istep
      call rzero   (fs_extk,10)
      call setabbd (fs_extk,fs_dtlag,nab,nbd)

!     only for now!!      
      if (ifmvbd) then
         nbdmsh = 1
         nabmsh = param(28)
         if (nabmsh.gt.istep .and. irst.le.0) nabmsh = istep
         call rzero   (abmsh,10)
         call setabbd (abmsh,fs_dtlag,nabmsh,nbdmsh)
      endif

!     Build RHS
!------------------------------       
      do i=1,n
        rhs_vx(i) = - damp*gridvx(i) - stiff*(gridx(i)-pos(i))
        rhs_x(i)  =  gridvx(i)
      enddo

!     Extrapolate rhs(vx)
!------------------------------       
      call copy(dmy1,rhs_vxlag(1,2),n)       
      call cmult(dmy1,fs_extk(3),n)
      call add2s2(dmy1,rhs_vxlag(1,1),fs_extk(2),n)
!     Add current time contribution      
      call add2s2(dmy1,rhs_vx,fs_extk(1),n)
!     Save lagged variables      
      call copy(rhs_vxlag(1,2),rhs_vxlag(1,1),n)
      call copy(rhs_vxlag(1,1),rhs_vx,n)
    

!     Extrapolate rhs(x)
!------------------------------       
      call copy(dmy2,rhs_xlag(1,2),n)       
      call cmult(dmy2,fs_extk(3),n)
      call add2s2(dmy2,rhs_xlag(1,1),fs_extk(2),n)
!     Add current time contribution      
      call add2s2(dmy2,rhs_x,fs_extk(1),n)
!     Save lagged variables
      call copy(rhs_xlag(1,2),rhs_xlag(1,1),n)
      call copy(rhs_xlag(1,1),rhs_x,n)

!     Backward difference terms
!------------------------------       
      call copy(rhs_vx,gridvx,n)
      call cmult(rhs_vx,fs_bd(2),n)

      call copy(rhs_x,gridx,n)
      call cmult(rhs_x,fs_bd(2),n)

      do ilag=2,nbd
      do i=1,n
        rhs_vx(i) = rhs_vx(i) + gridvx_lag(i,ilag-1)*fs_bd(ilag+1)
        rhs_x(i)  = rhs_x(i)  + gridx_lag(i,ilag-1)*fs_bd(ilag+1)
      enddo 
      enddo
      dti = 1.0/fs_dt
      call add2s2(dmy1,rhs_vx,dti,n)
      call add2s2(dmy2,rhs_x,dti,n)

!     Save lagged variables
      do ilag=2,2,-1
        call copy(gridvx_lag(1,ilag),gridvx_lag(1,ilag-1),n)
        call copy(gridx_lag(1,ilag),gridx_lag(1,ilag-1),n)
      enddo
      call copy(gridvx_lag,gridvx,n)
      call copy(gridx_lag,gridx,n)

!     Solution
      do i=1,n
        gridvx(i) = dmy1(i)*fs_dt/fs_bd(1)
        gridx(i)  = dmy2(i)*fs_dt/fs_bd(1)
      enddo  

      call copy(wx,gridvx,n)
      call copy(wy,gridx,n)

      return
      end subroutine 
!---------------------------------------------------------------------- 

      subroutine setdt_simple
c
c     Set the new time step. All cases covered.
c
      include 'SIZE'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'PARALLEL'

      common /scruz/ cx(lx1*ly1*lz1*lelt)
     $ ,             cy(lx1,ly1,lz1,lelt)
     $ ,             cz(lx1,ly1,lz1,lelt)

      common /cprint/ ifprint
      logical         ifprint
      common /udxmax/ umax
      REAL     DTOLD
      SAVE     DTOLD
      DATA     DTOLD /0.0/
      REAL     DTOpf
      SAVE     DTOpf
      DATA     DTOpf /0.0/
      logical iffxdt
      save    iffxdt
      data    iffxdt /.false./
C

      if (param(12).lt.0.or.iffxdt) then
         iffxdt    = .true.
         param(12) = abs(param(12))
         dt        = param(12)
         dtopf     = dt

      else if (param(84).ne.0.0) then
         if (dtold.eq.0.0) then
            dt   =param(84)
            dtold=param(84)
            dtopf=param(84)
            return
         else
            dtold=dt
            dtopf=dt
            dt=dtopf*param(85)
            dt=min(dt,param(12))
         endif
      endif

      return
      end subroutine 
!---------------------------------------------------------------------- 















