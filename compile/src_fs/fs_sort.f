!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for free surface mesh movement.
!
!----------------------------------------------------------------------
!====================================================================== 
      subroutine fs_map1

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'FS_ALE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

!     Temporary arrays
      real ta1,ta2,ta3
      integer ninseg
      logical ifseg
      common /scrmg/ ta1    (lt)
     $ ,             ta2    (lt)
     $ ,             ta3    (lt)
     $ ,             ninseg (lt)
     $ ,             ifseg  (lt)

!     Temporary arrays
      integer ind
      integer wki
      common /scrcg/ ind (lt)
     $ ,             wki (lt)

      real tb1,tb2
      integer ninseg1,ninseg2
      logical ifseg1,ifseg2
      common /scruz/ tb1     (lt)
     $ ,             tb2     (lt)
     $ ,             ninseg1 (lt)
     $ ,             ninseg2 (lt)
     $ ,             ifseg1  (lt)
     $ ,             ifseg2  (lt)

      real radius(lt)

      integer i,j,k,n,ntot1
      integer nseg
      
      integer local_num_u(lx1*lz1*lelt)   ! unique local numbering
      integer global_num_u(lx1*lz1*lelt)  ! unique Global numbering
      integer local_num(lx1*ly1*lz1*lelt) ! Local numbering
      integer global_num(lx1*ly1*lz1*lelt)! Global numbering
      integer nlocalu                     ! No. of unique local elements

      real rad,theta
      real y,z
      real pi

      ntot1 = lx1*ly1*lz1*nelv

      call copy(ta1,ym1,ntot1)
      call copy(ta2,zm1,ntot1)

!     Perform Local sort
!     Assuming data is continuous
!     Inputs:
!     ta1    - Primary field entries
!     ta2    - Secondary field entries
!     ntot1  - no of entries
!     ndim   - No of fields (2/3)
!     Outputs:
!     ta1    - sorted entries (Primary field)
!     ta2    - sorted entries (Secondary field)
!     ind    - permutation indicies
!     nseg   - no of segments
!     ninseg - entries in each segment
!     ifseg  - segment breaks
!     Work:
!     ta3    - real work array
!     wki    - integer work array
      call fs_local_tuple_sort(ta1,ta2,ind,ifseg,nseg,
     $                         ninseg,ntot1,ndim,ta3,wki)

!     Contract array to get locally unique points
      j = 0
      do i=1,ntot1
        if (ifseg(i)) then
          j = j + 1  
          tb1(j)              = ta1(i)
          tb2(j)              = ta2(i)
          local_num_u(j)      = j
        endif
        k = ind(i)            ! i - Position in current array
                              ! k - Position in original array 
        local_num(k) = j      ! j - Local numbering
      enddo
      nlocalu = j             ! Should be the same as nseg

      call copy(ta1,tb1,nlocalu)
      call copy(ta2,tb2,nlocalu)
      call fs_global_sort(tb1,tb2,global_num_u,nlocalu)

      do i=1,ntot1
        j = local_num(i)
        k = global_num_u(j)
        fs_gl_num(i) = k
      enddo

!     Gather-Scatter Setup      
      call fs_setupds(fs_gs_handle,ntot1,fs_gl_num)


      return
      end subroutine fs_map1
!---------------------------------------------------------------------- 

      subroutine fs_global_sort(fld1,fld2,gnum,n)

!     Get global numbering of entries

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      integer n               ! No. of entries.
      real fld1(n)            ! Input field 1
      real fld2(n)            ! Input field 2
      integer gnum(n)         ! Permutation indicies

      integer alen            ! Work array length
      integer iglsum

      integer gl_n            ! Total number of global entries to sort
      character*64 str


      alen = 2*lt

      gl_n = iglsum(n,1)

      if (gl_n.lt.alen) then
!       Everything fits within one processor memory
        call blank(str,64)
        write(str,'(A25,I5,A11,I6)') 'Single Processor Sort: N=',
     $      gl_n,', Arr. Len=', alen
        call mntr_log(fs_id,fs_log,str)
        
        call fs_global_sort_sp(fld1,fld2,gnum,n)
        return
      else
        call blank(str,64)
        write(str,'(A27,I5,A11,I6)') 'Multiple Processor Sort: N=',
     $      gl_n,', Arr. Len=', alen
        call mntr_log(fs_id,fs_log,str)
        call mntr_log(fs_id,fs_log,'Not Yet implemented. Exitting')
        call exitt
      endif


      return
      end subroutine fs_global_sort
!---------------------------------------------------------------------- 
      subroutine fs_global_sort_sp(fld1,fld2,gnum,n)

!     Get global numbering of entries
!     When the Single processor work memory is large enough
!     for all entries.

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      integer n               ! No. of entries.
      real fld1(n)            ! Input field 1
      real fld2(n)            ! Input field 2
      integer gnum(n)         ! Permutation indicies

!     Temporary arrays
      real wk1,wk2,wk3
      common /scrns/ wk1(lt*2),
     $               wk2(lt*2),
     $               wk3(lt*3)

!     Temporary arrays
      integer ind,wki
      common /scrvh/ ind (lt*2)
     $ ,             wki (lt*2)

!     Temporary arrays
      integer ninseg
      logical ifseg
      common /screv/ ninseg (lt*2)
     $ ,             ifseg  (lt*2)

      integer nseg
      integer nglobalu

      integer alen            ! Work array length
      integer iglsum
      integer igl_running_sum

      integer gl_n            ! Total number of global entries to sort
      integer gl_n0
      character*64 str

      integer i,j,k

      alen = 2*lt             ! Available array length

      gl_n  = iglsum(n,1)
      gl_n0 = igl_running_sum(n)
      gl_n0 = gl_n0 - n       ! entries before this nid.

      call rzero(wk1,alen)
      call rzero(wk2,alen)
      call rzero(wk3,3*lt)

      do i=1,n
        wk1(gl_n0+i) = fld1(i)
        wk2(gl_n0+i) = fld2(i)
      enddo

!     Collect all entries in each processor
      call gop(wk1,wk3,'+  ',gl_n)
      call gop(wk2,wk3,'+  ',gl_n)

!     Local sorting of the large array
!     Inputs:
!     wk1    - Primary field entries
!     wk2    - Secondary field entries
!     gl_n   - no of entries
!     ndim   - No of fields (2/3)
!     Outputs:
!     wk1    - sorted entries (Primary field)
!     wk2    - sorted entries (Secondary field)
!     ind    - permutation indicies
!     nseg   - no of segments
!     ninseg - entries in each segment
!     ifseg  - segment breaks
!     Work:
!     wk3    - real work array
!     wki    - integer work array
      call fs_local_tuple_sort(wk1,wk2,ind,ifseg,nseg,
     $                         ninseg,gl_n,ndim,wk3,wki)


!     Contract array to get Globally unique points
!     using wki as global numbering array (temporary)
      call izero(wki,gl_n)
      j = 0
      do i=1,gl_n
        if (ifseg(i)) then
          j = j + 1  
        endif
        k = ind(i)            ! i - Position in current array
                              ! k - Position in original array 
        wki(k) = j            ! j - Global numbering
      enddo
      nglobalu = j             ! Should be the same as nseg

      call icopy(gnum,wki(gl_n0+1),n)

      return
      end subroutine fs_global_sort_sp
!---------------------------------------------------------------------- 

      subroutine fs_local_tuple_sort(fld1,fld2,ind,ifseg,nseg,
     $                               ninseg,n,fdim,wkr,wki)

!     Perform Local sort
!     Assuming data is continuous
!     Inputs:
!     fld1   - Primary field entries
!     fld2   - Secondary field entries
!     n      - no of entries
!     fdim   - No of fields (2/3)
!     Outputs:
!     fld1   - sorted entries (Primary field)
!     fld2   - sorted entries (Secondary field)
!     ind    - permutation indicies
!     ifseg  - segment breaks
!     nseg   - no of segments
!     ninseg - entries in each segment
!     Work:
!     wkr    - real work array
!     wki    - integer work array

      implicit none   

      include 'SIZE'
      
      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      integer n               ! No. of entries to sort
      real fld1(n)            ! Input field 1
      real fld2(n)            ! Input field 2
      integer ind(n)          ! Permutation indicies
      integer ninseg(n)       ! No. of entries in each segment
      logical ifseg(n)        ! Segment breaks/Start of new segment
      integer nseg            ! No. of segments

!     Work arrays
      real wkr(n)             ! Work (reals)
      integer wki(n)          ! Work (integer)

      integer fdim            ! No. of fields

      integer ipass,iseg,jloop
      integer i,j
      integer lda

      integer nkey
      integer key(1)
      real tol

      do i=1,n
        ind(i)    = i      ! Initialize keys
        wki(i)    = 0
      enddo 

!     Segments (initialization)        
      do i=1,n
        ifseg(i)  = .false.
      enddo  

      ifseg(1)    = .true.
      ninseg(1)   = n
      nseg        = 1
      tol         = 1.0e-14

      jloop = fdim-1
      lda   = 1         ! leading dimension of array
!     Local sorting
      do ipass=1,2
        do j=1,jloop
          key(1) = 1  
          nkey   = 1
          i      = 1
          do iseg = 1,nseg
            if (j.eq.1) then
              call tuple_sort(fld1(i),lda,ninseg(iseg),key,nkey,
     $                        wki,wkr)
              call iswap_ip(ind(i),wki,ninseg(iseg))
              call swap_ip(fld2(i),wki,ninseg(iseg))
            elseif (j.eq.2) then
              call tuple_sort(fld2(i),lda,ninseg(iseg),key,nkey,
     $                        wki,wkr)
              call iswap_ip(ind(i),wki,ninseg(iseg))
              call swap_ip(fld1(i),wki,ninseg(iseg))
            endif
            i = i + ninseg(iseg) 
          enddo

          do i=2,n
            if (j.eq.1) then
              if ((abs(fld1(i)-fld1(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif  
            elseif (j.eq.2) then
              if ((abs(fld2(i)-fld2(i-1)).gt.tol)) then
                ifseg(i) = .true.
              endif    
            endif
          enddo

!         Count up number of different segments
          nseg = 0
          do i=1,n
            if (ifseg(i)) then
              nseg = nseg + 1
              ninseg(nseg) = 1
            else
              ninseg(nseg) = ninseg(nseg) + 1
            endif
          enddo
        enddo     ! j
      enddo       ! ipass


      return
      end subroutine fs_local_tuple_sort
!---------------------------------------------------------------------- 

      subroutine msgwaitstat(imsg,status)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer ierr,imsg
c
      call mpi_wait (imsg,status,ierr)
c
      return
      end

!---------------------------------------------------------------------- 

      subroutine lcopy(a,b,n)

      logical a(1),b(1)
      integer i,n

      do i=1,n
        a(i)=b(i)
      enddo

      return
      end
!----------------------------------------------------------------------
 
!      integer stat(mpi_status_size)




