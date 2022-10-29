# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 

#      include("JNek_PARALLEL.jl")         # JNek_PARALLEL
      
      include("/home/prabal/workstation/git/julia_scripts/nek_scripts/JNek_IO.jl")               # JNek_IO

#      import JNek_IO
      
      using MPI

      nid0  = 0

      if (MPI.Initialized() == false)
        MPI.Init()
      end  
        
      comm = MPI.COMM_WORLD
      rank = MPI.Comm_rank(comm)

      rea   = "box.rea"
      re2   = "three.re2"
      map   = "three.ma2"
      file1   = "prcthree0.f00001"
      file2   = "prcthree0.f00002"
      file3   = "prcthree0.f00003"
      file4   = "prcthree0.f00004"
      file5   = "prcthree0.f00005"
      file6   = "prcthree0.f00006"
      file7   = "prcthree0.f00007"
      file8   = "extthree0.f00001"
      file9   = "tmpthree0.f00001"
      file10  = "tmpthree0.f00002"
      file11  = "tmpthree0.f00003"
      file12  = "tmpthree0.f00004"
      file13  = "tmpthree0.f00005"
      file14  = "excthree0.f00001"
      file15  = "lapthree0.f00001"
      file16  = "lapthree0.f00002"
      file17  = "cylthree0.f00001"

#      wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)
      re2fld1 = JNek_IO.read_re2_struct(re2,nid0)

#      fld1 = JNek_IO.read_fld_struct(file1,MPI,nid0)
#
#      fld2 = JNek_IO.read_fld_struct(file2,MPI,nid0)
#
#      fld3 = JNek_IO.read_fld_struct(file3,MPI,nid0)

#      fld4 = JNek_IO.read_fld_struct(file4,MPI,nid0)

#      fld5 = JNek_IO.read_fld_struct(file5,MPI,nid0)
#
#      fld6 = JNek_IO.read_fld_struct(file6,MPI,nid0)
#
#      fld7 = JNek_IO.read_fld_struct(file7,MPI,nid0)
#
#      fld8 = JNek_IO.read_fld_struct(file8,MPI,nid0)

#      fld9 = JNek_IO.read_fld_struct(file9,MPI,nid0)
#
#      fld10 = JNek_IO.read_fld_struct(file10,MPI,nid0)
#
#      fld11 = JNek_IO.read_fld_struct(file11,MPI,nid0)
#
#      fld12 = JNek_IO.read_fld_struct(file12,MPI,nid0)

#      fld13 = JNek_IO.read_fld_struct(file13,MPI,nid0)

#      fld14 = JNek_IO.read_fld_struct(file14,MPI,nid0)

#      fld15 = JNek_IO.read_fld_struct(file15,MPI,nid0)

#      fld16 = JNek_IO.read_fld_struct(file16,MPI,nid0)

      fld17 = JNek_IO.read_fld_struct(file17,MPI,nid0)


#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x2,y2,z2,u2,v2,w2,p2,t2 = JNek_IO.read_fld(file12,MPI,nid0)

#
#       plot(yvec,pvec7)

      close("all")

      iy    = 1
      el    = 1:3
      x17   = fld17.x[:,iy,1,el]
      x     = x17[:]
      u17   = fld17.u[:,iy,1,el]
      v17   = fld17.v[:,iy,1,el]
      plot(x,u17[:])
#      plot(x,v17[:])

      if rank == 0
        println("Done")
      end  

#      MPI.Finalize()





