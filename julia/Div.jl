"""
    Div(U::SymPy.Sym,dir::String)

    Divergence of U in Cylindrical Coordinates
    dir defines on which side the unit vectors are placed.
    dir = left/right  
"""
function Div(U::SymPy.Sym,dir::String)

# Divergence in cylindrical coordinates

  rp1   = (Sym.diff(er,r),0)
  rp2   = (Sym.diff(er,θ),eθ)
  rp3   = (Sym.diff(eθ,r),0)
  rp4   = (Sym.diff(eθ,θ),-er)

  rp5   = (et*et,1.0)
  rp6   = (et*er,0.0)
  rp7   = (et*eθ,0.0)
  rp8   = (er*er,1.0)
  rp9   = (er*eθ,0.0)
  rp10  = (er*et,0.0)
  rp11  = (eθ*eθ,1.0)
  rp12  = (eθ*er,0.0)
  rp13  = (eθ*et,0.0)


  if lowercase(dir)=="left" 
    du    = er*Sym.diff(U,r) + 1/r*eθ*Sym.diff(U,θ) + et*Sym.diff(U,t)
  elseif lowercase(dir)=="right"  
    du    = Sym.diff(U,r)*er + 1/r*Sym.diff(U,θ)*eθ + Sym.diff(U,t)*et
  else
    println("Unknown Placement of unit vectors in Div.jl, dir=$dir")
  end  

# Unit vector rotations
  tmp = du.xreplace(Dict([rp1 rp2 rp3 rp4]))

  du  = Sym.expand(tmp)
#
  tmp2 = du.subs(Dict([rp5 rp6 rp7 rp8 rp9 rp10 rp11 rp12 rp13]))

  divu = tmp2

  return divu
end 
