"""
    Dot(V::SymPy.Sym, U::SymPy.Sym)

    Dot product between V and U
    The placement of unit vectors is assumed to be handled by the user
"""
function Dot(V,U)

# Apparently there should be no space between """ and the function definition

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


  tmp    = V*U

# Unit vector rotations
#  tmp = du.xreplace(Dict([rp1 rp2 rp3 rp4]))

  du  = Sym.expand(tmp)
#
  tmp2 = du.subs(Dict([rp5 rp6 rp7 rp8 rp9 rp10 rp11 rp12 rp13]))

  dotvu = tmp2

  return dotvu
end 
