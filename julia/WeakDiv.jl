function WeakDiv(U,v,xin...)

# Divergence in cylindrical coordinates
#

  dir = "left"

  if isempty(xin)
    dir = "left"
  else
    dir = xin[1]
  end  

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
    du    = -Sym.diff(v,r)*er*U - 1/r*Sym.diff(v,θ)*eθ*U - Sym.diff(v,t)*et*U 
  else  
    du    = -U*er*diff(v,r) - U*eθ*Sym.diff(v,θ)*1/r  - U*et*Sym.diff(v,t)
  end  

# Unit vector gradients
  tmp = du.xreplace(Dict([rp1 rp2 rp3 rp4]))

  du  = Sym.expand(tmp)
#
  tmp2 = du.subs(Dict([rp5 rp6 rp7 rp8 rp9 rp10 rp11 rp12 rp13]))

  divu = tmp2

#  divu = du

  return divu
end 
