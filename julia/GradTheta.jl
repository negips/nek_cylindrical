function GradTheta(U,r,θ,t,er,eθ,et,xin...)
       
# Calculate the result as a tensor
# if the input is a vector
 
  if isempty(xin)
    dir = "left"
  else
    dir = xin[1]
  end

  tmp  = Sym.diff(U,θ)
  dudθ = Sym.simplify(ReplaceUnitVecDeriv(tmp))

 # Get er Coefficients
  rp1  = (er,1.0)
  rp2  = (eθ,0.0)
  rp3  = (et,0.0)

  tmp  = dudθ.subs(Dict([rp1 rp2 rp3]))
  dudθ_er = tmp

# Get eθ coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,1.0)
  rp3  = (et,0.0)

  tmp  = dudθ.subs(Dict([rp1 rp2 rp3]))
  dudθ_eθ = tmp

# Get et coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,0.0)
  rp3  = (et,1.0)

  tmp  = dudθ.subs(Dict([rp1 rp2 rp3]))
  dudθ_et = tmp

  dudθ_vec = (1/r).*[dudθ_et dudθ_er dudθ_eθ]

  if lowercase(dir)=="left"
    dudθ     = eθ*(1/r)*dudθ
  else
    dudθ  = Sym.diff(U,θ)*(1/r)*eθ
  end 

# Gradients
  mdudθ = Sym.simplify(dudθ_vec).as_mutable()

  return dudθ,mdudθ
end



