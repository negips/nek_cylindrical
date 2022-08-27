function GradT(U,r,θ,t,er,eθ,et,xin...)
       
# Calculate the result as a tensor
# if the input is a vector
 
  if isempty(xin)
    dir = "left"
  else
    dir = xin[1]
  end

  tmp  = Sym.diff(U,t)
  dudt = Sym.simplify(ReplaceUnitVecDeriv(tmp))

 # Get er Coefficients
  rp1  = (er,1.0)
  rp2  = (eθ,0.0)
  rp3  = (et,0,0)

  tmp  = dudt.subs(Dict([rp1 rp2 rp3]))
  dudt_er = tmp

# Get eθ coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,1.0)
  rp3  = (et,0.0)

  tmp  = dudt.subs(Dict([rp1 rp2 rp3]))
  dudt_eθ = tmp

# Get et coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,0.0)
  rp3  = (et,1.0)

  tmp  = dudt.subs(Dict([rp1 rp2 rp3]))
  dudt_et = tmp

  dudt_vec = [dudt_et dudt_er dudt_eθ]

  if lowercase(dir)=="left"
    dudt     = et*dudt
  else
    dudt  = Sym.diff(U,θ)*et
  end 

# Gradients
  mdudt = Sym.simplify(dudt_vec).as_mutable()

  return dudt,mdudt
end



