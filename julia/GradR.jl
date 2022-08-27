function GradR(U, r,θ,t,er,eθ,et,xin...)

# Calculate the result as a tensor
# if the input is a vector
 
  if isempty(xin)
    dir = "left"
  else
    dir = xin[1]
  end  

  tmp  = Sym.diff(U,r)
  dudr = Sym.simplify(ReplaceUnitVecDeriv(tmp))

# Complicated Manipulations to put things in a matrix/Vector Form
# Should be an easier way

# Get er Coefficients
  rp1  = (er,1.0)
  rp2  = (eθ,0.0)
  rp3  = (et,0.0)

  tmp  = dudr.subs(Dict([rp1 rp2 rp3]))
  dudr_er = tmp

# Get eθ coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,1.0)
  rp3  = (et,0.0)

  tmp  = dudr.subs(Dict([rp1 rp2 rp3]))
  dudr_eθ = tmp

# Get et coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,0.0)
  rp3  = (et,1.0)

  tmp  = dudr.subs(Dict([rp1 rp2 rp3]))
  dudr_et = tmp

#  dudr_vec = Sym.MatrixSymbol("dUdR", 1,3)
#  dudr_vec[1,1] = dudr_er
#  dudr_vec[1,2] = dudr_eθ
#  dudr_vec[1,3] = dudr_et

  dudr_vec = [dudr_et dudr_er dudr_eθ]

  if lowercase(dir)=="left"
    dudr     = er*dudr
  else
    dudr     = dudr*er
  end 

# Put it in a Matrix

  mdudr = Sym.simplify(dudr_vec).as_mutable()
  return dudr,mdudr
end 
