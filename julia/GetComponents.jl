function GetComponents(U)

# Complicated Manipulations to put things in a matrix/Vector Form
# Should be an easier way

# Get er Coefficients
  rp1  = (er,1.0)
  rp2  = (eθ,0.0)
  rp3  = (et,0.0)

  tmp  = U.subs(Dict([rp1 rp2 rp3]))
  U_er = tmp

# Get eθ coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,1.0)
  rp3  = (et,0.0)

  tmp  = U.subs(Dict([rp1 rp2 rp3]))
  U_eθ = tmp

# Get et coefficients  
  rp1  = (er,0.0)
  rp2  = (eθ,0.0)
  rp3  = (et,1.0)

  tmp  = U.subs(Dict([rp1 rp2 rp3]))
  U_et = tmp

#  dudr_vec = Sym.MatrixSymbol("dUdR", 1,3)
#  dudr_vec[1,1] = dudr_er
#  dudr_vec[1,2] = dudr_eθ
#  dudr_vec[1,3] = dudr_et

  U_vec = [U_et; U_er; U_eθ]

# Put it in a Matrix

  Uvec = Sym.simplify(U_vec).as_mutable()
  return Uvec
end 
