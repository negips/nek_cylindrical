#!/usr/bin/julia

#    Derivation of Cylindrical coordinates using symbolic Julia

println("Symbolic Math in Julia using SymPy (Python Package)")
println("Derivation for Cylindrical Coordinates")

using PyCall
using SymPy

include("ReplaceUnitVecDeriv.jl")
include("ReplaceUnitVecDot.jl")
include("GetComponents.jl")
include("GradR.jl")
include("GradTheta.jl")
include("GradT.jl")
include("Grad.jl")
include("Div.jl")


Sym = pyimport("sympy")

# Cartesian Coordinates
x = Sym.Symbol("x")
y = Sym.Symbol("y")
z = Sym.Symbol("z")

# Cylindrical Coordinates
t = Sym.Symbol("t")          # x in the cylindrical coordinates
r = Sym.Symbol("r", positive=true)
θ = Sym.Symbol("θ")

x = t
y = r*cos(θ)
z = r*sin(θ)


# x
dxdr = Sym.diff(x,r)
dxdθ = Sym.diff(x,θ)
dxdt = Sym.diff(x,t)

# y
dydr = Sym.diff(y,r)
dydθ = Sym.diff(y,θ)
dydt = Sym.diff(y,t)

# z
dzdr = Sym.diff(z,r)
dzdθ = Sym.diff(z,θ)
dzdt = Sym.diff(z,t)

row_t = [dxdt  dydt dzdt]
row_r = [dxdr  dydr dzdr]
row_θ = [dxdθ  dydθ dzdθ]

A     = Sym.simplify([row_t; row_r; row_θ])

A_inv = Sym.simplify(Sym.inv_quick(A))
#
dtdx = A_inv[0,0] 
dtdy = A_inv[0,1]  
dtdz = A_inv[0,2]
drdx = A_inv[1,0] 
drdy = A_inv[1,1]  
drdz = A_inv[1,2]
dθdx = A_inv[2,0] 
dθdy = A_inv[2,1]  
dθdz = A_inv[2,2]

# Unit Vectors in Clindrical coordinates
er = Sym.Function("er", commutative=false)(r,θ)
eθ = Sym.Function("eθ", commutative=false)(r,θ)
et = Sym.Function("et", commutative=false)()

# Commutatitive symbols for collection
erc = Sym.Function("erc", commutative=true)(r,θ)
eθc = Sym.Function("eθc", commutative=true)(r,θ)
etc = Sym.Function("etc", commutative=true)()

# Unit Vectors in Cartesian coordinates
ex = Sym.Symbol("ex")
ey = Sym.Symbol("ey")
ez = Sym.Symbol("ez")


er1 = cos(θ)*ex + sin(θ)*ey
eθ1 = -sin(θ)*ex + cos(θ)*ey
et1 = ez

derdθ = Sym.diff(er,θ)
deθdθ = Sym.diff(eθ,θ)  
#replace(derdθ,-ex*sin(θ) + ey*cos(θ),eθ)
#replace(deθdθ,-ex*cos(θ) - ey*sin(θ),-er)

# Tangent Vectors in Clindrical coordinates
Ur = Sym.Function("Ur",commutative=false)(r,θ,t)
Uθ = Sym.Function("Uθ",commutative=false)(r,θ,t)
Ut = Sym.Function("Ut",commutative=false)(r,θ,t)

U_right  = Ur*er + Uθ*eθ + Ut*et
U_left   = er*Ur + eθ*Uθ + et*Ut

#gradU = [dUdr.__pyobject__[0] dUdr.__pyobject__[1] dUdr.__pyobject__[2]; 
#         RinvdUdθ.__pyobject__[0] RinvdUdθ.__pyobject__[1] RinvdUdθ.__pyobject__[2];
#         dUdt.__pyobject__[0] dUdt.__pyobject__[1] dUdt.__pyobject__[2]]

#gradU = Sym.simplify(gradU)

gradU  = Grad(U_right,r,θ,t,"left")

#gradUM = Sym.simplify(gradUM)

#parsed_expr = Sym.parse_expr(gradU,evaluate=False)

Err   = Sym.diff(Ur,r)
Eθr   = r/2*Sym.diff(1/r*Uθ,r) + (1/2)*(1/r)*Sym.diff(Ur,θ)
Etr   = 1/2*Sym.diff(Ur,t) + 1/2*Sym.diff(Ut,r)

Erθ   = copy(Eθr)
Eθθ   = (1/r)*Sym.diff(Uθ,θ) + 1/r*Ur
Etθ   = (1/2)*(1/r)*Sym.diff(Ut,θ) + 1/2*Sym.diff(Uθ,t)

Ert   = copy(Etr)
Eθt   = copy(Etθ)
Ett   = Sym.diff(Ut,t)

Erj   = (er*Err + eθ*Eθr + et*Etr)*2
Eθj   = (er*Erθ + eθ*Eθθ + et*Etθ)*2
Etj   = (er*Ert + eθ*Eθt + et*Ett)*2

divEr = Div(Erj,"left")
divEθ = Div(Eθj,"left")
divEt = Div(Etj,"left")

#VecdivEr = Div(Erj*er,"left")
#VecdivEθ = Div(Eθj*eθ,"left")
#VecdivEt = Div(Etj*et,"left")
#
#alldiv   = VecdivEr + VecdivEθ + VecdivEt

# Lets try

gradur,mgradur = GradR(U_right,r,θ,t,"left") 
graduθ,mgraduθ = GradTheta(U_right,r,θ,t,"left") 
gradut,mgradut = GradT(U_right,r,θ,t,"left") 

mGU = [mgradur[1] mgradur[2] mgradur[3]; 
       mgraduθ[1] mgraduθ[2] mgraduθ[3]; 
       mgradut[1] mgradut[2] mgradut[3]]

#Eij = mGU .+ 0.
#
#for i=1:3, j=1:3
#  global E
#  Eij[i,j]   = 1/2*(mGU[i,j] + mGU[j,i])
#end
#
#

# er = Sym.Function("er", commutative=false)()
# eθ = Sym.Function("eθ", commutative=false)()
# et = Sym.Function("et", commutative=false)()


GjR   = er*mGU[1,1] + eθ*mGU[2,1] + et*mGU[3,1]
GRj   = er*mGU[1,1] + eθ*mGU[1,2] + et*mGU[1,3]

Gjθ   = er*mGU[1,2] + eθ*mGU[2,2] + et*mGU[3,2]
Gθj   = er*mGU[2,1] + eθ*mGU[2,2] + et*mGU[2,3]

Gjt   = er*mGU[1,3] + eθ*mGU[2,3] + et*mGU[3,3]
Gtj   = er*mGU[3,1] + eθ*mGU[3,2] + et*mGU[3,3]

#GRj2  = ReplaceUnitVecDot(GRj)

#ERj   = (er*mGU[1,1] + er*mGU[1,1]) + (eθ*mGU[2,1] + eθ*mGU[1,2]) + (et*mGU[3,1] + et*mGU[1,3])
ER    = GRj + GjR
Eθ    = Gθj + Gjθ
Et    = Gtj + Gjt

#ERj   = ERj/2

#ERj2  = ReplaceUnitVecDot(ERj)

divER = Div(ER*er,"left")
divEθ = Div(Eθ*eθ,"left")
divEt = Div(Et*et,"left")

divall = divER + divEθ + divEt

divcom = GetComponents(divall)

#divGR = Div(GjR,"left")
#divGθ = Div(Gjθ,"left")
#divGt = Div(Gjt,"left")


#divU   = Div(U_left,"left")
#
#dRdivU = Sym.diff(divU,r)
#dθdivU = (1/r)*Sym.diff(divU,θ)
#dtdivU = Sym.diff(divU,t)

#tmp = divERj.subs(Dict([rp1]))

println("Done")








