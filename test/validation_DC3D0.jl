include("./DC3Df90.jl")

# -----------------------------
# parameters
# -----------------------------
X = -24.5:1.0:24.5
Y = X
nx = length(X)
ny = length(Y)

z = -30.0
depth = 50.0
dip = 0.0
# -----------------------------

# -----------------------------
# calc
# -----------------------------
Uf90 = zeros(Float64,ny,nx,4)
for i = 1:ny
for j = 1:nx
    Uf90[i,j,1] = DC3D0f90(X[j], Y[i], z, depth, dip, 1e3, 0.0, 0.0, 0.0)[3]
    Uf90[i,j,2] = DC3D0f90(X[j], Y[i], z, depth, dip, 0.0, 1e3, 0.0, 0.0)[3]
    Uf90[i,j,3] = DC3D0f90(X[j], Y[i], z, depth, dip, 0.0, 0.0, 1e3, 0.0)[3]
    Uf90[i,j,4] = DC3D0f90(X[j], Y[i], z, depth, dip, 0.0, 0.0, 0.0, 1e3)[3]
end
end

# -----------------------------

# -----------------------------
# Figure
# -----------------------------
using Plots
#plotlyjs()
pyplot()

plt1 = plot(X, Y, Uf90[:,:,1], linetype=:surface, c=:balance, clims=(-0.20,0.20), zlims=(-0.15,0.15))
plt2 = plot(X, Y, Uf90[:,:,2], linetype=:surface, c=:balance, clims=(-0.20,0.20), zlims=(-0.15,0.15))
plt3 = plot(X, Y, Uf90[:,:,3], linetype=:surface, c=:balance, clims=(-0.50,0.50), zlims=( 0.00,0.55))
plt4 = plot(X, Y, Uf90[:,:,4], linetype=:surface, c=:balance, clims=(-0.10,0.10), zlims=( 0.00,0.10))

plts = plot(plt1, plt2, plt3, plt4, layout=(2,2), size=(1000,600))
# -----------------------------
