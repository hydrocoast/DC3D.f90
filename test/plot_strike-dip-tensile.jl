include("./DC3Df90.jl")

# -----------------------------
# parameters
# -----------------------------
X = -98.0:4.0:98.0
Y = X
nx = length(X)
ny = length(Y)

L1 = -50.0
L2 = 50.0
W1 = -50.0
W2 = 50.0

z = -30.0
depth = 50.0
#dip = 90.0
# -----------------------------

# -----------------------------
# calc
# -----------------------------
Uf90s = zeros(Float64,ny,nx,3)
Uf90d = zeros(Float64,ny,nx,3)
Uf90t = zeros(Float64,ny,nx,3)
for i = 1:ny
for j = 1:nx
    Uf90s[i,j,:] = DC3Df90(X[j], Y[i], z, depth, 90.0, L1, L2, W1, W2, 1e2, 0.0, 0.0)[1:3] # Δu₁ = 100.0 (cm)
    Uf90d[i,j,:] = DC3Df90(X[j], Y[i], z, depth, 45.0, L1, L2, W1, W2, 0.0, 1e2, 0.0)[1:3] # Δu₂ = 100.0 (cm)
    Uf90t[i,j,:] = DC3Df90(X[j], Y[i], z, depth, 90.0, L1, L2, W1, W2, 0.0, 0.0, 1e2)[1:3] # Δu₃ = 100.0 (cm)
end
end
# -----------------------------

# -----------------------------
# Figure
# -----------------------------
using Plots
plotlyjs()
#pyplot()

xp = permutedims(repeat(X, inner=(1,ny)),[2 1]) # km
yp = repeat(Y, inner=(1,nx))

# Note: in x and y coordinates, values of the dislocation u are thousand times magnified

# Strike-slip
u = 1e-2.*Uf90s # cm -> m
x = xp .+ u[:,:,1]
y = yp .+ u[:,:,2]
z = u[:,:,3]
plts = plot(x, y, z, linetype=:surface, c=:balance, xlabel="(km)", ylabel="(km)", zlabel="(m)",
            zlims=(-0.025,0.025), clims=(-0.02,0.02), colorbar=false, guidefont=10, titlefont=font("Times New Roman",16), title="Strike-slip, δ=90°")

# Dip-slip
u = 1e-2.*Uf90d
x = xp .+ u[:,:,1]
y = yp .+ u[:,:,2]
z = u[:,:,3]
pltd = plot(x, y, z, linetype=:surface, c=:balance, xlabel="(km)", ylabel="(km)", zlabel="(m)",
            zlims=(-0.75,0.75), clims=(-0.5,0.5), colorbar=false, guidefont=10, titlefont=font("Times New Roman",16), title="Dip-slip, δ=45°", viewangle=(130,3))

# Tensile
u = 1e-2.*Uf90t
x = xp .+ u[:,:,1]
y = yp .+ u[:,:,2]
z = u[:,:,3]
pltt = plot(x, y, z, linetype=:surface, c=:balance, xlabel="(km)", ylabel="(km)", zlabel="(m)",
            zlims=(0.0,0.20), clims=(0.0,0.15), colorbar=false, guidefont=10, titlefont=font("Times New Roman",16), title="Tensile, δ=90°")

# -----------------------------
