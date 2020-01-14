include("./DC3Df.jl")
include("./DC3Df90.jl")

# -----------------------------
# parameters
# -----------------------------
X = -99.0:2.0:99.0
Y = X
nx = length(X)
ny = length(Y)

L1 = -50.0
L2 = 50.0
W1 = -50.0
W2 = 50.0

z = -30.0
depth = 50.0
dip = 70.0
Δu₁ = 200.0
Δu₂ = -150.0
Δu₃ = 100.0
# -----------------------------

# -----------------------------
# calc
# -----------------------------
Uf = zeros(Float64,ny,nx,12)
Uf90 = zeros(Float64,ny,nx,12)
for i = 1:ny
for j = 1:nx
    Uf[i,j,:] = DC3Df(X[j], Y[i], z, depth, dip, L1, L2, W1, W2, Δu₁, Δu₂, Δu₃)
    Uf90[i,j,:] = DC3Df90(X[j], Y[i], z, depth, dip, L1, L2, W1, W2, Δu₁, Δu₂, Δu₃)
end
end
Δ = Uf90 .- Uf
Δmax = findmax(abs.(Δ))
# -----------------------------

# -----------------------------
# Figure
# -----------------------------
using Plots
gr()
using LaTeXStrings

plts = [plot(X, Y, Δ[:,:,k], linetype=:heatmap, c=:balance,
            axis_ratio=:equal, xlims=extrema(X), ylims=extrema(Y)) for k=1:12]

pltu = plot(plts[1:3]..., clims=(-1e-5,1e-5), title=[L"u_x" L"u_y" L"u_z"], layout=(1,3), size=(900,200))

titlestr = reshape([L"\frac{\partial u_x}{\partial x}" L"\frac{\partial u_y}{\partial x}" L"\frac{\partial u_z}{\partial x}"
                    L"\frac{\partial u_x}{\partial y}" L"\frac{\partial u_y}{\partial y}" L"\frac{\partial u_z}{\partial y}"
                    L"\frac{\partial u_x}{\partial z}" L"\frac{\partial u_y}{\partial z}" L"\frac{\partial u_z}{\partial z}"], (1,9))

pltd = plot(plts[4:12]..., clims=(-1e-6, 1e-6), title=titlestr, layout=(3,3), size=(1000,800))

savefig(pltu,"diffu.svg")
savefig(pltd,"diffu_derivative.svg")
# -----------------------------
