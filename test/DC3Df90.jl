function DC3Df90(x::T, y::T, z::T, depth::T, dip::T,
                 L1::T, L2::T, W1::T, W2::T,
                 Δu₁::T, Δu₂::T, Δu₃::T; α::T=2/3) where T <: Float64

    # args
    ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz = [zeros(Float64,1) for i=1:12]
    IRET = zeros(Int32,1)

    # call
    ccall((:dc3d_, "../DC3D.f90.so"), Nothing,
          ( # argtype
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, # input
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} # output
           ),
           α, x, y, z, depth, dip, L1, L2, W1, W2, Δu₁, Δu₂, Δu₃, # input
           ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz, IRET # output
          )

    ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz = map(A->A[1],
    [ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz])

    # return
    #return ux, uy, uz
    return [ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz]
end

function DC3D0f90(x::T, y::T, z::T, depth::T, dip::T,
                  pot1::T, pot2::T, pot3::T, pot4::T; α::T=2/3) where T <: Float64

    # args
    ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz = [zeros(Float64,1) for i=1:12]
    IRET = zeros(Int32,1)


    # call
    ccall((:dc3d0_, "../DC3D.f90.so"), Nothing,
          ( # argtype
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, # input
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32} # output
           ),
           α, x, y, z, depth, dip, pot1, pot2, pot3, pot4, # input
           ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz, IRET # output
          )

    ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz = map(A->A[1],
    [ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz])

    # return
    #return ux, uy, uz
    return [ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz]
end
