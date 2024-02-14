function vibration_mode_selection(n_vibration_modes::Int, Ox::Bool, Oy::Bool, Oz::Bool, Phix::Bool, Nelast::Int, omega2::Vector{Float64}, eigenmodes::Array{Float64,2})

Ndofs = size(eigenmodes,1)

Nelast = size(omega2,1)

P =zeros(Ndofs,4)


Ox == true  && (P[1,1] = 1.0)
Oy == true  && (P[2,2] = 1.0)
Oz == true  && (P[3,3] = 1.0)
Phix == true  && (P[4,4] = 1.0)

Z, p = findmax(abs.(transpose(P)*eigenmodes), dims = 1)

Z = dropdims(Z, dims = 1)

PREC = sqrt(eps(Float64))

phi_index = nonnuls((Z .> PREC).*[i for i in 1:Nelast])

omega2 = omega2[phi_index]
eigenmodes = eigenmodes[:,phi_index]

return omega2[1:n_vibration_modes], eigenmodes[:,1:n_vibration_modes]

end
