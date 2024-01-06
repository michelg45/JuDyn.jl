mutable struct SEMatrixSet

    Gyro::Bool
    SE_name::String
    SE_type::String
    Nrig::Int
    N_B::Int
    N_I::Int
    Mrig::Array{Float64,2}
    K::Array{Float64,2}
    M::Array{Float64,2}
    S_rot::Vector{Array{Float64,2}}
    A_rot::Vector{Array{Float64,2}}

 
end
