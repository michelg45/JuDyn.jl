mutable struct  VibrationModel

   Nrig::Int
   Nelast::Int
   node_coordinates_CM_frame::Array{Float64,2}
   Mrig::Array{Float64,2}
   U::Array{Float64,2}
   omega2::Vector{Float64}
   eigenmodes::Array{Float64,2}


   function VibrationModel(Nrig,nbr_vibration_modes,node_coordinates_CM, Mrig,
      U, omega2, eigenmodes)
      


      return new(Nrig, nbr_vibration_modes, node_coordinates_CM, Mrig,
      U, omega2, eigenmodes)

   end

   
 end
