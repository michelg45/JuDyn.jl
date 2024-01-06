__precompile__()

"""
    BoundaryConditions
"""
module BoundaryConditions

using LinearAlgebra
using ..MyAlgebra
using ..Utils

dir = "boundary_conditions/"

include(dir*"set_node_BC.jl")
include(dir*"end_BC.jl")


export set_node_BC
export end_BC

end # module
