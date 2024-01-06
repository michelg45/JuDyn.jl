__precompile__()

"""
    InitialConditions

    Module that builds the set of initial conditions
    of the problem under analysis.

"""
module InitialConditions

using ..MyAlgebra
using ..Utils

dir = "initial_conditions/"

include(dir*"initial_conditions.jl")
include(dir*"get_initial_position.jl")
include(dir*"set_initial_displacement.jl")
include(dir*"get_initial_displacements.jl")
include(dir*"set_initial_velocity.jl")
include(dir*"end_initial_conditions.jl")
include(dir*"print_initial_configuration.jl")

export initial_conditions
export get_initial_position
export set_initial_displacement
export get_initial_displacements
export set_initial_velocity
export end_initial_conditions
export print_initial_configuration


end # module
