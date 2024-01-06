
__precompile__()

"""
    Assemble
"""
module Assemble

using ..Utils

dir = "assemble/"

include(dir*"assemble.jl")
include(dir*"get_struc_loc.jl")
include(dir*"get_initial_configuration.jl")


export assemble
export get_struc_loc
export get_initial_configuration


end # module
