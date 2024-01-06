
__precompile__()

module FORDYN
using JLD
using ..MyAlgebra

"""

    FORDYN is a preprocessor module that builds  JuDyn.FORDYN.model_container, JuDyn.FORDYN.node_container and JuDyn.FORDYN.element_container, data structures containing the general data of the model.

"""

dir = "fordyn/"

include(dir*"ModelArray.jl")

include(dir*"print_model.jl")

include(dir*"set_general_data.jl")

include(dir*"save_model_array.jl")

include(dir*"load_model_array.jl")

include(dir*"create_model.jl")

include(dir*"NodeArray.jl")

include(dir*"ElementArray.jl")

export ModelArray
export NodeArray
export ElementArray
export create_model
export print_model
export save_model_array
export load_model_array
export set_general_data

end # module