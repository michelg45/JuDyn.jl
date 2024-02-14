
__precompile__()

"""
    FORDYN

The *FORDYN* module creates the 3 main data containers used throuhout the execution of the  *JUDyn* module: 
* _model__container_, a data container  of _ModelArray_ type and   containing general data such as number of nodes, number of elements, number of elements of each type, etc. The list of entries in the _model__container_ array can be read from the source code of `ModelArray.jl`. 

* _node__container_, a data container  of _NodeArray_ type and   containing the data  of  the node set. The list of entries in the _model__container_ array can be read from the source code of `NodeArray.jl`. 

* _element__container_, a data container  of _ElementArray_ type and   containing the data  of element set. The list of entries in the _element__container_ can be read from the source code of `ElementArray.jl`. 

These arrays being _global_, they can be consulted after the execution of the *JUDyn*  module as showh in the examples below.

> Examples

````{verbatim}
julia> model_container.Rigid_bodies
1

julia> node_container.node_numbers
2-element Vector{Int64}:
 10
 20
 
julia> element_container.element_types
2-element Vector{String}:
 "rigid_body"
 "frame_link"
````

"""
module FORDYN
using JLD
using ..MyAlgebra

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