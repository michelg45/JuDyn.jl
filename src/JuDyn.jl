__precompile__()

"""
    JuDyn

"JuDyn.jl" is  is the main module of the package. It automatically includes all the sub-modules [`Assemble`](@ref), [`BoundaryConditions`](@ref), [`BuildElements`](@ref), [`FORDYN`](@ref), 
[`Frames`](@ref),  [`InitialConditions`](@ref),  [`InputFunctions`](@ref),  [`MyAlgebra`](@ref),  [`Nodes`](@ref),  [`LinearBeam`](@ref), [`SetElements`](@ref),  [`Solve`](@ref), [`Utils`](@ref) of the package and their external dependencies (_LinearAlgebra_, _SparseArrays, _Arpack_, _JSON_, _JLD_, _HDF5_, _Revise_, _Reexport_, etc. ). 

"""
module JuDyn

using Revise
using Reexport
using JLD
using JSON

include("MyAlgebra.jl")

@reexport using .MyAlgebra

include("FORDYN.jl")
@reexport using .FORDYN

include("Nodes.jl")
@reexport using .Nodes

include("Utils.jl")
@reexport using .Utils

include("LinearBeam.jl")
@reexport using .LinearBeam

include("SetElements.jl")
@reexport using .SetElements

include("BoundaryConditions.jl")
@reexport using .BoundaryConditions

include("Assemble.jl")
@reexport using .Assemble

include("InitialConditions.jl")
@reexport using .InitialConditions

include("Frames.jl")
@reexport using .Frames

include("BuildElements.jl")
@reexport using .BuildElements

include("InputFunctions.jl")
@reexport using .InputFunctions

include("Solve.jl")
@reexport using .Solve



end