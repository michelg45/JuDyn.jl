__precompile__()

"""
    Solve

    The module "Solve"  provides the environment to solve diffent types of structural problems: 
        * dynamic transient response,
        * static response,
        * static response under constraints.
"""
module Solve



using LinearAlgebra
using Dates
using JSON
using HDF5
using SparseArrays
using JuMP
using Ipopt
using Arpack
using ..MyAlgebra
using ..Utils
using ..Frames
using ..BuildElements

dir = "solve/"

include(dir*"solve.jl")
include(dir*"create_h5_file_dynamic.jl")
include(dir*"create_h5_file_static.jl")
include(dir*"dynamic_element_forces.jl")
include(dir*"dynamic_element_system.jl")
include(dir*"dynamic_solve.jl")
include(dir*"eig_solve.jl")
include(dir*"eigenvalue_system_dynamic.jl")
include(dir*"sparse_matrix_structure.jl")
include(dir*"null_sparse_matrix.jl")
include(dir*"record_on_h5_file_dynamic.jl")
include(dir*"record_on_h5_file_static.jl")
include(dir*"sparse_matrix_index.jl")
include(dir*"sparse_matrix_prepro.jl")
include(dir*"sparse_submatrix.jl")
include(dir*"static_constrained_solve.jl")
include(dir*"static_element_forces.jl")
include(dir*"static_element_system.jl")
include(dir*"static_solve.jl")
# 

export solve
export create_h5_file_dynamic
export create_h5_file_static
export dynamic_element_forces
export dynamic_element_system
export dynamic_solve
export eig_solve
export eigenvalue_system_dynamic
export sparse_matrix_structure
export null_sparse_matrix
export record_on_h5_file_dynamic
export record_on_h5_file_staic
export sparse_matrix_index
export sparse_matrix_prepro
export sparse_submatrix
export static_constrained_solve
export static_element_forces
export static_element_system
export static_solve


end # module
