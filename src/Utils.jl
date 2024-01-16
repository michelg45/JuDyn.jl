__precompile__()

module Utils
using JLD
using Plots

dir= ("utils/")

include(dir*"print_matrix.jl")
export print_matrix

include(dir*"vector_mean.jl")
export vector_mean

include(dir*"nonnuls.jl")
export nonnuls

include(dir*"save_topology.jl")
include(dir*"save_results.jl")

export save_topology, save_results

include(dir*"read_topology.jl")
include(dir*"read_results.jl")

export read_topology, read_results, read_results_dynamic, read_results_static, read_results_constraints

include(dir*"read_mesh_nodes.jl")
include(dir*"read_mesh_elements.jl")

export read_mesh_nodes, read_mesh_elements

include(dir*"struct_to_dict.jl")
export struct_to_dict

include(dir*"read_string.jl")

export read_string

include(dir*"read_positive_number.jl")
include(dir*"find_positive_integers.jl")
export read_positive_number
export find_positive_integers

include(dir*"read_number.jl")
export read_number

include(dir*"read_eigenvalues.jl")
export read_eigenvalues

include(dir*"waterfall_plot.jl")
export waterfall_plot

include(dir*"cylinder.jl")
export cylinder

include(dir*"read_strain_stresses.jl")
export read_strain_stresses

include(dir*"shape_interpol.jl")
include(dir*"lower_lim.jl")

export shape_interpol,lower_lim

include(dir*"save_shapes.jl")
include(dir*"save_shape.jl")
include(dir*"read_shape.jl")

export save_shape,save_shapes,read_shape

end # module