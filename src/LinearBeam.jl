__precompile__()

"""
    LinearBeam
"""
module LinearBeam

using LinearAlgebra
using ..MyAlgebra
using ..Utils

dir = "linear_beam/"

include(dir*"bar_mass_2D.jl")
include(dir*"bar_stiffness_2D.jl")
include(dir*"beam_mass_2D.jl")
include(dir*"beam_stiffness_2D.jl")
include(dir*"beam_pure_bending_mass_3D_local.jl")
include(dir*"beam_pure_bending_stiffness_3D_local.jl")
include(dir*"beam_mass_3D_local.jl")
include(dir*"beam_stiffness_3D_local.jl")
include(dir*"linear_beam_element.jl")
include(dir*"linear_beam_element_pure_bending.jl")
include(dir*"beam_gyr_3D_local.jl")
include(dir*"super_beam_matrix_kernel.jl")

# ===================================================================================



# ============================================================================

export linear_beam_element
export linear_beam_element_pure_bending
export super_beam_matrix_kernel
export bar_mass_2D
export bar_stiffness_2D
export beam_mass_2D
export beam_stiffness_2D
export beam_mass_3D_local
export beam_stiffness_3D_local
export beam_mass_pure_bending_3D_local
export beam_stiffness_pure_bending_3D_local
export beam_gyr_3D_local
export super_beam_matrix_kernel

end # module
