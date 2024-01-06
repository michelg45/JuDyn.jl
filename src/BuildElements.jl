__precompile__()
"""
BuildElements

    Module constructing the  residual vectors and iteration matrices of the elements.

    The elements included are:

        rigid_mass
        rigid_body
        beam
        super_beam
        node_link
        node_force
        frame_link
        frame_spring
        ground_hinge
        hinge
        ground_spring_damper

    (frame_spring and hinge have not yet been fully tested).

"""
module BuildElements

using VMLS
using LinearAlgebra
using ..MyAlgebra
using ..Utils
using ..Nodes
using ..SetElements
using JLD

include("build_elements/rigid_mass.jl")

include("build_elements/rigid_body.jl")

include("build_elements/beam.jl")

include("build_elements/shell.jl")

include("build_elements/super_beam.jl")

include("build_elements/node_link.jl")

include("build_elements/node_force.jl")

include("build_elements/node_torque.jl")

include("build_elements/node_displacement.jl")

include("build_elements/frame_link.jl")

include("build_elements/frame_spring.jl")

include("build_elements/ground_hinge.jl")

include("build_elements/hinge.jl")

include("build_elements/ground_spring_damper.jl")

include("build_elements/push_element_sparse.jl")
include("build_elements/push_matrix_sparse.jl")
include("build_elements/condensed_element_matrix.jl")

include("build_elements/element_constraints.jl")
include("build_elements/linear_constraint.jl")
include("build_elements/linear_constraint_force.jl")
include("build_elements/init_element_constraints.jl")
include("build_elements/point_to_axis.jl")
include("build_elements/point_along_direction.jl")
include("set_elements/distance_to_point.jl")


export rigid_body
export rigid_body_force
export rigid_mass
export rigid_mass_force
export node_link
export node_link_force
export node_force
export node_torque
export node_displacement
export frame_link 
export frame_link_force
export frame_spring
export ground_hinge
export ground_hinge_force
export hinge
export hinge_force
export beam
export beam_force
export shell
export shell_force
export super_beam
export ground_spring_damper
export ground_spring_damper_force
export push_element_sparse
export push_matrix_sparse
export element_constraints
export linear_constraint,linear_constraint_force
export init_element_constraints
export point_to_axis
export point_along_direction
export distance_to_point


end # module
