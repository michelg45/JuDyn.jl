
__precompile__()

"""
    SetElements

    The module "SetElements.jl" initializes the data associated to the diffrent types of elemnts
    (except superelements, managed through the "SuperElments.jl" module).
    It constructs the "element_container" global array (ElementArray structure).

"""
module SetElements

using LinearAlgebra
using ..MyAlgebra
using ..Utils
using ..LinearBeam
using JSON
using JLD

dir = "set_elements/"

include(dir*"set_rigid_body.jl")
include(dir*"RigidBodyArray.jl")

include(dir*"set_rigid_mass.jl")
include(dir*"RigidMassArray.jl")

include(dir*"set_beam.jl")
include(dir*"BeamArray.jl")

include(dir*"set_super_beam.jl")
include(dir*"SuperBeamArray.jl")

include(dir*"set_shell.jl")
include(dir*"ShellArray.jl")
include(dir*"shell_mass_kernel.jl")
include(dir*"shell_stiffness_matrix.jl")

include(dir*"set_node_force.jl")
include(dir*"NodalForceArray.jl")

include(dir*"set_node_link.jl")
include(dir*"NodeLinkArray.jl")

include(dir*"set_frame_link.jl")
include(dir*"FrameLinkArray.jl")

include(dir*"set_ground_hinge.jl")
include(dir*"GroundHingeArray.jl")

include(dir*"set_hinge.jl")
include(dir*"HingeArray.jl")

include(dir*"set_inequality.jl")
include(dir*"InequalityArray.jl")

include(dir*"set_linear_constraint.jl")
include(dir*"LinConstrArray.jl")

include(dir*"set_node_displacement.jl")
include(dir*"NodalDispArray.jl")

include(dir*"set_node_torque.jl")
include(dir*"NodalTorqueArray.jl")

include(dir*"set_prismatic_joint.jl")
include(dir*"PrismaticJointArray.jl")

include(dir*"set_frame_spring.jl")
include(dir*"FrameSpringArray.jl")

include(dir*"set_ground_spring_damper.jl")
include(dir*"GroundSpringDamperArray.jl")

include(dir*"set_super_element.jl")
include(dir*"SEMatrixSet.jl")
include(dir*"read_SE_herting.jl")

include(dir*"append_element.jl")
include(dir*"read_beam_properties.jl")
include(dir*"read_shell_properties.jl")
include(dir*"print_element.jl")
include(dir*"pull_vectors.jl")
include(dir*"find_element_components.jl")
include(dir*"find_element_dof.jl")
include(dir*"end_elements.jl")
include(dir*"SuperElementArray.jl")


export RigidBodyArray
export RigidMassArray
export NodeLinkArray
export FrameLinkArray
export FrameSpringArray
export HingeArray
export InequalityArray
export GroundHingeArray
export GroundSpringDamperArray
export NodalForceArray
export NodalTorqueArray
export NodalDispArray
export PrismaticJointArray
export BeamArray
export ShellArray
export LinConstrArray
export end_elements
export SuperBeamArray
export SEMatrixSet
export SuperElementArray



export set_rigid_body
export set_beam
export set_super_beam
export set_rigid_mass
export set_node_link
export set_node_force
export set_node_torque
export set_node_displacement
export set_frame_link
export set_frame_spring
export set_ground_hinge
export set_hinge
export set_inequality
export set_ground_spring_damper
export set_prismatic_joint
export set_shell
export set_linear_constraint
export shell_stiffness_matrix
export shell_mass_kernel
export set_super_beam
export set_super_element

export append_element
export print_element
export pull_vectors
export find_element_components
export find_element_dof
export read_beam_properties
export read_shell_properties
export read_SE_herting

end # module
