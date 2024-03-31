|Definition of elements   |Datas tructures | 
|:-----------------------|:--------------------------|
| [`set_beam`](@ref) |  [`BeamArray`](@ref) |
| [`set_frame_link`](@ref) | [`FrameLinkArray`](@ref) | 
| [`set_frame_spring`](@ref) | [`FrameSpringArray`](@ref) |
| [`set_ground_hinge`](@ref) | [`GroundHingeArray`](@ref) |
| [`set_ground_spring_damper`](@ref) |[`GroundSpringDamperArray`](@ref) |
| [`set_hinge`](@ref) | [`HingeArray`](@ref) |
| [`set_linear_constraint`](@ref) | [`LinConstrArray`](@ref) |
| [`set_node_displacement`](@ref) | [`NodalDispArray`](@ref) |
|  [`set_node_force`](@ref) | [`NodalForceArray`](@ref) |
| [`set_node_link`](@ref) | [`NodeLinkArray`](@ref) |
| [`set_node_torque`](@ref) | [`NodalTorqueArray`](@ref) |
| [`set_prismatic_joint`](@ref) | [`PrismaticJointArray`](@ref) |
| [`set_spherical_joint`](@ref) | [`SphericalJointArray`](@ref)|
| [`set_rigid_body`](@ref) | [`RigidBodyArray`](@ref) |
| [`set_rigid_mass`](@ref) | [`RigidMassArray`](@ref) |
| [`set_shell`](@ref) | [`ShellArray`](@ref) |
| [`set_super_beam`](@ref) | [`SuperBeamArray`](@ref) |
| [`set_super_element`](@ref) | [`SEMatrixSet`](@ref) |


### Data structures

```@docs
BeamArray
FrameLinkArray
FrameSpringArray
GroundHingeArray
GroundSpringDamperArray
HingeArray
LinConstrArray
NodalDispArray
NodalForceArray
NodeLinkArray
NodalTorqueArray
PrismaticJointArray
SphericalJointArray
RigidBodyArray
RigidMassArray
ShellArray
SuperBeamArray
SEMatrixSet
```

### Element definition

```@docs
set_beam
set_frame_link
set_frame_spring
set_ground_hinge
set_ground_spring_damper
set_hinge
set_linear_constraint
set_node_displacement
set_node_force
set_node_link
set_node_torque
set_prismatic_joint
set_spherical_joint
set_rigid_body
set_rigid_mass
set_shell
set_super_beam
set_super_element
```

### Auxiliary functions

```@docs
end_elements
read_beam_properties
read_shell_properties
shell_stiffness_matrix
read_SE_herting
```