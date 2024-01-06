I want to make an hyperlink to Nodes

# JuDyn documentation

```@contents
Pages = ["index.md"]
Depth  = 3
```

link to [Nodes](@ref)
link to [MyAlgebra](@ref)

## Modules

```@docs
FORDYN
MyAlgebra
Nodes
SetElements
```

## FORDYN

```@docs
ModelArray
NodeArray
ElementArray
```

## MyAlgebra

```@docs
RV3
Vec3
Mat3
Quat
RP3
NodeFrame
crossp
rot
tang
Dtang
invtang
invrot
tilde
euler_to_RV
exp_SE3
log_SE3
invT_SE3
DinvT_SE3
frame_solve
frame_solve1
compute_strains
Adj
shape_functions_1D
shape_functions_2D
gauss_points
```

- link to [`tilde`](@ref)

## Nodes

```@docs
set_node
end_nodes
append_node
set_node_connection
print_node
nodeframe_loc
struct_loc
find_node_components
select_results_XYZ_nodes
select_results_nodes
```

## SetElements

### Data structures

```@docs
RigidBodyArray
RigidMassArray
NodeLinkArray
FrameLinkArray
FrameSpringArray
HingeArray
GroundHingeArray
GroundSpringDamperArray
NodalForceArray
NodalTorqueArray
NodalDispArray
PrismaticJointArray
BeamArray
ShellArray
```

### Element definition

```@docs
set_rigid_body
set_beam
set_super_beam
set_rigid_mass
set_node_link
set_node_force
set_node_torque
set_node_displacement
set_frame_link
set_frame_spring
set_ground_hinge
set_hinge
set_ground_spring_damper
set_prismatic_joint
set_shell
shell_stiffness_matrix
end_elements
read_beam_properties
read_shell_properties
```

## Model assembly

```@docs
assemble
get_struc_loc
get_initial_configuration
sparse_matrix_structure
```

## Initial conditions

```@docs
initial_conditions
get_initial_position
set_initial_displacement
get_initial_displacements
set_initial_velocity
end_initial_conditions
print_initial_configuration
```

## Element construction

```@docs
rigid_body
rigid_body_force
rigid_mass
rigid_mass_force
node_link
node_link_force
node_force
node_torque
node_displacement
frame_link 
frame_link_force
frame_spring
ground_hinge
ground_hinge_force
hinge
hinge_force
beam
beam_force
visco_beam_dyn
visco_beam_dyn_force
visco_beam_static
visco_beam_static_force
shell
shell_force
visco_shell_dyn
visco_shell_dyn_force
visco_shell_static
visco_shell_static_force
super_beam
ground_spring_damper
ground_spring_damper_force
push_element_sparse
push_matrix_sparse
element_constraints
init_element_constraints
point_to_axis
point_along_direction
distance_to_point
```

## Frames

```@docs
CurrentFrame
init_frames
increment_frames
```



Reference{#reference}