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
Assemble
BoundaryConditions
BuildElements
FORDYN
Frames
InitialConditions
InputFunctions
JuDyn
MyAlgebra
Nodes
LinearBeam
SetElements
Solve
Utils
```

## FORDYN

```@docs
FORDYN
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
levicivita
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
RigidBodyArray
RigidMassArray
ShellArray
SuperBeamArray
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
set_rigid_body
set_rigid_mass
set_shell
set_super_beam
```

### Auxiliary functions

```@docs
end_elements
read_beam_properties
read_shell_properties
shell_stiffness_matrix
```

### Linear beam

```@docs
linear_beam_element
linear_beam_element_pure_bending
super_beam_matrix_kernel
beam_mass_2D
beam_stiffness_2D
beam_mass_3D_local
beam_stiffness_3D_local
beam_mass_pure_bending_3D_local
beam_stiffness_pure_bending_3D_local
beam_gyr_3D_local
super_beam_matrix_kernel
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
shell
shell_force
super_beam
super_beam_force
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
increment_frames
init_frames
```




Reference{#reference}