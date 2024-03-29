"""
    print_model

        function printing the general data of the model contained in "Main.model.container".

        calling sequence: print_model(model_container)
"""
function print_model(model)

    println("problem name: ", model.Name)
    println("nbr of nodes: ", model.Nodes)
    println("nbr of elements: ", model.Elements)
    println("nbr of node forces: ", model.Nodal_forces)
    println("nbr of node torques: ", model.Nodal_torques)
    println("nbr of node imposed displacements: ", model.Nodal_imposed_displacements)
    println("nbr of superelements: ", model.SuperElements)
    println("nbr of rigid bodies: ", model.Rigid_bodies)
    println("nbr of rigid masses: ", model.Rigid_masses)
    println("nbr of node_links: ", model.Node_links)
    println("nbr of frame_links: ", model.Frame_links)
    println("nbr of ground_hinges: ", model.Ground_hinges)
    println("nbr of hinges: ", model.Hinges)
    println("nbr of prismatic joints: ", model.Prismatic_joints)
    println("nbr of spherical joints: ", model.Spherical_joints)
    println("nbr of Frame_springs: ", model.Frame_springs)
    println("nbr of beams: ", model.Beams)
    println("nbr of shells: ", model.Shells)
    println("nbr of super_beams: ", model.SuperBeams)
    println("nbr of ground_spring_dampers: ", model.Ground_spring_dampers)
    println("nbr of inequality constraints: ", model.Inequalities)
    println("nbr of linear constraints: ", model.Linear_constraints)
    println("total nbr of dofs: ", model.Ndofs)
    println("nbr of frame dofs: ", model.Ndofs_x)
    println("nbr of internal  dofs: ", model.Ndofs_int)
    println("nbr of kinematic  dofs: ", model.Ndofs_q)
    println("nbr of velocity dofs: ", model.Ndofs_v)
    println("nbr of multipliers: ", model.Ndofs_mult)
    println("highest frame dof nbr: ", model.max_x)
    println("highest internal dof nbr: ", model.max_int)
    println("highest velocity dof nbr: ", model.max_v)
    println("highest multiplier nbr: ", model.max_mult)
    println("structural localization q unknowns: ", model.struc_loc_q)
    println("structural localization v unknowns: ", model.struc_loc_v)
    println("initial shape: ", model.initial_shape)
    println("gravity: ", model.gravity)
    println("uniform rotation: ", model.uniform_rotation)
    println("initial rotation speed: ", model.rotation_speed)
    println("matrix_update: ", model.matrix_update)
    println("init_file: ", model.init_file)
end