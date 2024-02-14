function herting_reduced_model(JLD_file,SE_name, boundary_nodes,boundary_node_components)

    file = jldopen(JLD_file)
    initial_model = read(file,"initial_model")
    vibration_model = read(file,"vibration_model")
    close(file)


    stiffness_matrix = initial_model.stiffness_matrix
    mass_matrix = initial_model.mass_matrix
    S_matrix = initial_model.S_matrix
    Ndofs   = initial_model.Ndofs
    Nodes = initial_model.Nodes
    ip_nodes = initial_model.ip_nodes
    node_components = initial_model.node_components


    Mrig = vibration_model.Mrig
    Nrig = vibration_model.Nrig
    
    node_coordinates_CM_frame = vibration_model.node_coordinates_CM_frame

    X_B = node_coordinates_CM_frame[:,boundary_nodes]
    U = vibration_model.U
    omega2 = vibration_model.omega2
    eigenmodes = vibration_model.eigenmodes
    N_I = vibration_model.Nelast

    N_bnodes = size(boundary_nodes, 1)
    N_B = sum(boundary_node_components)


    loc_boundary_dofs = Int[]



    for i = 1:N_bnodes
        b_node = boundary_nodes[i]
        node_components = boundary_node_components[i]
        loc_boundary_dofs = [loc_boundary_dofs; [ip_nodes[b_node]-1+j for j in 1:node_components]]
    end
        loc_internal_dofs  = setdiff([i for i in 1:Ndofs], loc_boundary_dofs)

    #
    # Computation of attachment modes
    #

    P = zeros(Ndofs+Nrig, N_B)
    P[loc_boundary_dofs,:] = eye(N_B)

    scal = norm(omega2)
    MU = mass_matrix*U
    K_extended = [stiffness_matrix scal*MU; (scal*MU)'  zeros(Nrig,Nrig)]

    G = (K_extended\P)[1:Ndofs,:]

    #
    # Deflation of rigid body modes
    #

    G = G - U*(Mrig\MU')*G


    #
    # scaling of  attachment modes
    #

    G = G *diagonal(1.0./(sqrt.(diag(G'*mass_matrix*G))))

    #
    # Construction of modal transformation matrix V
    #

    P = G[loc_internal_dofs,:]*(inv(G[loc_boundary_dofs,:]))

    V = [eye(N_B); P]

    if N_I > 0
        println("nbr_vibration_modes" ,N_I)
        println("nbr_boundary_dofs ", N_B)
        println("size(eigenmodes) ", size(eigenmodes))
        println("size(P) ", size(P))
        println("size(loc_internal_dofs) ", size(loc_internal_dofs))
        println("size(loc_boundary_dofs) ", size(loc_boundary_dofs))
        println("size(P) ", size(P))
        T = [zeros(N_B,N_I);  (eigenmodes[loc_internal_dofs,:] - P*eigenmodes[loc_boundary_dofs,:])]
        V = [V T]
    end

    new_loc = [loc_boundary_dofs; loc_internal_dofs]

    reduced_stiffness = V'*stiffness_matrix[new_loc,new_loc]*V
    reduced_mass = V'*mass_matrix[new_loc,new_loc]*V
    reduced_S = Matrix{Float64}[]
    push!(reduced_S,V'*S_matrix[1][new_loc,new_loc]*V)
    push!(reduced_S,V'*S_matrix[2][new_loc,new_loc]*V)
    push!(reduced_S,V'*S_matrix[2][new_loc,new_loc]*V)

    println("X_B", X_B)

    write_SE_herting(JLD_file,SE_name,N_bnodes, Nrig, N_B, N_I, boundary_node_components,
    X_B, Mrig, reduced_stiffness, reduced_mass, reduced_S)

end
