function read_initial_model(JLD_file)
  """
  mat_file is a .MAT file containing the inital (K,M) model
  and associated geometric/topological data.
  """

    file = jldopen(JLD_file)
    initial_model = read(file,"initial_model")

    Nodes = initial_model.Nodes
    Ndofs   = initial_model.Ndofs
    node_coordinates = initial_model.node_coordinates
    node_components = initial_model.node_components
    ip_nodes = initial_model.ip_nodes
    loc_vector = initial_model.loc_vector
    stiffness_matrix = initial_model.stiffness_matrix
    mass_matrix = initial_model.mass_matrix

    close(file)

  return Nodes, Ndofs, node_coordinates, node_components, ip_nodes,
    loc_vector, stiffness_matrix, mass_matrix
  end
