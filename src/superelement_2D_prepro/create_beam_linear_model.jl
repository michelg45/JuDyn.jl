function create_beam_linear_model(beam_label::String,Nel::Int,n_SE::Int,JLD_file::String,JSON_file::String)


    length, stiffness_properties, mass_properties = read_beam_properties(JSON_file,beam_label)


    length = length/n_SE
    Nodes = Nel+1
    Ndofs = Nodes*6
    model = LinearStructureModel(Nodes,Ndofs)

    model.node_coordinates[1,:] = [(i-1)*length/Nel for i in 1:Nodes]
    model.node_components .= 6
    model.loc_vector = [(i-1)*6+j for i in 1:Nodes for j in 1:6]
    model.ip_nodes = [(i-1)*6+1 for i in 1:Nodes]
    L = length / Nel
    K_el = beam_stiffness_3D_local(L,stiffness_properties)
    M_el = beam_mass_3D_local(L,mass_properties)
    S_el = beam_gyr_3D_local(L,mass_properties[1])
    for iel = 1:Nel
        locel = [[(iel-1)*6+i for i in 1:6]; [iel*6+i for i in 1:6]]
        model.stiffness_matrix[locel,locel] += K_el
        model.mass_matrix[locel,locel] += M_el
        for i = 1:3
            model.S_matrix[i][locel,locel] += S_el[i]  
        end
    end

    file = jldopen(JLD_file,"w")
    write(file,"initial_model",model)
    close(file)

return

end
