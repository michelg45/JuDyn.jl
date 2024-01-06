function print_element(iel::Int)

    ec = JuDyn.FORDYN.element_container

    println("element number: ",ec.element_numbers[iel])
    println("type: ", ec.element_types[iel])
    println(" nodes: ",ec.element_nodes[iel])
    println("localization of  frame variables: ", ec.loc_x[iel])
    println("localization of  internal variables: ", ec.loc_int[iel])
    println("localization of multipliers: ", ec.loc_mult[iel])
    println("localization of velocities: ", ec.loc_v[iel])
    println("inverse localization of  frame variables: ", ec.inv_loc_x[iel])
    println("inverse localization of  internal variables: ", ec.inv_loc_int[iel])
    println("inverse localization of multipliers: ", ec.inv_loc_mult[iel])
    println("inverse localization of velocities: ", ec.inv_loc_v[iel])


end
