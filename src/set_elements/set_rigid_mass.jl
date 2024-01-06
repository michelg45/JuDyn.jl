"""
    set_rigid_mass
    
        function to define the topology and inertia properties of a rigid mass element.

        constructs the  "Main.SetElements.rigid_mass_container" array  (::RigidMassArray type) collecting the intial data.

        Calling sequence: set_rigid_mass(nbr,node,mass)

        Input data:

        nbr::Int                                number of the element
        node::Int                               reference node of the element
        mass::Float64                           mass of the element

"""
function set_rigid_mass(nbr::Int,node::Int,mass::Float64)


    mc = Main.model_container
    nc = Main.node_container

    if mc.Rigid_masses == 0
        global rigid_mass_container=RigidMassArray()
    else
         rigid_mass_container = SetElements.rigid_mass_container
    end

    rmc = rigid_mass_container

    node_order=findfirst(x -> x==node, nc.node_numbers)
    append!(rmc.numbers,nbr)
    append!(rmc.node_orders,node_order)
    append!(rmc.masses,mass)
    loc = nc.locs[node_order][1:3]
    loc_x=copy(loc)
    loc_v= collect(mc.max_v +i for i=1:3)
    println("loc_v ", loc_v)
    loc_mult=Int[]
    loc_int=Int[]
    mc.max_v += 3
    mc.Rigid_masses += 1
    append_element(nbr,"rigid_mass",[node_order],loc_x,loc_int,loc_v,loc_mult)

    return
end