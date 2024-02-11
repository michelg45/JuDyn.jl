"""
    set_rigid_mass

Function to define the topology and inertia properties of a rigid mass element. It constructs the  `Main.SetElements.rigid_mass_container[iel]` array entry (`::RigidMassArray` type) collecting the intial data.
    
Calling sequence: 

````{verbatim}
set_rigid_mass(nbr,node,mass)
````

Input data:

| | |
|:------------------------------------- |:--------------------------------------- |
|  `nbr::Int` | number of the element | 
|  `node::Int` | reference node of the rigid body | 
| `mass::Float64` |  mass of the rigid body | 

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