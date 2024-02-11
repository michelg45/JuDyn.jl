"""
    set_rigid_body
    
Function to define the topology and inertia properties of a rigid body element. It constructs the  `Main.SetElements.rigid_body_container[iel]` array entry (`::RigidBodyArray` type) collecting the intial data.

Calling sequence: 

````{verbatim}
set_rigid_body(nbr,node,mass,inertia)
````
Input data:

| | |
|:------------------------------------- |:--------------------------------------- |
|  `nbr::Int` | number of the element | 
|  `node::Int` | reference node of the rigid body | 
| `mass::Float64` |  mass of the rigid body | 
|  `inertia::Mat3`  or  `inertia::Vec3` | inertia tensor of the rigid body | 

"""
function set_rigid_body(nbr::Int,node::Int,mass::Float64,inertia::Mat3)

    

    mc = Main.model_container
    nc = Main.node_container


    if mc.Rigid_bodies == 0
        global rigid_body_container=RigidBodyArray()
    else
         rigid_body_container = SetElements.rigid_body_container
    end

    rbc = rigid_body_container

    node_order=findfirst(x -> x==node, nc.node_numbers)
    append!(rbc.numbers,nbr)
    append!(rbc.node_orders,node_order)
    append!(rbc.masses,mass)
    push!(rbc.inertia_tensors,inertia)
    loc = nc.locs[node_order]
    loc_x=copy(loc)
    loc_int=Int[]
    loc_v= collect(mc.max_v +i for i=1:6)
    loc_mult=Int[]
    mc.max_v +=6
    mc.Rigid_bodies +=1
    append_element(nbr,"rigid_body",[node_order],loc_x,loc_int,loc_v,loc_mult)

    return
end

function set_rigid_body(nbr::Int,node::Int,mass::Float64,inertia::Vec3)

    inert=Mat3(diagm(inertia.v))
    return set_rigid_body(nbr,node,mass,inert)
end

