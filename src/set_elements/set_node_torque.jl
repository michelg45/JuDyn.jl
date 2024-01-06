
"""
    set_node_torque

    function introducing a torque at a given node. 

        Input:

            nbr::Int = element number
            node::Int = node at which torque is applied
            axe::Vec3 = torque direction
            ncomp::Int = node component
            frame::String = "absolute" or "material"
            time_function::String = name of the function
                                    defined in Dict "input_functions"
            params::Vector{Float64} = set of parameters for time_function

        Output: entry in JUDYN.SetElements.torque_container::NodalTorqueArray

        Calling sequences:

            set_node_torque(nbr,node,axe,frame,time_function,params)
            set_node_torque(nbr,node,ncomp,frame,time_function,params)

"""
function set_node_torque(nbr::Int, node::Int, axe::Vec3, frame::String, func::String,params::Vector{Float64})


    
    
    
        mc=Main.model_container
        nc=Main.node_container
    
    
        if mc.Nodal_torques == 0
            global torque_container = NodalTorqueArray()
        else
             torque_container = SetElements.torque_container
        end
    
        tc = torque_container
    
    findfirst(x -> x==node, tc.node) !== nothing ? error(" only one torque can be defined at node ", node) : node_order=findfirst(x -> x==node, nc.node_numbers)
    append!(tc.number,nbr)
    append!(tc.node,node)
    push!(tc.axe,axe)
    (frame != "absolute" && frame != "material") && (error("torque  element ", nbr, " incorrect frame definition"))
    push!(tc.frame,frame)
    push!(tc.time_function,func)
    push!(tc.params,params)
    
    loc_x=nc.locs[node_order][4:6]
    loc_int=Int[]
    loc_v= Int[]
    loc_mult=Int[]
    mc.Nodal_torques +=1
    append_element(nbr,"node_torque",[node],loc_x,loc_int,loc_v,loc_mult)
    
    return
    
    end
    
    function set_node_torque(nbr::Int, node::Int, ncomp::Int, frame::String, func::String,params::Vector{Float64})
    
    if (ncomp == 0 || ncomp > 3)
        error("wrong torque definition at node ", node)
    end
    
    axe = Vec3(zeros(3))
    
    axe[ncomp] = 1.0
    
    set_node_torque(nbr,node,axe,frame,func,params)
    
    end
    