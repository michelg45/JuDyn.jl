"""
    set_node_force

        function to define a time-dependent forcing function in a prescribed direction
            
        calling sequences:

            set_node_force(nbr,node,axe,func,params)

            or

            set_node_force(nbr,node,ncomp,func,params)

        Input:

            nbr::Int = element number
            node::Int = node at which force is applied
            axe::Vec3 = force direction
            or
            ncomp::Int = nodal component (1, 2 or 3)
            time_function::String = name of the function
                                defined in Dict "input_functions"
            params::Vector{Float64} = set of parameters for time_function

        Output: entry in JUDYN.SetElements.force_container::NodalForceArray

"""
function set_node_force(nbr::Int, node::Int, axe::Vec3, func::String,params::Vector{Float64})

    mc=Main.model_container
    nc=Main.node_container


    if mc.Nodal_forces == 0
        global force_container = NodalForceArray()
    else
         force_container = SetElements.force_container
    end

         fc = force_container
    node_order=findfirst(x -> x==node, nc.node_numbers)
    findfirst(x -> x == node_order, fc.node_order) !== nothing && (error(" only one force can be defined at node ", node))  
    append!(fc.number,nbr)
    append!(fc.node_order,node_order)
    push!(fc.axe,axe)
    push!(fc.time_function,func)
    push!(fc.params,params)

    loc_x=nc.locs[node_order][1:3]
    loc_int=Int[]
    loc_v= Int[]
    loc_mult=Int[]
    mc.Nodal_forces +=1
    append_element(nbr,"node_force",[node_order],loc_x,loc_int,loc_v,loc_mult)

return

end

function set_node_force(nbr::Int, node::Int, ncomp::Int, func::String,params::Vector{Float64})

    if (ncomp == 0 || ncomp > 3)
        error("wrong force definition at node ", node)
    end

    axe = Vec3(zeros(3))

    axe[ncomp] = 1.0

    set_node_force(nbr,node,axe,func,params)

end
