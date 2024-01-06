"""
    set_prismatic_joint

    function defing a prismatic joint.

    Input:

        nbr::Int = element number
        connecting_nodes::Int = cnodes connecteed
    
        time_function::String = name of the function
                                defined in Dict "input_functions"
        params::Vector{Float64} = set of parameters for time_function

    Output: entry in JUDYN.SetElements.prism_container::PrismaticJointArray

"""
function set_prismatic_joint(nbr::Int, connecting_nodes::Vector{Int},mode::String,func::String,params::Vector{Float64},scale_factor::Float64)

        mc=Main.model_container
        nc=Main.node_container
    
    
        if mc.Prismatic_joints == 0
            global prism_container = PrismaticJointArray()
        else
             prism_container = SetElements.prism_container
        end
    
             pjc = prism_container
    
    n1=findfirst(x -> x==connecting_nodes[1], nc.node_numbers)[1]
    n2=findfirst(x -> x==connecting_nodes[2], nc.node_numbers)[1]
    node_orders = [n1;n2]
    x_1=copy(nc.init_positions[n1])
    x_2=copy(nc.init_positions[n2])
    l_0 = norm2(x_2-x_1)
    phi_1 = copy(nc.init_orientations[n1])
    axis = rot(-phi_1,(x_2-x_1))
    axis = 1.0/norm2(axi)s*axis
    
    append!(pjc.number,nbr)
    push!(pjc.node_orders,node_orders)
    append!(pjc.l_0,l_0)
    push!(pjc.axis,axis)
    append!(pjc.scale_factor,scale_factor)
    push!(pjc.time_function,func)
    push!(pjc.params,params)
    
    loc_x = [copy(nc.locs[node_orders[1]]); copy(nc.locs[node_orders[2]])]
    loc_v= Int[]
    if mode == "driven"
        loc_int = Int[]
    else
        loc_int = [(mc.max_int+1)]
        mc.max_int += 1
    end
    loc_mult= collect(mc.max_mult +i for i=1:6)
    mc.max_mult += 6
    mc.Prismatic_joints +=1
    append_element(nbr,"prismatic_joint",node_orders,loc_x,loc_int,loc_v,loc_mult)
    
    return
    
    end
    
    