    """
        set_node_displacement
    
        function defining a displacement at node
    
        Input:
    
            nbr::Int = element number
            node::Int = node at which displacement is applied
            axe::Vec3 = displacement direction
            ncomp::Int = node component
            time_function::String = name of the function
                                    defined in Dict "input_functions"
            params::Vector{Float64} = set of parameters for time_function

        Calling sequences: 

            set_node_displacement(nbr,node,axe,time_function,params,scale_factor)  

            set_node_displacement(nbr,node,ncomp,time_function,params,scale_factor) 

        Output: entry in JUDYN.SetElements.disp_container::NodalDispArray
    
    """
function set_node_displacement(nbr::Int, node::Int, axe::Vec3, func::String,params::Vector{Float64},scale_factor::Float64)

   
        mc=Main.model_container
        nc=Main.node_container
    
    
        if mc.Nodal_imposed_displacements == 0
            global disp_container = NodalDispArray()
        else
             disp_container = SetElements.disp_container
        end
    
             dc = disp_container
    node_order=findfirst(x -> x==node, nc.node_numbers)
    findfirst(x -> x == node_order, dc.node_order) !== nothing && (error(" only one displacement can be defined at node ", node))  
    append!(dc.number,nbr)
    append!(dc.node_order,node_order)
    push!(dc.axe,axe)
    push!(dc.scale_factor,scale_factor)
    push!(dc.time_function,func)
    push!(dc.params,params)
    
    loc_x=nc.locs[node_order][1:3]
    loc_int=Int[]
    loc_v= Int[]
    loc_mult= collect(mc.max_mult +i for i=1:3)
    mc.max_mult += 3
    mc.Nodal_imposed_displacements +=1
    append_element(nbr,"node_displacement",[node_order],loc_x,loc_int,loc_v,loc_mult)
    
    return
    
    end
    
    function set_node_displacement(nbr::Int, node::Int, ncomp::Int, func::String,params::Vector{Float64},scale_factor::Float64)
   
    if (ncomp == 0 || ncomp > 3)
        error("wrong displacement definition at node ", node)
    end
    
    axe = Vec3(zeros(3))
    
    axe[ncomp] = 1.0
    
    set_node_displacement(nbr,node,axe,func,params,scale_factor)
    
    end
    