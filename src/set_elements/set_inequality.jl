
"""
    set_inequality

        function defining an inequality constraint at a node with current position (x_0 + u), , which can be described  in two different ways: 
            
            - "point_to_axis" : distance between the constrained node  and an axis defined by the (p1,p2) node set.

                    distance(x_0+u,axis) >= c  or  distance(node,axis) <= c

            - "point_along_direction" : projection 'h' of  the constrained node position  on a direction  'direction'. 

                    h(x,direction)  >= c  or  h(x,direction) <= c
            
            - "linear" :  distance from the constrained node along a fixed direction defined by the "direction" parameter.

                    x_0 + dot(u,direction) >= c  or  x_0 + dot(u,direction) >= c
            
        Input data:

            nbr::Int                        number of the element
            node::Int                       node of appication of the constraint
            p1::Vec3,p2::Vec3               set of nodes to define an axis
            direction::Vec3                 specified direction
            type::String                    "point_to_axis" , "point_along_direction" or "linear"
            bound_type::String              "upper" or "lower"
            time_function::String           time function from 'input_functions' directory.
            bound_params::Vector{Float64}   parameters of the time function.

        Calling sequences: 

            set_inequality(nbr,node,p1,p2,"point_to_axis",bound_type,time_function, bound_params)

            or

            set_inequality(nbr,node,direction,"linear",bound_type,time_function, bound_params)            


"""

function set_inequality(nbr::Int, node::Int, p1::Vec3,p2::Vec3, type::String, bound_type::String,time_function::String, bound_params::Vector{Float64})

    mc=Main.model_container
    nc=Main.node_container


    if mc.Inequalities == 0
        global inequality_container = InequalityArray()
    else
         inequality_container = SetElements.inequality_container
    end

    inc = inequality_container

    direction = Vec3()

    node_order=findfirst(x -> x==node, nc.node_numbers)
    append!(inc.number,nbr)
    append!(inc.node,node)

    push!(inc.points,[p1,p2])
    push!(inc.type,type)
    push!(inc.bound_type,bound_type)
    push!(inc.time_function,time_function)
    push!(inc.bound_params,bound_params)
    push!(inc.direction,direction)
    loc_x = nc.locs[node_order][1:3]
    loc_int=Int[]
    loc_v= Int[]
    loc_mult=Int[]
    mc.Inequalities +=1
    type_el = "inequality"
    append_element(nbr,type_el,[node],loc_x,loc_int,loc_v,loc_mult)

return

end

function set_inequality(nbr::Int, node::Int, direction::Vec3, type::String,bound_type::String, time_function::String, bound_params::Vector{Float64})

    
    
        mc=Main.model_container
        nc=Main.node_container
    
    
        if mc.Inequalities == 0
            global inequality_container = InequalityArray()
        else
             inequality_container = SetElements.inequality_container
        end
    
        inc = inequality_container
    
        direction = 1.0/norm2(direction)*direction
        p = Vec3()
    
        node_order=findfirst(x -> x==node, nc.node_numbers)
        append!(inc.number,nbr)
        append!(inc.node,node)
        push!(inc.direction,direction)
        push!(inc.points,[p])
        push!(inc.type,type)
        push!(inc.bound_type,bound_type)
        push!(inc.time_function,time_function)
        push!(inc.bound_params,bound_params)
    
        loc_x = nc.locs[node_order][1:3]
        loc_int=Int[]
        loc_v= Int[]
        loc_mult=Int[]
        mc.Inequalities +=1
        type_el="inequality"
        append_element(nbr,type_el,[node],loc_x,loc_int,loc_v,loc_mult)
    
    return
    
    end
    
    function set_inequality(nbr::Int, node::Int, component::Int, type::String, bound_type::String,time_function::String, bound_params::Vector{Float64})
    
        direction = Vec3()
        direction[component] = 1.0
        set_inequality(nbr,node,direction,type,bound_type,time_function,bound_params)
    
        return 
    end

"""    function set_inequality(nbr::Int, node::Int, p1::Vec3, type::String, bound_type::String,time_function::String, bound_params::Vector{Float64})
    
        p2 = Vec3()

        type != "distance_to_point" && (error("set_inequality - type must be distance_to_point"))

        set_inequality(nbr,node,p1,p2,type,bound_type,time_function,bound_params)
    
        return 
    end"""

