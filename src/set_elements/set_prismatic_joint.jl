"""
    set_prismatic_joint


Function defining a prismatic joint between two nodes 'nodes[1]' and 'nodes[2]'.

* It sets the relative motion direction and orientation between the nodes.
* It defines a set of 6 Lagrange multipliers to express the connection.
* It stores the initial joint geometry  and other proerties into the dataset entry `prismatic_joint_container[iel]` (`PrismaticJointArray` type.). 
* The prismatic joint can have different physical properties defined by the `mode::string` variable: 
    * "free" : no actuation force.
    * "force": actuated by a force defined by a time-dependent input function `func::string` and its parameters `params::Vector{Float64}`.
    * "driven": imposed displacement defined by a time-dependent input function `func::string` and its parameters `params::Vector{Float64}`.
    * "spring": viscoelastic pring characterized by a stiffness ``k``, a viscous constant ``c            `` and an initial parameter ``d_0`` collected in vector `params::Vector{Float64}` = ``[k,c,d_0]``.



Calling squences: 
````{verbatin}
        set_prismatic_joint(nbr,nodes,mode,func,params,scale_factor)

        set_prismatic_joint(nbr,nodes)

        set_prismatic_joint(nbr,nodes,scale_factor)

        set_prismatic_joint(nbr,nodes,"spring",params)

        set_prismatic_joint(nbr,nodes,"spring",params,scale_factor)
````
    
Input:

|   |  |
|:--------------|:---------------------------------------------|
| nbr::Int | element number |
| nodes::Vector{Int} | set of 2 nodes connected by the prismatic joint element. |
| mode::string | can take the values "free", "force", "driven", "spring". |
| func::string | time function from `ÌnputFunctions` describing an imposed force or displacement.|
| params::Vector{Float64} | set of parameters defining either the function `func` ("force" and "driven" modes) or  the viscoelastic spring ("spring" mode). | 
| scale_factor::Float64 | scale factor to the constraints (set to 1.0 if absent). |

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
    RV_1 = copy(nc.init_orientations[n1])
    RV_2 = copy(nc.init_orientations[n2])
    RV_line = frame_on_line(x_1,x_2)
    RV_0 = RV3(RV_1,RV_line)
    RV_rel = RV3(-RV_1,RV_2)
    axis = rot(-RV_1,(x_2-x_1))
    axis = 1.0/norm2(axis)*axis
    
    append!(pjc.number,nbr)
    push!(pjc.node_orders,node_orders)
    append!(pjc.l_0,l_0)
    push!(pjc.RV_line,RV_line)
    push!(pjc.mode,mode)
    push!(pjc.RV_0,RV_0)
    push!(pjc.RV_rel,RV_rel)
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
    
    function set_prismatic_joint(nbr::Int, connecting_nodes::Vector{Int})

        mode = "free"
        scale_factor = 1.0
        func = " "
        params = []

        set_prismatic_joint(nbr,connecting_nodes,mode,func, params,scale_factor)

    end

    function set_prismatic_joint(nbr::Int, connecting_nodes::Vector{Int},scale_factor::Float64)

        mode = "free"
        func = " "
        params = []

        set_prismatic_joint(nbr,connecting_nodes,mode,func, params,scale_factor)

    end


    function set_prismatic_joint(nbr::Int, connecting_nodes::Vector{Int},mode::String,params::Vector{Float64}, scale_factor::Float64)

        mode != "free" && (error("set_prismatic_joint calling sequence: mode must be free."))
        func = " "

        set_prismatic_joint(nbr,connecting_nodes,mode,func, params,scale_factor)

    end

    function set_prismatic_joint(nbr::Int, connecting_nodes::Vector{Int},mode::String,params::Vector{Float64})

        mode != "free" && (error("set_prismatic_joint calling sequence: mode must be free."))
        func = " "
        scale_factor = 1.0

        set_prismatic_joint(nbr,connecting_nodes,mode,func, params,scale_factor)

    end