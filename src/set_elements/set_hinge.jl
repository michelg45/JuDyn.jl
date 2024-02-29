"""
    set_hinge

        Function defining the topology of a  hinge connection between nodes 'n1' and 'n2'
        Takes into account the relative position  and orientation of the  rotation axis with respect to n1 and n2
        Defines a set of 6 Lagrange multipliers
        constructs the array hinge_container[iel]

        Calling sequence: 
            set_hinge(nbr::Int,node::Vector{Int},pos::Vector{Vec3},or::Vector{RV3},axe::Vec3,mode::String,func::String)

"""
function set_hinge(nbr::Int,node1::Int,node2::Int,pos1::Vec3,or1::RV3,pos2::Vec3,or2::RV3,axe::Vec3,scale_factor::Float64,mode::String,func::String,params::Vector{Float64})
    

    
        mc = Main.model_container
        nc = Main.node_container
    
        if mc.Hinges == 0
            global hinge_container = HingeArray()
        else
             hinge_container = SetElements.hinge_container
        end
    
        hc = hinge_container
        
        append!(hc.numbers,nbr)
        append!(hc.scale_factor,scale_factor)
    
        inodes = [findfirst(x -> x==node1, nc.node_numbers)[1];findfirst(x -> x==node2, nc.node_numbers)[1]]
    
        push!(hc.node_orders,inodes)
    
        axis=RV3(1.0/norm(axe.v)*axe.v)
    
    
        push!(hc.positions,[pos1; pos2])
        push!(hc.orientations,[or1; or2])
        push!(hc.axis,axis)
        push!(hc.mode,mode)
        push!(hc.time_function,func)
        push!(hc.params,params)
    
    
        mc.Hinges  += 1
    
        loc_x = [copy(nc.locs[inodes[1]]); copy(nc.locs[inodes[2]])]
        loc_v = Int[]
    
    
        if mode == "driven"
            loc_int = Int[]
        else
            loc_int = [(mc.max_int+1)]
            mc.max_int += 1
        end
        loc_mult= collect(mc.max_mult +i for i=1:6)
        mc.max_mult += 6
        append_element(nbr,"hinge",inodes,loc_x,loc_int,loc_v,loc_mult)
    
    end
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,axe::Vec3,scale_factor::Float64,mode::String,func::String,params::Vector{Float64})
    
    """
    Hinge with no offset
    """
        pos_0 = Vec3(zeros(3))
        rot_0 = RV3(zeros(3))
    
        set_hinge(nbr,node1,node2,pos_0,rot_0,pos_0,rot_0,axe,scale_factor,mode,func,params)
    
    end
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,axe::Vec3,mode::String,func::String,params::Vector{Float64})
        scale_factor = 1.0
        set_hinge(nbr,node1,node2,axe,scale_factor,mode,func,params)
    end
    
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,axe::Vec3,scale_factor::Float64)
    
    """
    Free hinge with no offset
    """
        mode = "force"
        func = "null_torque"
        params = [0.0]
    
        set_hinge(nbr,node1,node2,axe,scale_factor,mode,func,params)
    
    end
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,axe::Vec3)
        scale_factor = 1.0
        set_hinge(nbr,node1,node2,axe,scale_factor)
    end
    
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,pos1::Vec3,or1::RV3,pos2::Vec3,or2::RV3,axe::Vec3,scale_factor::Float64)
    
    """
    Free hinge with  offset
    """
        mode = "force"
        func = "null_torque"
        params = [0.0]
    
        set_hinge(nbr,node1,node2,pos1,or1,pos2,or2,axe,scale_factor,mode,func,params)
    
    end
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,pos1::Vec3,or1::RV3,pos2::Vec3,or2::RV3,axe::Vec3)
        scale_factor = 1.0
        set_hinge(nbr,node1,node2,pos1,or1,pos2,or2,axe,scale_factor)
    end
    
    function set_hinge(nbr::Int,node1::Int,node2::Int,pos1::Vec3,or1::RV3,pos2::Vec3,or2::RV3,axe::Vec3,mode::String,func::String,params::Vector{Float64})
        scale_factor = 1.0
        set_hinge(nbr,node1,node2,pos1,or1,pos2,or2,axe,scale_factor,mode,func,params)
    end

    function set_hinge(nbr::Int,node1::Int,node2::Int,axe::Vec3,mode::String,params::Vector{Float64})

        scale_factor = 1.0
        func = " "
        pos = Vec3()
        or = RV3()

        set_hinge(nbr,node1,node2,pos,or,pos,or,axe,scale_factor,mode,func,params)
    end