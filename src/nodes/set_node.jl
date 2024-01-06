"""
    set_node(...)
        function allowing to define a new node in the node set of the model.
            Nodes can be of 3 types: 
                "frame"     defined by its postion  and orientation.
                "linked"    has same orientation as its parent "frame" node.
                "point_3D"  has only translation degrees of freedom.
            To each node corresponds an entry in the  Main.node_container::NodeArray.
            The Nodes fields data are: 
                nbr::Int        node number
                type::String    "frame", "linked" or "point_3D"
                npar::Int       parent node of a linked node.
                init_pos::Vec3  initial position
                init_or::Vec3   initial orientation
            All nodes must be defined before defining the elements.
            The function has 4 different intances:

                function set_node(nbr::Int,init_pos::Vec3,init_or::RV3,type::String)
                    Node must be of "frame" or "point_3D" type. if "point_3D", 
                    its rotation dofs are fixed.

                function set_node(nbr::Int,init_pos::Vec3,type::String)
                    Node must be of "frame" or "point_3D" type. if "point_3D" type, 
                    its rotation dofs are fixed. if "frame"  type, its initial
                    orientation is set to 0.

                function set_node(nbr::Int,init_pos::Vec3)
                    Node is of "point_3D" type.

                function set_node(nbr::Int,npar::Int,init_pos::Vec3,type::String)
                    Node must be of "linked" type. its orientation is the same as
                    frame node "npar".
"""    
function set_node(nbr::Int,init_pos::Vec3,init_or::RV3,type::String)
    nc = Main.node_container


    if  (type != "frame") & (type != "point3D")
        error("node ", nbr, " wrong type ",type)
    end
    val = findfirst(x -> x==nbr, nc.node_numbers)
    if   typeof(val) == Nothing
        val=0
    end
    if val > 0
      error("node ", nbr, " already defined")
    end
    npar=0
    loc=nodeframe_loc(nbr,npar)
    if type == "point3D"
        loc[4:6]=[0,0,0]
    end
    append_node(nbr::Int, npar::Int, loc, type::String, init_pos::Vec3, init_or::RV3)
    return
end

# =========================================================================================

function set_node(nbr::Int,init_pos::Vec3,type::String)

    init_or=RV3([0.,0.,0.])
    return set_node(nbr,init_pos,init_or,type)
    end

# =========================================================================================

function set_node(nbr::Int,init_pos::Vec3)

    type="point3D"
    init_or=RV3([0.,0.,0.])
    return set_node(nbr,init_pos,init_or,type)

end


# =========================================================================================

function set_node(nbr::Int,npar::Int,init_pos::Vec3,type::String)
    nc = Main.Main.node_container


    val = findfirst(x -> x==nbr, nc.node_numbers)
    if   typeof(val) == Nothing
        val=0
    end
    if val > 0
        error("node ", nbr, " already defined")
    end
    idx=findfirst(x -> x==npar, nc.node_numbers)
    if (npar > 0) & (typeof(idx) == Nothing)
        error("wrong node definition ", nbr, "parent node does not exist")
    elseif (npar > 0) & (type != "linked")
        error("wrong node definition ", nbr,  " has no parent node")
    else
        init_or=nc.init_orientations[idx]
        loc=nodeframe_loc(nbr,npar)
        println(" set node linked")
        println(nbr, " ", npar, " ", loc)
        append_node(nbr, npar, loc, type, init_pos, init_or)
    end
    return
end