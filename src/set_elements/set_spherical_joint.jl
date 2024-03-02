"""
    set_frame_link


Function which 

* defines the topology of a spherical connection between nodes 'node1' and 'node2'.
* Defines a set of 6 Lagrange multipliers to express the connection.
* Stores the relative positions  of the connection node into the element entry `SetElements.spherical_joint_container[iel]`. 

Calling squences: 
````{verbatin}
        set_sherical_joint(nbr,node1,node2)

        set_sherical_joint(nbr,node1,pos1, node2, pos2)

        set_sherical_joint(nbr,node1,node2,k)

        set_sherical_joint(nbr,node1,pos1, node2, pos2,k)
````
    
Input:

|               |                                            |    
|:--------------|:---------------------------------------------|
| `nbr` | element number |
| `node1, node2` | nodes connected by the spherical element. |
| `pos1, pos2` | relative positions of the connecting node in node frames 1 and 2. |
| ``k`` | constraint scaling factor. |

"""
function set_spherical_joint(nbr::Int,node1::Int,pos1::Vec3,node2::Int,pos2::Vec3, scale_factor::Float64)

  
        mc = Main.model_container
        nc = Main.node_container
    
    
    
    
        if mc.Frame_links == 0
            global spherical_joint_container=SphericalJointArray()
        else
            spherical_joint_container = SetElements.frame_link_container
        end
    
        sjc = spherical_joint_container
    
        append!(sjc.numbers,nbr)
        append!(sjc.scale_factor,scale_factor)
    
        n1=findfirst(x -> x == node1, nc.node_numbers)[1]
        n2=findfirst(x -> x == node2, nc.node_numbers)[1]
        push!(sjc.node_orders,[n1,n2])
        push!(sjc.relative_positions,[pos1, pos2])

        mc.Spherical_joints  += 1
    
        loc_x = [copy(nc.locs[n1]);copy(nc.locs[n2])]
        loc_mult= collect(mc.max_mult +i for i=1:3)
        mc.max_mult +=3
        loc_int = Int[]
        loc_v = Int[]
        append_element(nbr,"spherical_joint",[n1,n2],loc_x,loc_int,loc_v,loc_mult)
    
    end
    
function set_spherical_joint(nbr::Int,node1::Int,node2::Int)
        
        scale_factor = 1.0
        pos1 = Vec3()
        pos2 = Vec3()
        set_spherical_joint(nbr,node1,pos1,node2,pos2,scale_factor)
    
    end

    function set_spherical_joint(nbr::Int,node1::Int,pos1::Vec3,node2::Int,pos2::Vec3)
        
        scale_factor = 1.0

        set_spherical_joint(nbr,node1,pos1,node2,pos2,scale_factor)
    
    end

    function set_spherical_joint(nbr::Int,node1::Int,node2::Int,scale_factor::Float64)

        pos1 = Vec3()
        pos2 = Vec3()
        set_spherical_joint(nbr,node1,pos1,node2,pos2,scale_factor)
    
    end