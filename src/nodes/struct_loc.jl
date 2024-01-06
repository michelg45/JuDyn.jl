"""
    struct_loc

        function constructing the "Loc" vector containing
        the structural set of node dofs from the initial data 
        in the "Main.node_container" array. 
        calling sequence: struct_loc()
"""
function struct_loc()
    Main.node_container = Main.node_container
    nbr_nodes=size(Main.node_container.node_numbers,1)
    Loc=Vector{Int}[]
    for i in 1:nbr_nodes
        loc_node = Main.node_container.locs[i]
        for j in 1:size(loc_node,1)
            if loc_node[j] > 0
                a=findfirst(ix->ix==loc_node[j],Loc)
                if typeof(a) == Nothing
                    Loc::Vector{Int}=[Loc;loc_node[j]]
                end
            end
        end
    end
    return Loc
end