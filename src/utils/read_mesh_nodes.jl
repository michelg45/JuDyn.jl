function read_mesh_nodes(fname)
    f1 = open(fname)
    line = " "
    while line[1] != "\$Nodes"
        line = split(readline(f1), " ")
    end
    line = split(readline(f1), " ")
    Nblocks = parse(Int,line[1])
    Nnodes = parse(Int,line[2])
    NumNodes = zeros(Int,Nnodes)
    loc_node_num = 0
    loc_node_val = 0
    NodeCoordinates = zeros(3,Nnodes)
    for nb = 1:Nblocks
        line = split(readline(f1), " ")
        n_nodes = parse(Int,line[4])
        for i = 1:n_nodes
            loc_node_num += 1
            line = split(readline(f1), " ")
            NumNodes[loc_node_num] = parse(Int,line[1])
        end
        for i = 1:n_nodes
            loc_node_val += 1
            line = split(readline(f1), " ")
            NodeCoordinates[1,loc_node_val] = parse(Float64,line[1])
            NodeCoordinates[2,loc_node_val] = parse(Float64,line[2])
            NodeCoordinates[3,loc_node_val] = parse(Float64,line[3])
        end
    end

    close(f1)

    return NumNodes, NodeCoordinates
end
