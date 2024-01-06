function read_mesh_elements(fname)
    f1 = open(fname)
    line = " "
    while line[1] != "\$Elements"
        line = split(readline(f1), " ")
    end
    line = split(readline(f1), " ")
    Nblocks = parse(Int,line[1])
    Nelements = parse(Int,line[2])
    NumElements = zeros(Int,Nelements)
    TypeElements = zeros(Int,Nelements)
    loc_el_num = 0
    Connected_nodes = zeros(Int,9,Nelements)
    for nb = 1:Nblocks
        line = split(readline(f1), " ")
        n_el = parse(Int,line[4])
        el_typ = parse(Int,line[3])
        for i = 1:n_el
            loc_el_num += 1
            line = split(readline(f1), " ")
            NumElements[loc_el_num] = parse(Int,line[1])
            TypeElements[loc_el_num] = el_typ
            if el_typ == 1
                Nnodes = 2
            elseif el_typ == 15
                Nnodes = 1
            elseif el_typ == 2
                Nnodes = 3
            elseif el_typ == 3
                Nnodes = 4
            elseif el_typ == 5
                Nnodes = 8
            elseif el_typ == 8
                Nnodes = 3
            elseif el_typ == 9
                Nnodes = 6
            elseif el_typ == 10
                Nnodes = 9
            else
                error("unknown element type:", el_typ)
            end
            Connected_nodes[1:Nnodes, loc_el_num] = [parse(Int,line[i+1]) for i in 1:Nnodes]
        end
    end

    close(f1)

    return NumElements, TypeElements, Connected_nodes
end