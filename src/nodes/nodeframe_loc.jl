"""
    nodeframe_loc
        function constructing the localization of a frame node.
        calling sequence: loc_x = nodeframe_loc(nbr::Int,npar::Int)
"""
function nodeframe_loc(nbr::Int,npar::Int)
    if nbr <= 0
        error("NodeFrame: node number must be > 0")
    else
            loc_x = [(nbr-1)*6+i for i in 1:6]
    end
    npar > 0 &&   for i in 4:6 loc_x[i]= (npar-1)*6+i end

    return loc_x
end
