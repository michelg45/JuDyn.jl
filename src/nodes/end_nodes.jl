"""
    end_nodes
        function call that closes the node set input (compulsory)
        calling sequence: end_nodes()
"""
function end_nodes()
    nc = Main.node_container
    mc = Main.model_container
    mc.end_of_nodes=true
    mc.max_x=findmax(reduce(vcat, nc.locs))[1]
    return
end