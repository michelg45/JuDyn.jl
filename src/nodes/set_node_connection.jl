"""
    set_node_connection

    function allowing identity connection between specific dofs 
    of nodes "n1" and "n2". The connected dofs are  specified in
    the 6x1  array 'bool_set': 1 = connected, 0 = not connected.
    calling sequence: set_node_connection(n1::Int,n2::Int,bool_set::Array{Int,1})

    Example

    Definition of a spherical joint between "node1" and "node2"

    set_node_connection(node1,node2,[1,1,1,0,0,0])

"""
function set_node_connection(n1::Int,n2::Int,bool_set::Array{Int,1})

nc = Main.Main.node_container
if Main.model_container.end_of_nodes == true
error("node connection not allowed after closure of node input")
end
if size(bool_set,1) != 6
error("node_connection : incorrect size of Boolean set")
end
node_order_1 = findfirst(x -> x== n1, nc.node_numbers)
node_order_1 === nothing ? error("node_connexion : undefined node ", n1) : loc1 = nc.locs[node_order_1]
node_order_2 = findfirst(x -> x== n2, nc.node_numbers)
if node_order_2 === nothing
error("node_connection : undefined node ", n2)
else
for ncomp = 1:6
if bool_set[ncomp] == 1
  nc.locs[node_order_2][ncomp] = loc1[ncomp]
  nc.inv_locs[node_order_2][ncomp] = loc1[ncomp]
end
end
end
end