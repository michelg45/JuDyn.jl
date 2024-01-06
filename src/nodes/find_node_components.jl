"""
	find_node_components
		function that finds the position of selected components for node "node"
		the in structural set
			str = "all" -> all components
			str = "pos" -> translation components
			str = "rot" -> translation components
		calling sequence:
			inv_loc = find_node_components(node::Int,str::String)

	Examples

	inv_Loc = find_node_components(2,"pos")
	3-element Vector{Int64}:
 		1
 		2
 		3

	inv_Loc = find_node_components(2,"rot")
	3-element Vector{Int64}:
 		4
		5
		0

"""
function find_node_components(node::Int,str::String)

 	nc = Main.node_container 
	inode=findfirst(ix -> ix == node,nc.node_numbers)
	inv_loc  = nc.inv_locs[inode]
	if str == "all"
		return inv_loc
	elseif str == "pos"
		return inv_loc[1:3]
	elseif str == "rot"
		return inv_loc[4:6]
	else
		error("string ", str, "does not exist")
	end
end

function find_node_components(node::Int,icomp::Int)
	
 	nc = Main.node_container
	inode=findfirst(ix -> ix == node,nc.node_numbers)
	if (icomp > 0 && icomp <= 6)
		return inv_loc  = nc.inv_locs[inode][icomp]
	else
		error("component ", icomp, " does not exist")
	end
end
