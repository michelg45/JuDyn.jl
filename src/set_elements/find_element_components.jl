"""
    find_node_components

        function to find the position in the structural set of selected components for node 'node' of element 'element' in structural set
            str1 = "frame" -> Frame components
            str1 = "velocities" -> velocity components
            str2 = "all" -> all components
            str2 = "pos" -> translation  components
            str2 = "rot" -> rotation components

        calling sequence: 

            inv_loc = find_element_components(el_nbr,node_nbr,str1,str2)


"""
function find_element_components(element::Int,node::Int,str1::String,str2::String)



	element_container = Main.SetElements.element_container

	iel = findfirst(ix -> ix == element, element_container.element_numbers)
    # println(iel)
    inode = findfirst(ix -> ix == node, element_container.element_nodes[iel])
    # println(inode)
    idx = [i for i in 1:6] .+(inode-1)*6

    if str1 == "velocities"
	    inv_loc  = element_container.inv_loc_v[iel][idx]
    elseif str1 == "frame"
        inv_loc  = element_container.inv_loc_x[iel][idx]
    else
        error("string ", str1, "does not exist")
	end
	if str2 == "all"
		return inv_loc
	elseif str2 == "pos"
		return inv_loc[1:3]
	elseif str2 == "rot"
		return inv_loc[4:6]
	else
		error("string ", str2, "does not exist")
	end

end

function find_element_components(element::Int,str1::String,str2::String)

	"""
	function find_node_components(element::Int,node::Int,str1::String,str2::String)

		finds position of internal degrees of freedom of  element 'element' in structural set
            str1 = "internal" -> internal degrees of freedom
	            str2 = "states" -> all state components
                str2 = "velocities" -> all velocity   components

	"""

	element_container = element_container

if str1 != "internal"
    error("string str1 ", str1, " not equal to internal ")
end

	iel = findfirst(ix -> ix == element, element_container.element_numbers)
    if element_container.n_int[iel] == 0
        error("element ", element, " has no internal degrees of freedom")
    end

    if str2 == "states"
	    inv_loc  = element_container.inv_loc_int[iel]
    elseif str2 == "velocities"

        n_nodes = size(element_container.element_nodes[iel],1)
        n_int = element_container.n_int[iel]
        n_x = element_container.n_x[iel]
        idx= [i for i in 1:n_int] .+ n_x
        inv_loc  = element_container.inv_loc_v[iel][idx]
    else
        error("string ", str2, " does not exist")
	end

    return inv_loc
end

function find_element_components(element::Int,str::String)

	"""
	function find_node_components(element::Int,node::Int,str::String)

		finds position of multipliers of  element 'element' in structural set
            str = "mult" -> multipliers

	"""

	element_container = Main.element_container

if str != "mult"
    error("string str ", str, " not equal to mult ")
end

	iel = findfirst(ix -> ix == element, element_container.element_numbers)
    if element_container.n_mult[iel] == 0
        error("element ", element, " has no multipliers")
    end

    return inv_loc  = element_container.inv_loc_mult[iel]

end
