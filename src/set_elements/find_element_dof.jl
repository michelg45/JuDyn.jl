"""
	find_element_dof

		function to find position of  component 'icomp' of type 'dof_type' 
		for element 'el_nbr' in structural set

		calling sequence: 

		position = find_element_dof(el_nbr,dof_type,icomp)

		input:
			el_nbr::Int			element number
			dof_type::String	type of dof: "mult", "int", "vit"
			icomp::Int			number the component

"""
function find_element_dof(el_nbr::Int,dof_type::String,icomp::Int)
	element_container = Main.element_container



	iel=findfirst(ix -> ix == el_nbr,element_container.element_numbers)

	if dof_type == "mult"
		inv_loc_mult = element_container.inv_loc_mult[iel]
		n_mult = size(inv_loc_mult)[1]
		if icomp > n_mult
			error("element ", el_nbr, " has no mult. component ", icomp)
		else
			return position = inv_loc_mult[icomp]
		end
	elseif dof_type == "int"
		inv_loc_int = element_container.inv_loc_int[iel]
		n_int = size(inv_loc_int)[1]
		if icomp > n_int
			error("element ", el_nbr, " has no internal component ", icomp)
		else
			return position = inv_loc_int[icomp]
		end
	elseif dof_type == "vit"
		inv_loc_v = element_container.inv_loc_v[iel]
		n_vit = size(inv_loc_v)[1]
		if icomp > n_vit
			error("element ", el_nbr, " has no velocity component ", icomp)
		else
			return position = inv_loc_v[icomp]
		end
	else
		error("element ", el_nbr, " has no  ",dof_type,  " component ", icomp)
	end


end