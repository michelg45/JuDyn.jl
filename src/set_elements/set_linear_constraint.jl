"""
    set_linear_constraint
  
        function to define a linear constraint

            Input:

            * nbr::Int                  element number
            * nodes::Vector{Int}        nodes involved in the constraint
            * components::Vector{Int}   node components involved in the constraint.
            * coefs::Vector{Float64}    coefficients of the constraint
            * val::Float64              constant rhs of the constraint
            * scale_factor::Float64     scaling factor.

            Calling sequences:
                set_linear_constraint(nbr,nodes,components,coefs,val,scale_factor)
                set_linear_constraint(nbr,nodes,components,coefs,val)


"""
    function set_linear_constraint(nbr::Int,nodes::Vector{Int}, components::Vector{Int},coefs::Vector{Float64},val::Float64,scale_factor::Float64)

    mc = Main.model_container
    nc = Main.node_container
    ec = Main.element_container


    if mc.Linear_constraints == 0
        global lin_constr_container = LinConstrArray()
    else
        lin_constr_container = SetElements.lin_constr_container
    end

    lc = lin_constr_container

    append!(lc.numbers,nbr)
    append!(lc.scale_factor,scale_factor)
    append!(lc.val,val)
    push!(lc.coefs,coefs)

    Nnodes = size(nodes,1)
    Ncomp = size(components,1)
    Ncoefs = size(coefs,1)

    Nav = Int(floor(Nnodes+Ncomp+Ncoefs)/3) 

    (Nnodes != Nav || Ncomp != Nav || Ncoefs != Nav ) && ("nodes, components and coefs must be vectors of same size")

   global dofs = zeros(Int,Nnodes)
   global node_orders = zeros(Int,Nnodes)

    for i = 1:Nnodes
        node = nodes[i]
        icomp= components[i]
        inode = findfirst(ix -> ix == node,nc.node_numbers)
        node_orders[i] = inode
        dofs[i] = nc.locs[inode][icomp]
    end

    push!(lc.dofs,dofs)
    push!(lc.node_orders,node_orders)
    push!(lc.components,components)

    loc_x = copy(dofs)
    loc_mult= [mc.max_mult+1]
    loc_int=Int[]
    loc_v= Int[]
    mc.max_mult += 1
    type_el  = "lin_constr"
    mc.Linear_constraints += 1
    append_element(nbr,type_el,node_orders,loc_x,loc_int,loc_v,loc_mult)

return
    
end

function set_linear_constraint(nbr::Int,nodes::Vector{Int}, components::Vector{Int},coefs::Vector{Float64},val::Float64)

    scale_factor = 1.0
    set_linear_constraint(nbr,nodes,components,coefs,val,scale_factor)

end