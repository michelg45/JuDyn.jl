"""
    linear_constraint

        function expressing a linear constraint between a set of dofs.

            The associated data are collected from Main.SetElements.lin_constr_container.

            calling sequence: linear_constraint(nbr,y,Dy,res)
"""
function linear_constraint(nbr::Int,y::Vector,Dy::Vector,res::Vector,matrix::Bool)

    ec =  Main.element_container
    lc = Main.SetElements.lin_constr_container
    cf =  Main.Frames.current_frames

    iel = findfirst(x -> x == nbr, lc.numbers)[1]
    nodes = lc.node_orders[iel]
    components = lc.components[iel]
    coefs = lc.coefs[iel]
    val = lc.val[iel]
    scale_factor = lc.scale_factor[iel]
  

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel2]
    inv_loc_mult = ec.inv_loc_mult[iel2]
    inv_loc = [inv_loc_x;inv_loc_mult]
    Ndim = size(inv_loc,1) 

    S_el = Array{Float64,2}(zeros(Ndim,Ndim))
    res_el  = Vector{Float64}(zeros(Ndim))

    (s,Ds) = pull_vectors(inv_loc,y,Dy)

    Nx = size(inv_loc_x,1)

    x = zeros(Nx)
       
    for i = 1:Nx
        ip = nodes[i]        
        icomp = components[i]
        if icomp  > 3 
             psi = cf[ip].psi
             Dpsi = RV3()
             Dpsi[icomp-3] = Ds[i]
             x[i] = RV3(psi,Dpsi)[icomp-3]
        else
            x[i] = s[i] + Ds[i]
        end
    end

    mult = s[Ndim] + Ds[Ndim]

    S_el[Ndim,1:Nx] = scale_factor*coefs
    S_el[1:Nx,Ndim] = S_el[Ndim,1:Nx]
    res_el[1:Nx] = -scale_factor*mult*coefs
    res_el[Ndim] = scale_factor*(val - coefs'*x)

    if matrix == false 
        push_element_sparse(res,iel2,inv_loc,S_el,res_el)
    else 
        return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end

    return
end

