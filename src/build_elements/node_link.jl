"""
    node_link
"""
function node_link(nbr::Int,Dy::Vector{Float64},res::Vector{Float64},matrix::Bool)

    S_el = Array{Float64,2}(zeros(12,12))
    res_el  = Vector{Float64}(zeros(12))

    nlc = Main.SetElements.node_link_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, nlc.numbers)[1]
    inode1 = nlc.node_orders[iel][1]
    inode_ref = nlc.node_orders[iel][2]
    X = nlc.relative_position[iel]
    k = nlc.scale_factor[iel]

    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]

    inv_loc_mult = ec.inv_loc_mult[iel]
    inv_loc = [inv_loc_x;inv_loc_mult]

    s = Vec3(pull_vectors(inv_loc,Dy))

    

    ix_1 =[i for i in 1:3]
    ix_ref = ix_1 .+ 3
    i_theta = ix_ref .+3
    i_mult = i_theta .+ 3

    mult = Vec3(s[i_mult])
    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_ref = cf[inode_ref].x + Vec3(s[ix_ref])
    Dpsi = RV3(s[i_theta])
    psi = RV3(cf[inode_ref].psi,Dpsi)

    T = tang(Dpsi)

    res_el[ix_1] = k*mult.v
    rtmult=rot(-psi,mult)
    res_el[i_theta] = - k*(crossp(X,rtmult)).v
    res_el[ix_ref] = -k*mult.v
    res_el[i_mult] = k*( x_1 - rot(psi,X) -x_ref).v


    S_el[ix_1,i_mult] = -k*eye(3)
    S_el[i_mult,ix_1] = - k*eye(3)
    S_el[ix_ref,i_mult ]=  k*eye(3)
    S_el[i_mult,ix_ref] =  k*eye(3)
    S_el[i_theta,i_theta] =  (k*tilde(X)*tilde(rtmult)).mat
    A =   tilde(X)*rot(-psi)

    S_el[i_theta,i_mult] = k*A.mat
    S_el[i_mult,i_theta] = k*(transpose(A)*T).mat
    inv_loc = [inv_loc_x; inv_loc_mult]

    if matrix == false 
        push_element_sparse(res,iel,inv_loc,S_el,res_el)
    else 
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc,S_el)
    end

end

"""
    node_link_force
"""
function node_link_force(nbr::Int,Dy::Vector{Float64},res::Vector{Float64})

    Ndim = 12
    res_el  = Vector{Float64}(zeros(Ndim))

    nlc = Main.SetElements.node_link_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container


    iel = findfirst(x -> x == nbr, nlc.numbers)[1]
    inode1 = nlc.node_orders[iel][1]
    inode_ref = nlc.node_orders[iel][2]
    X = nlc.relative_position[iel]
    k = nlc.scale_factor[iel]

    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]
    inv_loc_mult = ec.inv_loc_mult[iel]
    inv_loc = [inv_loc_x;inv_loc_mult]

    s = Vec3(pull_vectors(inv_loc,Dy))

    ix_1 =[i for i in 1:3]
    ix_ref = ix_1 .+ 3
    i_theta = ix_ref .+3
    i_mult = i_theta .+ 3

    mult = Vec3(s[i_mult])

    x_1 = cf[inode1].x + Vec3(s[ix_1])
    x_ref = cf[inode_ref].x + Vec3(s[ix_ref])
    Dpsi = RV3(s[i_theta])
    psi = RV3(cf[inode_ref].psi,Dpsi)

    res_el[ix_1] = k*mult.v
    rtmult=rot(-psi,mult)
    res_el[i_theta] = - k*(crossp(X,rtmult)).v
    res_el[ix_ref] = -k*mult.v
    res_el[i_mult] = k*( x_1 - rot(psi,X) -x_ref).v

    inv_loc = [inv_loc_x; inv_loc_mult]

    for i = 1:Ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] = +res_el[i])
    end

    return

end

