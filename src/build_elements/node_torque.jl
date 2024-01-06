"""
    node_torque 
"""
function node_torque(nbr::Int,Dy::Vector,res::Vector,p::Vector,itime::Int,h::Float64,matrix::Bool)


    tc = SetElements.torque_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, tc.number)[1]

    str_time_function = tc.time_function[iel]

    time_function = input_functions[str_time_function]
    axe = tc.axe[iel]
    frame = tc.frame[iel]
    params = tc.params[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inode = ec.element_nodes[iel2][1]

    psi_n = cf[inode].psin

    psi = RV3(cf[inode].psin, cf[inode].Dpsi)

    frame == "absolute" ?   (axe_n = rot(-psi_n,axe); axe = rot(-psi,axe)) : axe_n = axe

    inv_loc_x= ec.inv_loc_x[iel2]

    torque_n = time_function(itime-1,h,params)*axe_n

    torque = time_function(itime,h,params)*axe

    p[inv_loc_x] += torque.v

    res_el = torque.v

    S_el = zeros(3,3)

    frame == "absolute" ? S_el = -(tilde(torque)*tang(cf[inode].Dpsi)).mat : S_el = zeros(3,3)

    ext_work_el = 0.5*transpose(Dy[inv_loc_x])*(torque + torque_n).v

    if matrix == false 
        push_element_sparse(res,iel2,inv_loc_x,S_el,res_el)
        return ext_work_el
    else 
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc_x,S_el)
    end
  

    end

    function node_torque(nbr::Int,Dy::Vector{Float64},p,itime,h)
        N = size(Dy,1)
        res = zeros(N)
        matrix = true
        println("node_torque")
        return inv_loc, S_el = node_torque(nbr,Dy,res,p,itime,h,matrix)
    end
    
    function node_torque(nbr::Int,Dy::Vector,res::Vector,p::Vector,itime::Int,h::Float64)
        matrix = false
        node_torque(nbr,Dy,res,p,itime,h,matrix)
    end
