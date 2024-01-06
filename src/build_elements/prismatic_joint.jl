"""
    prismatic_joint
"""
function prismatic_joint(nbr::Int,Dy::Vector,y_n,res::Vector,itime::Int,h::Float64,matrix::Bool)


    pjc = SetElements.prism_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, pjc.number)[1]

    str_time_function = pjc.time_function[iel]

    params = pjc.params[iel]
    scale_factor = pjc.scale_factor[iel]
    l_0 =  pjc.l_0[iel]
    axis = pjc.axis[iel]
    nodes =  pjc.node_orders[iel]
    

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    mode == "driven" ? ndim = 18 : ndim =19


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    ext_work_el = 0.0

    nodes = pjc.connecting_nodes[iel]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult= ec.inv_loc_mult[iel2]
    inv_loc_int = ec.inv_loc_int[iel2]

    mode == "driven" ? inv_loc_s = inv_loc_mult : inv_loc_s = [inv_loc_int; inv_loc_mult]

    inv_loc = [inv_loc_x; inv_loc_s]

    disp = input_functions[str_time_function](itime,h,params)


    if  mode != "driven"
        s = pull_vectors(inv_loc_s, Dy)
        if niter == 1
             Dy[inv_loc_int] .= 0.0
             y_n[inv_loc_int] .= 0.0
        end
    else
        s_n, s = pull_vectors(inv_loc_s, y_n, Dy)
    end

    i_x = [i for i=1:12]
    i_x1 = i_x[1:3]
    i_x2 = i_x[7:9]
    i_theta1 = i_x[4:6]
    i_theta2 = i_x[10:12]

    if mode != "driven"
        i_int = 13
        imult= [i for i=14:19]
        x_int = s[1]
        mult = s[2:7]
    else
        imult = [i for i=13:18]
        mult = s[1:6]
        mult_theta_av = 0.5*(s_n[4:6]+s[4:6])
    end

    imult_x = imult[1:3]
    imult_theta = imult[4:6]
    mult_x = Vec3(mult[1:3])
    mult_theta = Vec3(mult[4:6])




    x1 = cf[inode1].x
    x2 = cf[inode2].x

    psi1  = cf[inode1].psi
    psi2  = cf[inode2].psi

  

    x =  cf[inode].x
    mult = pull_vectors(inv_loc_mult, Dy)
    mult_n = pull_vectors(inv_loc_mult, y_n)

    pos_imposed = nbc.init_positions[inode]  + disp  

    res_el[imult] = scale_factor*(x - pos_imposed).v
    res_el[ix] = scale_factor*mult
    S_el[ix, imult] = -scale_factor*eye(3)
    S_el[imult,ix] = S_el[ix, imult]'

 
#    println("S_el ", S_el)
#   println("inv_loc ", inv_loc)

    ext_work_el = 0.5*transpose(disp.v)*(mult + mult_n)

    if matrix == false 
        push_element_sparse(res,iel2,inv_loc,S_el,res_el)
        return ext_work_el
    else 
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc,S_el)
    end
  

    end

    function prismatic_joint(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},itime::Int,h::Float64)
        N = size(Dy,1)
        res = zeros(N)
        matrix = true
        println("prismatic_joint")
        return inv_loc, S_el = node_torque(nbr,Dy,y_n, res,itime,h,matrix)
    end
    
    function prismatic_joint(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},res::Vector{Float64},itime::Int,h::Float64)
        matrix = false
        prismatic_joint(nbr,Dy,y_n,res,itime,h,matrix)
    end
