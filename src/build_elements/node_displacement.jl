
"""
    node_displacement
"""
function node_displacement(nbr::Int,Dy::Vector,y_n,res::Vector,itime::Int,h::Float64,matrix::Bool)


    dc = SetElements.disp_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, dc.number)[1]

    str_time_function = dc.time_function[iel]

    time_function = input_functions[str_time_function]
    axe = dc.axe[iel]
    params = dc.params[iel]
    scale_factor = dc.scale_factor[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inode = ec.element_nodes[iel2][1]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_mult= ec.inv_loc_mult[iel2]

    inv_loc = [inv_loc_x; inv_loc_mult]

    disp = time_function(itime,h,params)*axe

    S_el = zeros(6,6)

    res_el = zeros(6)

    ix = [1;2;3]
    imult = [4;5;6]

    inv_loc = [inv_loc_x; inv_loc_mult]

    Ds = pull_vectors(inv_loc, Dy)
    mult = Ds[imult]
    Dx = Vec3(Ds[ix])
    x =  cf[inode].x + Dx
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

    function node_displacement(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},itime::Int,h::Float64)
        N = size(Dy,1)
        res = zeros(N)
        matrix = true
        println("node_displacement")
        return inv_loc, S_el = node_torque(nbr,Dy,y_n, res,itime,h,matrix)
    end
    
    function node_displacement(nbr::Int,Dy::Vector{Float64},y_n::Vector{Float64},res::Vector{Float64},itime::Int,h::Float64)
        matrix = false
        node_displacement(nbr,Dy,y_n,res,itime,h,matrix)
    end
