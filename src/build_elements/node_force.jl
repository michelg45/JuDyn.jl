
"""
    node_force
"""
function node_force(nbr::Int,Dy::Vector,res::Vector,p::Vector,itime::Int,h::Float64)


    fc = SetElements.force_container
    nbc = Main.node_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    input_functions = Main.InputFunctions.input_functions


    iel = findfirst(x -> x == nbr, fc.number)[1]

    str_time_function = fc.time_function[iel]

    time_function = input_functions[str_time_function]
    axe = (fc.axe[iel]).v
    params = fc.params[iel]

    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    
    inv_loc_x= ec.inv_loc_x[iel2]

    amp = time_function(itime,h,params)

    p[inv_loc_x] += axe*amp

    res[inv_loc_x] += axe*amp

    ext_work_el = 0.5*(amp + time_function(itime-1,h,params))*transpose(Dy[inv_loc_x])*axe

    return ext_work_el

    end
