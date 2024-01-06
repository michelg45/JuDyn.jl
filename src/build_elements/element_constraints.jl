"""
    element_constraints 
"""
function element_constraints(Ndim::Int64 ,Nel::Int64, Nbounds::Int64,bounds::Vector{Float64},y_n::Vector{Float64},Dy::Vector{Float64},itime::Int,h::Float64)

    ec =  Main.element_container
    inc = Main.SetElements.inequality_container
    cf =  Main.Frames.current_frames
    input_functions = Main.InputFunctions.input_functions

    global B = zeros(Ndim,Nbounds)
    global g_n = zeros(Nbounds)


    ibound = 0
    for iel = 1:Nel
        el_type = ec.element_types[iel]
        if el_type == "inequality" 
            ibound +=1
            nbr = ec.element_numbers[iel]
            loc_x = ec.inv_loc_x[iel]
            iel2 = findfirst(x -> x == nbr, inc.number)[1]
            time_function = input_functions[inc.time_function[iel2]]
            bounds[ibound] = time_function(itime,h,inc.bound_params[iel2])
            type = inc.type[iel2]
            if type == "linear"
                B[loc_x,ibound]= inc.direction[iel2].v
                g_n[ibound] = B[loc_x,ibound]'*(y_n+Dy)[loc_x]
            elseif type == "point_along_direction"
                B[loc_x,ibound]= inc.direction[iel2].v
                g_n[ibound] = B[loc_x,ibound]'*(y_n+Dy)[loc_x]
            elseif type == "point_to_axis"
                p1 = inc.points[iel2][1]
                p2 = inc.points[iel2][2]
                x= Vec3((y_n)[loc_x])
                # x= cf[inc.node[iel2]].x
                f, dfdx = point_to_axis(x,p1,p2)
                B[loc_x,ibound]= dfdx
                x= Vec3((y_n+Dy)[loc_x])
                f, dfdx = point_to_axis(x,p1,p2)
                # B[loc_x,ibound]= dfdx
                g_n[ibound] = f
            elseif type == "distance_to_point"
                p1 = inc.points[iel2][1]
                x= Vec3((y_n)[loc_x])
                f, dfdx = distance_to_point(x,p1)
                B[loc_x,ibound]= dfdx
                x= Vec3((y_n+Dy)[loc_x])
                f, dfdx = distance_to_point(x,p1)
                # B[loc_x,ibound]= dfdx
                g_n[ibound] = f
                
                
            end
        end
    end

    return B, g_n,  bounds

end
