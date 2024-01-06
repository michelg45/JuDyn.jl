"""
    eigenvalue_system_dynamic
"""
function  eigenvalue_system_dynamic(Nel::Int64, element_numbers::Vector{Int64},element_types::Vector{String},Ndofs::Int64)    

    global A = spzeros(Ndofs,Ndofs)
    global B = spzeros(Ndofs,Ndofs)

    tol = eps(Float64)

    itime = 0
    h = 0.0
    niter = 1

    matrix = true

    theta_p = 0.0
    alpha_stiff = 1.0

    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]


        if el_type == "rigid_body"
            inv_loc, S_el =  rigid_body(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "beam"
            inv_loc, S_el =  beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "visco_beam"
            inv_loc, S_el =  visco_beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h,niter)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "shell"
            inv_loc, S_el =  shell(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,niter,h)
            A[inv_loc,inv_loc] += S_el
        end
    end 
   
    droptol!(A,tol)

    theta_p = 1.0
    alpha_stiff = 0.0


    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]
    
        if el_type == "rigid_body"
            inv_loc, S_el =   rigid_body(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "beam"
            inv_loc, S_el =   beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "visco_beam"
            inv_loc, S_el =   visco_beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h,niter)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "shell"
            inv_loc, S_el =   shell(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,niter,h)
            B[inv_loc,inv_loc] += S_el
        end
    end

    droptol!(B,tol)

    return A, B

end