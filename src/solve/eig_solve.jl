"""
    eig_solve
"""
function eig_solve(y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},max_vals::Int,itime::Int,h::Float64,niter::Int)

    element_container = Main.element_container
    model_container = Main.model_container
    element_numbers = element_container.element_numbers
    element_types = element_container.element_types
    Nel = model_container.Elements

    Ndofs_q = model_container.Ndofs_q
    Ndofs_x = model_container.Ndofs_x
    Ndofs_v = model_container.Ndofs_v
    Ndofs   = model_container.Ndofs

    global A = spzeros(Ndofs,Ndofs)
    global B = spzeros(Ndofs,Ndofs)

    tol = eps(Float64)
    matrix = true

    alpha_stiff = 1.0
    theta_p = 0.0

    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]


        if el_type == "beam"
            inv_loc, S_el = beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "shell"
            inv_loc, S_el = shell(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,niter,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "super_beam"
            inv_loc, S_el = super_beam(nbr,y,Dy,ydot,"stiffness",itime,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "rigid_body"
            inv_loc, S_el = rigid_body(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            A[inv_loc,inv_loc] += S_el
        elseif el_type == "superelement"
             inv_loc, S_el =superelement_herting(nbr,y,Dy,ydot,"stiffness")
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "ground_hinge"
             inv_loc, S_el = ground_hinge(nbr,Dy,y_n,res,p,niter,itime,h,matrix)
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "rigid_mass"
             inv_loc, S_el = rigid_mass(nbr,y,Dy,ydot,res,p,alpha_stiff,theta_p,matrix,itime,h)
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "node_link"
             inv_loc, S_el = node_link(nbr,Dy,res,matrix)
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "frame_link"
             inv_loc, S_el = frame_link(nbr,Dy,res,matrix)
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "hinge"   
             inv_loc, S_el = hinge(nbr,Dy,y_n,res,p,niter,itime,h,matrix)
             A[inv_loc,inv_loc] += S_el
        elseif el_type == "ground_spring_damper"
             inv_loc, S_el = ground_spring_damper(nbr,Dy,ydot_np1,ydot_n,res,alpha_stiff,theta_p,matrix)
             A[inv_loc,inv_loc] += S_el 
        end

    end 
   
    droptol!(A,tol)

    alpha_stiff = 0.0
    theta_p = 1.0

    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]
    
        if el_type == "beam"
            inv_loc, S_el = beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "shell"
            inv_loc, S_el = shell(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,niter,h)
            B[inv_loc,inv_loc] += S_el
        elseif  el_type == "super_beam"
            inv_loc, S_el = super_beam(nbr,y,Dy,ydot,"mass",itime,h)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "rigid_body"
            inv_loc, S_el = rigid_body(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "superelement"
            inv_loc, S_el = superelement_herting(nbr,y,Dy,ydot,"mass")
            B[inv_loc,inv_loc] += S_el
        elseif el_type == "rigid_mass"
             inv_loc, S_el = rigid_mass(nbr,y,Dy,ydot,res,p,alpha_stiff,theta_p,matrix,itime,h)
             B[inv_loc,inv_loc] += S_el
        elseif el_type == "ground_spring_damper"
             inv_loc, S_el = ground_spring_damper(nbr,Dy,ydot_np1,ydot_n,res,alpha_stiff,theta_p,matrix)
             B[inv_loc,inv_loc] += S_el
        end
        
    end


    vals  = eigs(B,A; nev = max_vals, tol = sqrt(eps(Float64)), which=:LM)
    eig_vals = [1/vals[1][i] for i in 1:max_vals]

    println("eigenvalues at step ",itime)

    println("real(values) ", real(eig_vals))
    println("imag(values) ", imag(eig_vals))


    return eig_vals
end