"""
    ground_spring_damper
"""
function ground_spring_damper(nbr::Int,Dy::Vector{Float64},ydot::Vector{Float64},ydot_n::Vector{Float64},res::Vector{Float64}, alpha_stiff::Float64, theta_p::Float64,matrix::Bool)


    gsdc = Main.SetElements.ground_spring_damper_container
    ec =  Main.element_container
    cf =  Main.Frames.current_frames
    mc =  Main.model_container


    iel = findfirst(x -> x == nbr, gsdc.number)[1]
    inode = gsdc.node_order[iel]
    params = gsdc.params[iel]
    x_0 = gsdc.position[iel].v

    rotation = mc.uniform_rotation
    Omega_d = mc.rotation_speed
    vit_rot = norm2(Omega_d)
    (rotation == true && vit_rot > 0.0) && (tOmega_d = tilde(Omega_d).mat ; 
    k = (Omega_d/vit_rot).v; dir = ones(3) - dot(ones(3),k)*k) 


    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]

    inv_loc = ec.inv_loc_x[iel2]

    x = cf[inode].x
   
    Dx, xdot, xdot_n = pull_vectors(inv_loc, Dy, ydot, ydot_n)


    elong = x.v + Dx -x_0
 
    if (rotation == true && vit_rot > 0.0) 
         f_elast = params[1]*dir.*elong 
         f_dissip = params[4]*dir.*(xdot + tOmega_d*elong)
         S_el = alpha_stiff*params[1]*diagm(dir) + params[4]*diagm(dir)*(theta_p*eye(3)+ alpha_stiff*tOmega_d)
    else
         f_elast = params[1:3].*elong
         f_dissip = params[4:6].*xdot
         S_el = alpha_stiff*diagm(params[1:3]) + theta_p*diagm(params[4:6])
    end
  
    res_el = - f_elast  - f_dissip

    if matrix == false 
     push_element_sparse(res,iel2,inv_loc,S_el,res_el)

     str_el = 0.5*transpose(elong)*f_elast
     ext_work_el = - 0.5*transpose(Dx)*(f_dissip+params[4:6].*xdot_n)

     return str_el, ext_work_el
    else
     return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc,S_el)
    end

end

"""
    ground_spring_damper_force
"""
function ground_spring_damper_force(nbr::Int,Dy::Vector{Float64},ydot::Vector{Float64},res::Vector{Float64})


     gsdc = Main.SetElements.ground_spring_damper_container
     ec =  Main.element_container
     cf =  Main.Frames.current_frames
     mc =  Main.model_container
 
 
     iel = findfirst(x -> x == nbr, gsdc.number)[1]
     inode = gsdc.node_order[iel]
     params = gsdc.params[iel]
     x_0 = gsdc.position[iel].v
 
     rotation = mc.uniform_rotation
     Omega_d = mc.rotation_speed
     vit_rot = norm2(Omega_d)
     (rotation == true && vit_rot > 0.0) && (tOmega_d = tilde(Omega_d).mat ; 
     k = (Omega_d/vit_rot).v; dir = ones(3) - dot(ones(3),k)*k) 
 
 
     iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]
 
     inv_loc = ec.inv_loc_x[iel2]
 
     x = cf[inode].x
 
     Dx, xdot = pull_vectors(inv_loc, Dy, ydot)
 
 
     elong = x.v + Dx - x_0
  
     if (rotation == true && vit_rot > 0.0) 
          f_elast = params[1]*dir.*elong 
          f_dissip = params[4]*dir.*(xdot + tOmega_d*elong)
     else
          f_elast = params[1:3].*elong
          f_dissip = params[4:6].*xdot
     end
   
     res_el = - f_elast  - f_dissip

     str_el = 0.5*transpose(elong)*f_elast
 
     for i = 1:3
          iloc = inv_loc[i] 
          iloc > 0  &&  (res[iloc] = +res_el[i])
     end
 
      return str_el, ext_work_el
 
 end