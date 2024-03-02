"""
    static_element_system
"""
function static_element_system(Nel::Int64, element_numbers::Vector{Int},element_types::Vector{String},y_n::Vector{Float64},Dy::Vector{Float64},ydot_np1::Vector{Float64},res::Vector{Float64},       p::Vector{Float64},ext_work::Float64,itime::Int,h::Float64,niter::Int, theta_p)
    
        pot_energy = 0.0
        str_energy = 0.0

        res .= 0.0
        p   .= 0.0

        matrix = false
        
        alpha_stiff = 1.0

    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]

        if el_type == "rigid_body"
            kin_el, pot_el = rigid_body(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            pot_energy += pot_el
        elseif el_type == "beam"
            kin_el, str_el, pot_el = beam(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h,niter)
            str_energy += str_el
            pot_energy += pot_el
        elseif el_type == "frame_link"
            frame_link(nbr,Dy,res,matrix)
        elseif el_type == "spherical_joint"
            spherical_joint(nbr,Dy,res,matrix)
        elseif el_type == "node_link"
            node_link(nbr,Dy,res,matrix)
        elseif el_type == "node_force"
            ext_work_el = node_force(nbr,Dy,res,p,itime,h)
            global ext_work += ext_work_el
        elseif el_type == "node_torque"
            ext_work_el = node_torque(nbr,Dy,res,p,itime,h)
            ext_work += ext_work_el
        elseif el_type == "shell"
            kin_el, str_el, pot_el = shell(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,niter,h)
            str_energy += str_el
            pot_energy += pot_el
        elseif el_type == "ground_spring_damper"
            str_el, ext_work_el = ground_spring_damper(nbr,Dy,ydot_np1,ydot_n,res,alpha_stiff,theta_p,matrix)
            ext_work += ext_work_el
            str_energy += str_el 
        elseif el_type == "rigid_mass"
            kin_el, pot_el = rigid_mass(nbr,y_n,Dy,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            kin_energy += kin_el
            pot_energy += pot_el
        elseif el_type == "ground_hinge"
            ext_work_el = ground_hinge(nbr,Dy,y_n,res,p,niter,itime,h,matrix)
            ext_work += ext_work_el
        elseif el_type == "hinge"
            ext_work_el = hinge(nbr,Dy,y_n,res,p,niter,itime,h,matrix)
            ext_work += ext_work_el  
        elseif el_type == "superelement"
            (pot_el,kin_el, str_el) = super_element(nbr,y_n,Dy,ydot_np1,res,p,theta_p,alpha_stiff,matrix)
            kin_energy += kin_el
            str_energy += str_el
            pot_energy += pot_el 
        elseif el_type == "node_displacement"
            ext_work_el = node_displacement(nbr,Dy,y_n,res,itime,h)
            ext_work += ext_work_el 
        elseif el_type == "lin_constr"
            linear_constraint(nbr,y_n,Dy,res,matrix)
        elseif el_type == "super_beam"
            (pot_el,kin_el, str_el) = super_beam(nbr,y_n,Dy,ydot_np1,res,p,theta_p,itime,h)
            str_energy += str_el
            pot_energy += pot_el 
        elseif el_type == "prismatic_joint"    
            str_el, ext_work_el = prismatic_joint(nbr,Dy,y_n,ydot_np1,res,p,alpha_stiff,theta_p,matrix,itime,h)
            ext_work += ext_work_el
            str_energy += str_el
            
"""        elseif el_type == "ground_spherical_joint"
            ground_spherical_joint(nbr,Dy,res,alpha_stiff,matrix)
        elseif el_type == "spherical_joint"
            spherical_joint(nbr,Dy,res,alpha_stiff,matrix)"""

        """




        elseif el_type == "frame_spring"
            str_el = frame_spring(nbr,res)
            str_energy += str_el






"""
        end

    end

    return  pot_energy, str_energy
end