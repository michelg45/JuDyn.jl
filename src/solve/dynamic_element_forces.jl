
"""
    dynamic_element_forces
"""
function dynamic_element_forces(Nel::Int64, element_numbers::Vector{Int},
         element_types::Vector{String},y_n::Vector{Float64},Dy::Vector{Float64},ydot_np1::Vector{Float64},
        res::Vector{Float64}, p::Vector{Float64},ext_work::Float64,itime::Int,h::Float64,niter::Int)
        
    
        # ext_work = ext_work_p
        kin_energy = 0.0
        pot_energy = 0.0
        str_energy = 0.0

        res .= 0.0
        p   .= 0.0

    for iel = 1:Nel

        nbr = element_numbers[iel]
        el_type = element_types[iel]

        if el_type == "rigid_body"
            kin_el, pot_el = rigid_body_force(nbr,y_n,Dy,ydot_np1,res,p,itime,h)
            kin_energy += kin_el
            pot_energy += pot_el
        elseif el_type == "beam"
            kin_el, str_el, pot_el = beam_force(nbr,y_n,Dy,ydot_np1,res,p,itime,h,niter)
            kin_energy += kin_el
            str_energy += str_el
            pot_energy += pot_el
        elseif el_type == "frame_link"
            frame_link_force(nbr,Dy,res)
        elseif el_type == "node_link"
            node_link_force(nbr,Dy,res)
        elseif el_type == "node_force"
            ext_work_el = node_force(nbr,Dy,res,p,itime,h)
            global ext_work += ext_work_el
        elseif el_type == "ground_hinge"
            str_el, ext_work_el = ground_hinge_force(nbr,Dy,y_n,ydot_np1,res,p,niter,itime,h)
            global ext_work += ext_work_el
            str_energy += str_el
        elseif el_type == "hinge"
            str_el, ext_work_el = hinge_force(nbr,Dy,y_n,ydot_np1,res,p,niter,itime,h)
            global ext_work += ext_work_el
            str_energy += str_el
        elseif el_type == "shell"
            kin_el, str_el, pot_el = shell_force(nbr,y_n,Dy,ydot_np1,res,p,itime,niter,h)
            str_energy += str_el
            pot_energy += pot_el
            kin_energy += kin_el
        elseif el_type == "rigid_mass"
            kin_el, pot_el = rigid_mass_force(nbr,y_n,Dy,ydot_np1,res,p,itime,h)
            kin_energy += kin_el
            pot_energy += pot_el 
        elseif el_type == "ground_spring_damper"
            str_el, ext_work_el = ground_spring_damper_force(nbr,Dy,ydot_np1,res)
            global ext_work += ext_work_el
            str_energy += str_el        
        elseif el_type == "superelement"
                (pot_el,kin_el, str_el) = super_element_force(nbr,y_n,Dy,ydot_np1,res,p)
                kin_energy += kin_el
                str_energy += str_el
                pot_energy += pot_el 
        elseif el_type == "node_displacement"
                ext_work_el = node_displacement(nbr,Dy,y_n,res,itime,h)
                global ext_work += ext_work_el 
        elseif el_type == "lin_constr"
                linear_constraint_force(nbr,y_n,Dy,res)
        elseif el_type == "super_beam"
                (pot_el,kin_el, str_el) = super_beam_force(nbr,y_n,Dy,ydot_np1,res,p)
                kin_energy += kin_el
                str_energy += str_el
                pot_energy += pot_el
        elseif el_type == "prismatic_joint"    
                str_el, ext_work_el = prismatic_joint_force(nbr,Dy,y_n,ydot_np1,res,p,itime,h)
                ext_work += ext_work_el
                str_energy += str_el
                       
"""        elseif el_type == "ground_spherical_joint"   
            ground_spherical_joint_force(nbr,Dy,res)
        elseif el_type == "spherical_joint"
            spherical_joint_force(nbr,Dy,res,)"""
        """




        elseif el_type == "frame_spring"
            str_el = frame_spring(nbr,res)
            str_energy += str_el



        elseif el_type == "node_torque"
            ext_work_el = node_torque(nbr,Dy,res,p,itime,h)
            ext_work += ext_work_el
       
        """
        end

    end

    # return  kin_energy, pot_energy, str_energy, ext_work
    return  kin_energy, pot_energy, str_energy
end