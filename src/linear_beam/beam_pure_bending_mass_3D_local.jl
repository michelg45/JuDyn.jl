"""
    beam_pure_bending_mass_3D_local

    Function computing the mass matrix of a pure bending (i.e. omission of extension and torsion) 3D beam element in local coordinates.  The degrees of freedom of the element are thus: 
            (0 u_y1 u_z1 0 phi_y1 phi_z1 0 u_y2 u_z2 0 phi_y2 phi_z2)

        Calling sequence:
            M = beam_pure_bending_mass_3D(l,mass_properties)

        input:

            2 cases:
            
            if size(mass_properties) == 2
            
                data for beam without rotary inertia of cross sections in bending
            
                m = mass_properties[1] (mass per unit length)
                Jxx = mass_properties[2] (rotary inertia in torsion per unit length)
            
            elseif size(mass_properties) == 4
            
                data for beam without rotary inertia of cross sections in bending
            
                m = mass_properties[1] (mass per unit length)
                Jyy = mass_properties[2] (rotary inertia about Y  axis per unit length)
                Jzz = mass_properties[2] (rotary inertia about Z axis per unit length)
                end
            
        Output:
            
              M::Array{Float64,2}   (12 x 12) mass matrix of beam element

"""
function  beam_pure_bending_mass_3D_local(l::Float64,mass_properties::Vector{Float64})

        if size(mass_properties,1) == 2
            m = mass_properties[1]
            Jyy = 0.
            Jzz = 0.
        elseif size(mass_properties,1) == 4
            m = mass_properties[1]
            Jyy = mass_properties[3]
            Jzz =mass_properties[4]
        else
            error("wrong size of mass properties vector")
        end


        M=zeros(12,12)

        M_bend_y=beam_mass_2D(l,m,Jzz)
        loc_bend_y=[2 6 8 12]
        M[loc_bend_y,loc_bend_y]=M_bend_y
        M_bend_z=beam_mass_2D(l,m,Jyy)
        loc_bend_z=[3 5 9 11]
        One=eye(4)
        One[2,2]=-1
        One[4,4]=-1
        M_bend_z=One*M_bend_z*One
        M[loc_bend_z,loc_bend_z]=M_bend_z

        return M
    end