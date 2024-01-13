"""
    beam_mass_3D

        function computing the mass matrix of a 3D beam element in local coordinates.

            Degrees of freedom:
                (u_x1 u_y1 u_z1 phi_x1 phi_y1 phi_z1 u_x2 u_y2 u_z2 phi_x2 phi_y2 phi_z2)

            calling sequence: 
                M = beam_mass_3D(l,mass_properties)

            input:

                l::Float64                          element length
                mass_properties::Vector(Float64)	mass properties of the element

                2 cases:
            
                    if size(mass_properties) == 2
            
                        data for beam without rotary inertia of cross sections in bending
            
                        m::Float64 = mass_properties[1]     mass per unit length
                        Jxx::Float64 = mass_properties[2]   rotary inertia in torsion per unit length
            
                    elseif size(mass_properties) == 4
            
                        data for beam without rotary inertia of cross sections in bending
            
                        m::Float64 = mass_properties[1]     mass per unit length
                        Jxx::Float64 = mass_properties[2]   rotary inertia in torsion per unit length
                        Jyy::Float64 = mass_properties[2]   rotary inertia about Y  axis per unit length
                        Jzz::Float64 = mass_properties[2]   rotary inertia about Z axis per unit length
                    end
            
                Output:
            
                  M::Array{Float64,2}   (12 x 12) mass matrix of THE beam element

"""
function  beam_mass_3D_local(l::Float64,mass_properties::Vector{Float64})

        if size(mass_properties,1) == 2
            m = mass_properties[1]
            Jxx = mass_properties[2]
            Jyy = 0.
            Jzz = 0.
        elseif size(mass_properties,1) == 4
            m = mass_properties[1]
            Jxx = mass_properties[2]
            Jyy = mass_properties[3]
            Jzz =mass_properties[4]
        else
            error("beam_mass_3D_local - wrong size of mass properties vector")
        end

        M=zeros(12,12)
        M_ext=bar_mass_2D(l,m)
        loc_ext=[1 7]
        M[loc_ext,loc_ext]=M_ext
        M_tors=bar_mass_2D(l,Jxx)
        loc_tors=[4 10]
        M[loc_tors,loc_tors]=M_tors
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