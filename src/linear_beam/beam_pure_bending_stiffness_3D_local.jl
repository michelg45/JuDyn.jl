"""
    beam_pure_bending_stiffness_3D_local

    Function computing the mass matrix of a pure bending (i.e. omission of extension and torsion) 3D beam element in local coordinates.  The degrees of freedom of the element are thus: 
            (0 u_y1 u_z1 0 phi_y1 phi_z1 0 u_y2 u_z2 0 phi_y2 phi_z2)

        Calling sequence:
            K = beam_pure_bending_stiffness_3D(l,stiffness_properties)

        input:

		if size(stiffness_properties) == 6

			data for Timoshenko beam

			
			kAG_y = stiffness_properties[2] (shear stiffness  in y direction)
			kAG_z = stiffness_properties[3] (shear stiffness  in z direction)
	
			EI_yy = stiffness_properties[5] (bending stiffness about y axis)
			EI_zz = stiffness_properties[6] (bending stiffness about z axis)

		elseif size(stiffness_properties) == 4

			data for Bernoulli beam

			EI_yy = stiffness_properties[3] (bending stiffness about y axis)
			EI_zz = stiffness_properties[4] (bending stiffness about z axis)

		end
            
        Output:
            
              K::Array{Float64,2}   (12 x 12) stiffness matrix of beam element

"""
function beam_pure_bending_stiffness_3D_local(l::Float64,stiffness_properties::Vector{Float64})

		if size(stiffness_properties,1) == 6
			
			kAG_y = stiffness_properties[2]
			kAG_z = stiffness_properties[3]
			EI_yy = stiffness_properties[5]
			EI_zz = stiffness_properties[6]

			phi_z =   12*EI_yy/(kAG_z*l^2)
			phi_y =   12*EI_zz/(kAG_y*l^2)

		elseif size(stiffness_properties,1) == 4

			EI_yy = stiffness_properties[3]
			EI_zz = stiffness_properties[4]

			phi_z = 0.
			phi_y = 0.
		else
			error("wrong size of stiffness properties vector")
		end


        K=zeros(12,12)

        K_bend_y=beam_stiffness_2D(l,EI_zz,phi_y)
        loc_bend_y=[2 6 8 12]
        K[loc_bend_y,loc_bend_y]=K_bend_y

        K_bend_z=beam_stiffness_2D(l,EI_yy,phi_z)

        One=eye(4)
        One[2,2] = -1.0
        One[4,4] = -1.0

        K_bend_z = One*K_bend_z*One
        loc_bend_z=[3 5 9 11]
        K[loc_bend_z,loc_bend_z]=K_bend_z

		return K

    end
