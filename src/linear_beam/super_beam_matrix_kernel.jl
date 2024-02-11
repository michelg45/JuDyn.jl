"""
    super_beam_matrix_kernel

    Construction of the linear matrix kernel of for a "super_beam"  model  using Herting's method. The model results from the assembly of 2 linear beam models based on cubic interpolation of the bending deflection and linear interpolation of the shear force.

    Calling sequence: 
        mass, Jrot, K_elast, M_elast = super_beam_matrix_kernel(length, stiffness_properties, mass_properties)


    Input:
        length::Float64		  			        length of the beam element
        stiffness_properties::Vector{Float64}	stiffness parameters as input to  function "linear_beam_element".
        mass_properties::Vector{Float64}		mass parameters as  input to  function "linear_beam_element".

    Output:
        M_rig::Array{Float64,2}					6x6 inertia tensor of rigid body
        K_elast::Array{Float64,2}				stiffness matrix asssociated to pure elastic deformation
        M_elast::Array{Float64,2}				associated mass matrix
        S::Array{Float64,2}                     associated gyroscopic matrix.
"""
    function super_beam_matrix_kernel(length::Float64,stiffness_properties::Vector{Float64}, mass_properties::Vector{Float64})

        l = length / 2.0
        
        K_el, M_el = linear_beam_element_pure_bending(l, stiffness_properties, mass_properties)
        l_1 = [i for i in 1:6]
        l_2 = [i for i in 7:12]
        l_3 = [i for i in 13:18]
        locel = Vector{Vector{Int}}(undef,2)
        locel[1] = [l_1; l_3]
        locel[2] = [l_3; l_2]
        
        K = Array{Float64,2}(zeros(18,18))
        M = Array{Float64,2}(zeros(18,18))
        
        for nel = 1:2
            K[locel[nel], locel[nel]] += K_el
            M[locel[nel], locel[nel]] += M_el
        end

        l_ext = [1; 7; 13]
        l_tors = [4; 10; 16]

        EA = stiffness_properties[1]
        m = mass_properties[1]

        degree = "quadratic"

        K_el = bar_stiffness_2D(length,EA,degree)
        M_el = bar_mass_2D(length,m,degree)

        K[l_ext,l_ext] += K_el
        M[l_ext,l_ext] += M_el

        if size(stiffness_properties,1) == 6
			GJ = stiffness_properties[4]
		elseif size(stiffness_properties,1) == 4
			GJ = stiffness_properties[2]
		end
        J_xx =  mass_properties[2]

        K_el = bar_stiffness_2D(length,GJ,degree)
        M_el = bar_mass_2D(length,J_xx,degree)

        K[l_tors,l_tors] += K_el
        M[l_tors,l_tors] += M_el  
        
        val, modes = eigen(K,M)
        
        U = zeros(18,6)
        
        
        A = zeros(6,6)
        for i = 1:6
         A[i,i] = 1.0
        end
        A[2,6] = -length/2.0
        A[3,5] = +length/2.0
        
        U= modes[:,1:6]*(modes[1:6,1:6]\A)
        
        mass = transpose(U[:,1:3]) * M * U[:,1:3]
        
        Jrot = transpose(U[:,4:6]) * M * U[:,4:6]
        
        K_extended = [K M*U; transpose(U)*M zeros(6,6)]
        
        I_12 = zeros(12,12)
        
        for i = 1:12
            I_12[i,i] = 1.0
        end
        
        p = [I_12; zeros(12,12)]
        
        GG = K_extended\p
        
        T = M*U
        
        M_rig = U'*T
        
        GG[1:18,:] -= U*(M_rig\(transpose(T)*GG[1:18,:]))
        
        
        F_BB = GG[1:12,1:12]
        F_IB = GG[13:18,1:12]
        
        R = [I_12;  F_IB*(F_BB\I_12)]
        
        K_elast =  transpose(R)*K*R
        M_elast = transpose(R)*M*R

        mass = diag(M_rig)[1]
        Jrot = diag(M_rig)[4:6]

        # Construction of the gyroscopic contribution S such that 
        #           G(omega) = - Σ  omega[i] * S[i]
        #           s(q,r)[i] = [q'*S[i]*r]


        S = beam_gyr_3D_local(l,m,degree)
        for i = 1:3
            S[i] = R'*S[i]*R
        end

        return mass, Jrot, K_elast, M_elast, S

        end
        

