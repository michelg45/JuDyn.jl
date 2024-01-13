"""
    linear_beam_element

            function contructing the stiffness and mass matrices of a linear beam element in global or local coordinates.

                calling sequences: 
                    K,M = linear_beam_element(x_1,x_2,x_3,stiffness_properties,mass_properties)
                    K,M = linear_beam_element(length,stiffness_properties,mass_properties)
    

    
        input :
    
            nel::Int                                element number
            x_1, x_2, x_3 ::Vec3                    set of 3 Nnodes defining the beam axis and the beam reference plane.
            length::Float64                         beam length (for element in local coordinates)
            stiffness_properties::Vector{Float64}   set of beam stiffness properties 
            mass_properties::Vector{Float64}        set of beam stiffness properties

            note: stiffness_properties and mass_properties vectors are described in 'beam_stiffness_3D_local' and 'beam_mass_3D_local'.
    
    
        output:
    
            K::Array{Float64,2}               linear stiffness matrix in global axes
            M::Array{Float64,2}               linear mass matrix in global axes
"""
    function  linear_beam_element(x_1::Vec3,x_2::Vec3,x_3::Vec3,stiffness_properties::Vector{Float64},mass_properties::Vector{Float64})


    
    # Definition of beam geometry and construction of rotation operator R
    
        tol=sqrt(eps(Float64))
        l=norm2(x_2-x_1)
        d=norm2(x_3-x_1)
        if l < tol
            error("element ",nel, " has zero length")
        end
        n_1=(x_2-x_1)/l
        n_2=(x_3-x_1)/d
        sina=sin(acos(dotp(n_1,n_2)))
        if sina < tol
            error("element ",nel, " has colinear nodes")
        end
        n_3=crossp(n_1,n_2)/sina
        n_2=crossp(n_3,n_1)
        R=[n_1.v n_2.v n_3.v]
    
    
     # massic properties of cross section
    
    
    
    # construction of stiffness and mass matrices in local axes
    
        K=beam_stiffness_3D_local(l,stiffness_properties)
        M=beam_mass_3D_local(l,mass_properties)
    
    # rotation to global axes
    
        for k=1:4
            kin=(k-1)*3
            index_k=[kin+1, kin+2, kin+3]
                for j=1:4
                jin=(j-1)*3
                index_j=[jin+1, jin+2, jin+3]
                K[index_k,index_j]=R*K[index_k,index_j]*R'
                M[index_k,index_j]=R*M[index_k,index_j]*R'
            end
        end
    
        return K, M
    
    end
    
    function  linear_beam_element(length::Float64,stiffness_properties::Vector{Float64},mass_properties::Vector{Float64})
    
    # construction of stiffness and mass matrices in local axes
    
        K=beam_stiffness_3D_local(length,stiffness_properties)
        M=beam_mass_3D_local(length,mass_properties)
    
        return K, M
    
    end