"""
    linear_beam_element_pure_bending

      
        generation of a beam element in global or local axes in whixh extension and torsional terms are omitted. The function is used exclusively by 'super_beam_matrix_kernel.jl'  to construct the bending contribution to the super_beam element.

            Calling sequences:
    
                linear_beam_element_pure_bending(x_1,x_2,x_3,stiffness_properties,mass_properties)
                linear_beam_element_pure_bending(length,stiffness_properties,mass_properties)

        input :

        nel::Int                               element number
        length::Float64                       element length
        x_1, x_2, x_3::Vec3                   node coordinates defining beam axis and beam reference plane.
        stiffness properties::Vector{Float64}    [0 kAG_y kAG_z 0 EI_yy EI_zz] or [0 0 EI_yy EI_zz]
        mass properties::Vector{Float64}         [m 0 Jyy Jzz] or [m 0]

        output:

            K::Array{Float64,2}      linear stiffness matrix in global / local  axes
            M::Array{Float64,2}      linear mass matrix in global /local  axes
"""
function  linear_beam_element_pure_bending(x_1::Vec3,x_2::Vec3,x_3::Vec3,stiffness_properties::Vector{Float64},mass_properties::Vector{Float64})
    
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
   
    # construction of stiffness and mass matrices in local axes
    
        K=beam_pure_bending_stiffness_3D_local(l,stiffness_properties)
        M=beam_pure_bending_mass_3D_local(l,mass_properties)
    
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
    
    function  linear_beam_element_pure_bending(length::Float64,stiffness_properties::Vector{Float64},mass_properties::Vector{Float64})
 
    # construction of stiffness and mass matrices in local axes
    
        K=beam_pure_bending_stiffness_3D_local(length,stiffness_properties)
        M=beam_pure_bending_mass_3D_local(length,mass_properties)
    
        return K, M
    
    end