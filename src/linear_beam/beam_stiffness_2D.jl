"""
    beam_stiffness_2D

        function constructing the  stiffness matrix of a 2D beam element in local coordinates.

            calling sequences: 

                M = beam_stiffness_2D(l,EI,phi)

            input:

                l::Float64    element length
                m::Float64    stiffness per unit length
                phi::Float64  shear parameter phi =   12*EI/(kAG*l^2)

            Output:

                K::Array{Float64,2}   (4x4)  stiffness matrix of bar element
"""
function beam_stiffness_2D(l::Float64, EI::Float64, phi::Float64)

    K=zeros(4,4)
    K[1,1]= 12.0*EI/l^3
    K[2,2]= EI/l*(4.0+phi)
    K[3,3]=K[1,1]
    K[1,3]=-K[1,1]
    K[3,1]=-K[1,1]
    K[4,4]=K[2,2]
    K[1,2]= 6.0*EI/l^2
    K[2,1]=K[1,2]
    K[1,4]=K[1,2]
    K[4,1]=K[1,4]
    K[2,3]= -K[1,2]
    K[3,2]= K[2,3]
    K[2,4]= EI/l*(2.0-phi)
    K[4,2]=K[2,4]
    K[4,3]=-K[1,2]
    K[3,4]=K[4,3]
    K[:,:]=K[:,:]/(1+phi)
    return K

end