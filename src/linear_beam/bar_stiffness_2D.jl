"""
    bar_stiffness_2D

        function constructing the linear or quadratic stiffness matrix of a 2D bar element in local coordinates.

            calling sequences: 

                M = bar_stiffness_2D(l,EA)
                M = bar_stiffness_2D(l,EA,"quadratic")

            input:

                l::Float64    element length
                EA::Float64   stiffness per unit length

            Output:

                K::Array{Float64,2}   (2x2) or (3x3) stiffness matrix of bar element
"""
function  bar_stiffness_2D(l::Float64,EA::Float64)

    K=zeros(2,2)
    K[1,1]= EA/l
    K[2,2]=K[1,1]
    K[1,2]= -K[1,1]
    K[2,1]=K[1,2]

    return K
end

function  bar_stiffness_2D(l::Float64,EA::Float64,degree::String)

  degree != "quadratic" && (error("bar_mass_2D:  degree must be quadratic "))
  
    K=zeros(3,3)
    K[1,1]= 7.0*EA/(3.0*l)
    K[2,2]=K[1,1]
    K[3,3]= 16.0*EA/(3.0*l)
    K[1,2]= EA/(3.0*l)
    K[2,1]=K[1,2]
    K[1,3]= -8.0*EA/(3.0*l)
    K[2,3]= K[1,3]
    K[3,1]= K[1,3]
    K[3,2]= K[2,3]

    return K
end
