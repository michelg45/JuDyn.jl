"""
    bar_mass_2D

        function constructing the linear or quadratic mass matrix of a 2D bar element in local coordinates.

            calling sequences: 

                M = bar_mass_2D(l,m)
                M = bar_mass_2D(l,m,"quadratic")

            input:

                l::Float64  element length
                m::Float64   mass per unit length

            Output:

                M::Array{Float64,2}   (2x2) or (3x3) mass matrix of bar element
"""
function  bar_mass_2D(l::Float64,m::Float64)



    M=zeros(2,2)
    M[1,1]= m*l/3.0
    M[2,2]=m*l/3.0
    M[1,2]= m*l/6.0
    M[2,1]=M[1,2]

    return M
end

function  bar_mass_2D(l::Float64,m::Float64,degree::String)

    if degree == "quadratic" 
        M=zeros(3,3)
        M[1,1]=2.0*m*l/15.0
        M[2,2]=M[1,1]
        M[3,3]=8.0*m*l/15.0
        M[1,2]= -m*l/30.
        M[2,1]=M[1,2]
        M[3,1]=m*l/15.0
        M[3,2]= M[3,1]
        M[1,3]= M[3,1]
        M[2,3]= M[3,2]
    else
        M = bar_mass_2D(l,m)
    end
    return M
end