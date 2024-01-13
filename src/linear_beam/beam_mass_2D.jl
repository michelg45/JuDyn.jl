"""
    beam_mass_2D

        function constructing the  mass matrix of a 2D beam element in local coordinates.

            calling sequences: 

                M = beam_mass_2D(l,m,J)

            input:

                l::Float64  element length
                m::Float64  mass per unit length
                J::Float64  rotary inertia per unit length

            Output:

                M::Array{Float64,2}   (4x4)  mass matrix of bar element
"""
function beam_mass_2D(l::Float64,m::Float64,J::Float64)

    M=zeros(4,4)
    M[1,1]= 156.0*m*l/420.0 + 6.0/5.0*J/l
    M[2,2]= 4.0*m*l^3/420.0 + 2.0/15.0*J*l
    M[3,3]=M[1,1]
    M[4,4]=M[2,2]
    M[1,2]= 22.0*m*l^2/420.0 + J/10.0
    M[2,1]=M[1,2]
    M[1,3]= 54.0*m*l/420.0 - 6.0/5.0*J/l
    M[3,1]=M[1,3]
    M[1,4]= -13.0*m*l^2/420.0 + J/10.0
    M[4,1]=M[1,4]
    M[2,3]= 13.0*m*l^2/420.0 - J/10.0
    M[3,2]= M[2,3]
    M[2,4]= -3.0*m*l^3/420.0-J/10.0*l
    M[4,2]=M[2,4]
    M[3,4]=-M[1,2]
    M[4,3]=M[3,4]

    return M
end