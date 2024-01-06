"""
    tilde(psi::Vec3)

    function constructing the skew-symmetric matrix (type ::Mat3) associated to a rotation vector (type ::Vec3).

       calling sequence: A = tilde(psi)
       
       example:

       psi = Vec3(5.0, 2.0, -1.0)
       A = tilde(psi)
       Mat3([0.000000e+00 1.000000e+00 2.000000e+00; -1.000000e+00 0.000000e+00 -5.000000e+00; 
       -2.000000e+00 5.000000e+00 0.000000e+00])

"""
function tilde(psi::Vector{Float64})

        tpsi=zeros(3,3)
        tpsi[1,2]=-psi[3]
        tpsi[2,1]=psi[3]
        tpsi[1,3]=psi[2]
        tpsi[3,1]=-psi[2]
        tpsi[2,3]=-psi[1]
        tpsi[3,2]=psi[1]
    
        return tpsi
    
    end