"""
    invrot(R::Mat3) 
    
        function calculating the rotation parameters of a (3x3) rotation operator.

        calling sequence: rv = invrot(R)

            input: rotation operator R::Mat3
            output: rotation vector psi::RV3
            
            for maximum accuracy the rotation vector is comuted from quaternion using
            the Cayley algorithm.
            PREC: precision set to machine precision sqrt(eps).
            if (PSI < PREC) psi is limited to the first-order expression: psi= 2*n*sin(PSI/2]
            
"""
function  invrot(R::Mat3)

    e=zeros(3)
    R = R.mat
    trR=R[1,1]+R[2,2]+R[3,3]
    G=R-transpose(R)
    S=R+transpose(R)
    cosP2=0.25*sqrt((1+trR)^2+G[3,2]^2+G[1,3]^2+G[2,1]^2)
    e[1]=0.25*sqrt(G[3,2]^2+ S[2,1]^2+S[3,1]^2+ (R[1,1]-R[2,2]-R[3,3]+1)^2)*sign(G[3,2])
    e[2]=0.25*sqrt(G[1,3]^2+ S[2,1]^2+S[3,2]^2+ (R[2,2]-R[1,1]-R[3,3]+1)^2)*sign(G[1,3])
    e[3]=0.25*sqrt(G[2,1]^2+ S[3,1]^2+S[3,2]^2+ (R[3,3]-R[2,2]-R[1,1]+1)^2)*sign(G[2,1])
    sinP2=norm(e)
    PSI=2*atan(sinP2,cosP2)
    PSI > PREC ? psi=PSI/sinP2*e : psi=2*e

    return  RV3(psi)
end

function RV3(R::Mat3)
    return invrot(R)
end
