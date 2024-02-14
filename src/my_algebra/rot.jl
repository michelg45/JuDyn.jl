"""
    rot(psi) 

Function computing the rotation operator `R::Mat3` from the Cartesian rotation vector `psi::RV3`

````math        
\\mathbf{R} = rot(\\mathbf{\\psi})
````
Calling sequence: 

````{verbatim}
R = rot(psi)
````

"""
    function rot(psi::RV3)
  
        I3 = Matrix(1.0I,3,3)
        psi = psi.v
        PSI2 = sum(@.psi^2)
        PSI = sqrt(PSI2)
        tps = [0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0]
        if PSI <= PREC  
            return Mat3(I3+tps)
        else
            return Mat3(I3*cos(PSI)+sin(PSI)/PSI*tps + ((1.0-cos(PSI))/PSI2*psi)*psi')
        end
    
    end
    
"""
    rot(psi,a)

Function computing the rotation of a vector `a::Vec3` by a Cartesian vector `psi::RV3`

````math        
\\mathbf{b} = rot(\\mathbf{\\psi},\\mathbf{a})
````

Calling sequence: 

````{verbatim}
b = rot(psi,a)
````

"""
function rot(psi::RV3,a::Vec3)

    psi = psi.v
    a = a.v
    PSI2 = sum(@.psi^2)
    PSI = sqrt(PSI2)

    if PSI > PREC
        return Vec3(cos(PSI)*a+(sin(PSI)/PSI)*cross(psi,a)+(((1.0-cos(PSI))/PSI2)*dot(a,psi))*psi)
    else
       return Vec3(a+cross(psi,a))
    end
        
    end


"""
    rot(phi,axe)

Function computing the rotation operator `R::Mat3` resulting from a rotation angle         `phi::Float64 about a coordinate axis  `axe::Int` (1,2 or 3).

````math        
\\mathbf{R} = rot(\\phi,\\mathbf{n})
````
Calling sequence: 

````{verbatim}
R = rot(phi,axe)
````

"""
function rot(phi::Float64,axe::Int)

    n=zeros(3)
    n[axe]=1.0
    R=Mat3(diagm(cos(phi)*(ones(3)-n)+n))+sin(phi)*tilde(Vec3(n))

    return R

end

"""
    rot(p)

Function computing the rotation operator  `R::Mat3` from the Euler-Rodrigues parameters  `p::RP3`

````math        
\\mathbf{R} = rot(\\mathbf{p})
````
Calling sequence: 

````{verbatim}
R = rot(p)
````

"""
function rot(p::RP3)

    I3 = Matrix(1.0I,3,3)
    P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
    tps = [0.0 -p[3] p[2]; p[3] 0.0 -p[1]; -p[2] p[1] 0.0]

    return R = Mat3(I3 + (P_0*eye(3) + 0.5*tps)*tps)

        end

        function rot2(psi::RV3)

            """
                function rot2(psi)

                input: rotation vector psi::RV3
                output: rotation operator R
                PREC: precision set to machine precision eps**1/2.
                if (PSI < PREC)  R is limited to its linearized expression.

            """
                PSI2=psi[1]*psi[1]+psi[2]*psi[2]+psi[3]*psi[3]
                tps=Mat3([0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0])
                addb = 0.5
                adda = 1.0
                beta = addb
                alpha = adda
                for i=1:5
                    addb = -addb*PSI2/((2*i+1)*(2*i+2))
                    adda = -adda*PSI2/((2*i)*(2*i+1))
                    beta += addb
                    alpha += adda
                end
                    R= I3+(alpha*I3+beta*tps)*tps

                return R
            end
