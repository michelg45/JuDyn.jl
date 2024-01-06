"""
    rot(psi::RV3) 
        function computing the rotation operator R::Mat3
        from the Cartesian rotation vector psi::RV3 

"""
    function rot(psi::RV3)
  
    psi = psi.v
    PSI=norm(psi)
    tps = [0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0]
    if PSI <= PREC  
        return Mat3(eye(3)+tps)
    else
        return Mat3(eye(3)*cos(PSI)+sin(PSI)/PSI*tps + ((1.0-cos(PSI))/PSI^2*psi)*psi')
    end
    
    end
    
"""
    rot(psi::RV3,a::Vec3)
        function computing the rotation of a vector a::Vec3 
        by a Cartesian vector psi::RV3

"""
    function rot(psi::RV3,a::Vec3)

        psi = psi.v
        a = a.v
        PSI = norm(psi)
        if PSI > PREC
            return Vec3(cos(PSI)*a+(sin(PSI)/PSI)*cross(psi,a)+(((1.0-cos(PSI))/PSI^2)*dot(a,psi))*psi)
        else
           return Vec3(a+cross(psi,a))
        end
        
    end


"""
    rot(phi::Float64,axe::Int)
        function computing the rotation operator R::Mat3 resulting from a rotation angle 
        phi::Float64 about a coordinate axis n =  axe::Int (1,2 or 3)

"""
    function rot(phi::Float64,axe::Int)


            n=zeros(3)
            n[axe]=1.0
            R=Mat3(diagm(cos(phi)*(ones(3)-n)+n))+sin(phi)*tilde(Vec3(n))


            return R
        end

"""
    rot(p::RP3)
        function computing the rotation operator R::Mat3
        from the Euler - Rodrigues parameters  p::RP3

"""
    function rot(p::RP3)


                P_0=sqrt(1.0 - (p[1]*p[1]+p[2]*p[2]+p[3]*p[3])/4.0)
                tps = [0.0 -p[3] p[2]; p[3] 0.0 -p[1]; -p[2] p[1] 0.0]

                return R = Mat3(eye(3) + (P_0*eye(3) + 0.5*tps)*tps)

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
