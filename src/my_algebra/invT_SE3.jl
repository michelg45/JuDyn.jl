"""
    invT_SE3
"""
function invT_SE3(p::Vector)
    size(p,1) != 6 && (error("Vector p must have size 6"))
    rv = RV3(p[4:6])
    u = Vec3(p[1:3])
    invT = invtang(rv).mat
    return invT_SE3 = [invT T_uam(u,rv).mat; zeros(3,3) invT]
end

function T_uap(u::Vec3,psi::RV3)

    PSI2 = transpose(psi.v)*psi.v
    PSI = sqrt(PSI2)
    psixu = transpose(psi.v)*u.v
    if (PSI > PREC)
        alpha = sin(PSI)/PSI
        beta = 2.0*(1 - cos(PSI))/PSI2
        T_uap = -0.5*beta*tilde(u) + (1.0-alpha)/PSI2*Mat3(psi.v*(u.v)'+ u.v*(psi.v)' -(2.0*psixu*eye(3))) +
                psixu/PSI2*((beta-alpha)*tilde(psi)+(0.5*beta -3.0*(1.0-alpha))/PSI2*Mat3(psi.v*(psi.v)'-PSI2*eye(3)))
    else
        bracket = Mat3(psi.v*(u.v)'+ u.v*(psi.v)' -(2.0*psixu*eye(3)))
        tpsi2 = Mat3(psi.v*(psi.v)'-PSI2*eye(3))
        T_uap = -0.5*tilde(u) +1.0/6.0*bracket - 1.0/24.0*(tilde(psi)*bracket + 
                tilde(u)*tpsi2) + 1.0/120.0*(bracket*tpsi2+tpsi2*bracket)
    end
    return T_uap
end

function T_uam(u::Vec3,psi::RV3)
    
        PSI2 = transpose(psi.v)*psi.v
        PSI = sqrt(PSI2)
        psixu = transpose(psi.v)*u.v
        PREC2 = 1.e-2
        if (PSI > PREC2)
            beta = 2.0*(1 - cos(PSI))/PSI2
            gamma = PSI/2.0*cot(PSI/2.0)
            T_uam = 0.5*tilde(u) + (1.0-gamma)/PSI2*Mat3(psi.v*(u.v)'+ u.v*(psi.v)' -(2.0*psixu*eye(3))) +
            psixu/(PSI2*PSI2)*(1.0/beta + gamma - 2.0)*Mat3(psi.v*(psi.v)'-PSI2*eye(3))
        else
            T_uam = 0.5*tilde(u)
            tmp = tilde(psi) * tilde(u) + tilde(u) * tilde(psi)
            T_uam += 1. / 12. *tmp
            tpw2 = tilde(psi) * tilde(psi)
            T_uam += 1. / (-30. * 24.) * (tpw2 * tmp + tmp * tpw2)
        end

    return T_uam
end


