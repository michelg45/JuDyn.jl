function A_SE3(p::Vector,y::Vec3)
    bracket(u,v) = tilde(u)*tilde(v)+tilde(v)*tilde(u)
    x_u = Vec3(p[1:3])
    x_omega = Vec3(p[4:6])
    PSI = norm2(x_omega)
    PSI2 = PSI*PSI
    if PSI > PREC 
        gamma = 2.0*cot(PSI/2.0)/PSI
        alpha = sin(PSI)/PSI
        beta = 2.0*(1- cos(PSI))/PSI2
    else
        gamma = 1.0
        alpha = 1.0
        beta = 0.0
    end
    c1 = (1.0 - gamma)/PSI2
    c2 = (gamma-1.0)*(gamma+2.0)/PSI2^2 + 0.25*PSI2
end

