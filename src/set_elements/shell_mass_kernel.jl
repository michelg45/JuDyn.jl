
"""
    shell_mass_kernel

        Function constructing the mass kernel of a shell element.

        Calling sequence:

            shell_mass_kernel(rho,thickness)

"""
function shell_mass_kernel(rho::Float64,thickness::Float64)

    I3 = diagm(ones(3))

    global M_k = zeros(6,6)
    global C = zeros(6,6)

    C[1:3,1:3] = I3

    
    global e3 = Vec3(0.0, 0.0,1.0)

    global xg, wg = gauss_points(1,3)
    
    for ig = 1:3
        beta = 0.5*thickness*xg[ig]
        C[1:3, 4:6] = -beta*tilde(e3).mat
        C[4:6, 1:3] =   transpose(C[1:3, 4:6]) 
        C[4:6, 4:6] = diagm([beta^2; beta^2; 0.0])
        M_k += wg[ig]*C
    end

    return M_k *= 0.5*thickness*rho

end