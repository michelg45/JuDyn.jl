
"""
    shell_stiffness_matrix

        Function constructing the stiffness kernel of a shell element.

        Calling sequence:

            shell_stiffness_matrix(stiffness_properties,thickness)

"""
function shell_stiffness_matrix(stiffness_properties::Vector{Float64}, thickness::Float64,F_0::Matrix)
        
        
    E = stiffness_properties[1]
    nu = stiffness_properties[2]

    K = zeros(12,12)

    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    lambda_bar = 2.0 * lambda * mu / (lambda + 2.0 * mu)

    ng1 = norm(F_0[1:3,1])
    ng2 = norm(F_0[1:3,2])

    e2o12 = thickness * thickness / 12.0
    tmp = 1.0 / (ng1 * ng1) * thickness
    K[1, 1] = tmp * (lambda_bar + 2.0 * mu)
    K[2, 2] = tmp * mu
    K[3, 3] = tmp * mu
    K[4, 4] = tmp * mu * e2o12
    K[5, 5] = tmp * (lambda_bar + 2.0 * mu) * e2o12

    tmp = 1. / (ng1 * ng2) * thickness
    K[1, 8] = tmp * lambda_bar
    K[8, 1] = K[1, 8]
    K[2, 7] = tmp * mu
    K[7, 2] = K[2, 7]
    K[4, 11] = -tmp * mu * e2o12
    K[11, 4] = K[4, 11]
    K[5, 10] = -tmp * lambda_bar * e2o12
    K[10, 5] = K[5, 10]

    tmp = 1. / (ng2 * ng2) * thickness
    K[7, 7] = tmp * mu
    K[8, 8] = tmp * (lambda_bar + 2.0 * mu)
    K[9, 9] = tmp * mu
    K[10, 10] = tmp * (lambda_bar + 2.0 * mu) * e2o12
    K[11, 11] = tmp * mu * e2o12

    return K
end