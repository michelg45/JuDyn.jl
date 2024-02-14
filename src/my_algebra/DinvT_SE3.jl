"""
    DinvT_SE3
"""
function DinvT_SE3(p::Vector, a::Vector)
    
    prec = 1.e-6
    pos = p[1:3]
    norm_x = norm(pos)

    norm_x <= prec &&  (norm_x = 1.0)

    D = zeros(6,6)
     
    for k = 1:6
        k <= 3 ? inc = prec*norm_x : inc = prec
        dp = zeros(6)
        dp[k] = inc
        v = 0.5/inc*(invT_SE3(p+dp) - invT_SE3(p-dp))*a
        D[:,k] = v[:]
    end

    return D

end


