"""
    Adj
"""
function Adj(H::NodeFrame)
    R = rot(H.p)
    return  [R.mat (tilde(H.x)*R).mat; zeros(3,3) R.mat]
    
end

function Adj(H::NodeFrame, a::Vector)
    a_u = Vec3(a[1:3])
    a_omega = Vec3(a[4:6])
    return  [(rot(H.p, a_u) + crossp(H.x, rot(H.p, a_omega))).v; (rot(H.p, a_omega)).v]
end

function Adj_inv(H::NodeFrame)
    Rt = transpose(rot(H.p))
    return  [Rt.mat -(Rt*tilde(H.x)).mat; zeros(3,3) Rt.mat]
    
end

