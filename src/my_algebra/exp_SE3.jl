"""
    exp_SE3
"""
function exp_SE3(p::Vector)
    (size(p,1) != 6 ) && (error("exp_SE3: Vector p and dp must have size 6"))
    rv = RV3(p[4:6])
    u = transpose(tang(rv))*Vec3(p[1:3])
    return H = NodeFrame(u,rv)
end