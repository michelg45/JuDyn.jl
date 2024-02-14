"""
    distance_to_point
"""
function distance_to_point(x::Vec3,p1::Vec3)
    r = (x - p1).v
    g = sqrt(r'*r)
    dgdx = r/g
    return g, dgdx
end