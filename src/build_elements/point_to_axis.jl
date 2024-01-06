"""
    point_to_axis

        function that calculates the distance 'g' between a node x and an axis(p1,p2) 
        and its derivative 'dgdx'

        Input:

        x::Vec3                 node position
        p1::Vec3,p2::Vec3       set of nodes defining the axis

        calling sequence:

            g, dgdx = point_to_axis(x,p1,p2)


"""
function point_to_axis(x::Vec3,p1::Vec3,p2::Vec3)
    r = (x - p1).v
    s = (p2 - p1).v
    a = r'*s
    b = sqrt(r'*r)
    c = sqrt(s'*s)
    g = 1/c*sqrt(b^2*c^2-a^2)
    dgda = -a/(c^2*g)
    dgdb = b/g
    dbdx = r/b
    dadx = s
    dgdx = dgda*s + dgdb*dbdx
    return g, dgdx
end