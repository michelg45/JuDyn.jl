"""
    point_along_direction

        function that calculates the projection 'h' of the position of a node x along a direction 'axe'
        and its derivative 'dddx'

        Input:

        x::Vec3                 node position
        axe::Vec3               direction

        calling sequence:

            h, dhdx = point_along_direction(x,axe)


"""
function point_along_direction(x::Vec3,axe::Vec3)
 
    h = dotp(x,h)
    dhdx = axe.v
    return h, dhdx
end