"""
    inv_loc

            function returning the location of a set of dofs "locel" into a structural set "Loc".  
            The location of missing dofs in the structural set is set to 0.

            Input: 
                locel::Vector::Vector{Int}  set of dof numbers (e.g,. element set)
                Loc::Vector{Int}            structural set of dofs.

            Output:
                invloc::Vector::Vector{Int} position of "locel" dofs into the structural set "Loc".

            calling sequence:
                invloc = inv_loc(locel,Loc)
"""
function inv_loc(locel::Vector{Int},Loc::Vector{Int})
    locn = copy(locel)
    sz = size(locn,1)
    for i in 1:sz
    a=findfirst(ix->ix==locn[i],Loc)
    typeof(a) == Nothing ? locn[i]=0 : locn[i]=a
    end
    return locn
end