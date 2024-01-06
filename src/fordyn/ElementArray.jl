"""
    ElementArray

        ElementArray is a data structure that contains the general data associated to elements. 
        It is defined in the SetElements module.

        the "Main.element_container" is created by the "create_model.jl" function.

        calling sequence:
            
            global element_container = ElementArray()

        It contains the following data fields:
            element_numbers::Vector{Int}            numbering of elements
            element_types::Array{String,1}          element types
            element_nodes::Array{Vector{Int},1}     nodes of elements
            n_x::Vector{Int}                        number of dofs of "displacement type"
            n_int::Vector{Int}                      number of internal dofs
            n_v::Vector{Int}                        number of velocity dofs                      
            n_mult::Vector{Int}                     number of Lagrange multipliers
            loc_x::Array{Vector{Int},1}             localization of "displacement" dofs
            loc_int::Array{Vector{Int},1}           localization of internal dofs
            loc_v::Array{Vector{Int},1}             localization of velocity dofs
            loc_mult::Array{Vector{Int},1}          localization of Lagrange multipliers
            inv_loc_x::Array{Vector{Int},1}         inverse localization of "displacement" dofs
            inv_loc_int::Array{Vector{Int},1}       inverse localization of internal dofs
            inv_loc_v::Array{Vector{Int},1}         inverse localization of velocity dofs
            inv_loc_mult::Array{Vector{Int},1}      inverse localization of Lagrange multipliers
        
            note: loc_... contains the dof numbers, while inv_loc_... provides the positiosn of the
            dofs in the structural set. 
"""
mutable struct  ElementArray

    element_numbers::Vector{Int}
    element_types::Array{String,1}
    element_nodes::Array{Vector{Int},1}
    n_x::Vector{Int}
    n_int::Vector{Int}
    n_v::Vector{Int}
    n_mult::Vector{Int}
    loc_x::Array{Vector{Int},1}
    loc_int::Array{Vector{Int},1}
    loc_v::Array{Vector{Int},1}
    loc_mult::Array{Vector{Int},1}
    inv_loc_x::Array{Vector{Int},1}
    inv_loc_int::Array{Vector{Int},1}
    inv_loc_v::Array{Vector{Int},1}
    inv_loc_mult::Array{Vector{Int},1}

end

"""
    ElementArray()

    function creating an array of ElementArray type.

        global element_container = ElementArray()
        
"""
    ElementArray()=ElementArray(Vector{Int}[],Array{String,1}[],Array{Vector{Int},1}[],Vector{Int}[],Vector{Int}[],Vector{Int}[],Vector{Int}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[],Array{Vector{Int},1}[])

