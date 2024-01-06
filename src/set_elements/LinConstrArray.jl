"""
LinConstrArray

        data structure defining a container for linear constraints

"""
mutable struct LinConstrArray

    numbers::Vector{Int}
    node_orders::Vector{Vector{Int}}
    components::Vector{Vector{Int}}
    dofs::Vector{Vector{Int}}
    coefs::Vector{Vector{Float64}}
    val::Vector{Float64}
    scale_factor::Vector{Float64}


    function LinConstrArray()

        numbers = []
        dofs = []
        node_orders = []
        components = []
        coefs = []
        val = []
        scale_factor = []


        return new(numbers, node_orders, components, dofs, coefs, val, scale_factor)

    end


end

lin_constr_container = LinConstrArray()