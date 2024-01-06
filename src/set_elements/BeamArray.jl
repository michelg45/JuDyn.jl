"""
    BeamArray

        Data structure for the "Main.SetElements.beam_container" array. Contains the following data
        for each beam element of the element set.

            numbers::Vector{Int}                                    number of the element
            node_orders::Vector{Vector{Int}}                        order of the element nodes in th structural node set
            local_node_orientations::Array{Vector{RV3},1}           node orientations relative to beam orientation
            length::Vector{Float64}                                 beam length
            stiffness_properties::Vector{Vector{Float64}}           element stiffness properties
            time_constants::Array{Vector{Float64},1}                viscoelastic time constants
            ratio_infty::Vector{Float64}                            fraction of elastic stress at infinity
            mass_properties::Vector{Vector{Float64}}                element mass_properties             
            constant_inertia::Vector{Bool}                          constant or linear inertia distribution
            stresses::Vector{Vector{Float64}}                       (6x1) stress vector
            visco_stresses::Vector{Array{Float64,2}}                viscous stresses  ("visco_QS" beam).
            strains::Vector{Array{Float64,2}}                       (6x1) strain vector at two successive times
            visco_type::Vector{String}                              viscoelastic type of the element
    
            creation sequence:
    
                global beam_container = BeamArray() 
"""
mutable struct BeamArray


    numbers::Vector{Int}
    node_orders::Array{Vector{Int},1}
    local_node_orientations::Array{Vector{RV3},1}
    length::Vector{Float64}
    stiffness_properties::Array{Vector{Float64},1}
    time_constants::Array{Vector{Float64},1}
    ratio_infty::Vector{Float64}
    mass_properties::Array{Vector{Float64},1}
    constant_inertia::Vector{Bool}
    stresses::Vector{Vector{Float64}}
    visco_strains::Vector{Array{Float64,2}}
    strains::Vector{Array{Float64,2}}
    visco_type::Vector{String}

    function BeamArray()
        
        numbers = []
        node_orders = []
        local_node_orientations = []
        length = []
        stiffness_properties = []
        time_constants = []
        ratio_infty = []
        mass_properties = []
        constant_inertia = []
        stresses = []
        visco_strains = []
        strains = []
        visco_type = []
        return new(numbers,node_orders,local_node_orientations,length,stiffness_properties,time_constants,ratio_infty,mass_properties,constant_inertia,stresses,visco_strains,strains,visco_type)
        
    end

end