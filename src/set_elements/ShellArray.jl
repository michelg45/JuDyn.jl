"""
    ShellArray

        Data structure for the "Main.SetElements.shell_container" array. Contains the following data
        for each shell element of the element set.

        numbers::Vector{Int}                                    number of the element
        npoints::Vector{Int}                                    number of nodes of the element
        ngauss_points::Vector{Int}                              number of Gauss points
        node_orders::Vector{Vector{Int}}                        order of the element nodes in the structural node set
        stiffness_properties::Vector{Vector{Float64}}           element stiffness properties
        time_constants::Array{Vector{Float64},1}                viscoelastic time constants
        ratio_infty::Vector{Float64}                            fraction of elastic stress at infinity
        mass_properties::Vector{Vector{Float64}}                element mass_properties             
        thickness::Vector{Float64}                              element thickness
        area::Vector{Float64}                                   surface of the element
        relative_rotations::Vector{Vector{RV3}}                 relative rotations at nodes with respect to shell plane
        gauss_reference_frames::Vector{Vector{NodeFrame}}       reference frames at Gauss points
        gauss_strains::Vector{Vector{Array{Float64,2}}}         strains at Gauss points
        gauss_jacobians::Vector{Vector{Array{Float64,2}}}       Jacobians at Gauss points
        gauss_base_vectors::Vector{Vector{Vector{Vec3}}}        base vectors at  Gauss points
        gauss_normals::Vector{Vector{Vec3}}                     normals at gauss points
        init_deformations::Vector{Vector{Array{Float64,2}}}     (6x2) matrix of initial deformations
        shapes::Vector{Vector{Vector{Float64}}}                 element shape functions
        shape_derivatives::Vector{Vector{Array{Float64,2}}}     derivatives of element shape functions
        gauss_weights::Vector{Vector{Float64}}                  weights at Gauss points
        x_gauss::Vector{Vector{Vector{Float64}}}                coordinates of Gauss points
        stiffness_matrix::Vector{Vector{Matrix}}                stiffness kernel of the element
        mass_kernel::Vector{Matrix}                             
        mass_matrix::Vector{Matrix}                             mass kernel of the element
        stresses::Vector{Vector{Float64}}                       vector of mean element stresses
        strains::Vector{Array{Float64,2}}                       vector of mean element strains at two successive times
        strains_g::Vector{Vector{Matrix}}                       strains at Gauss points
        visco_type::Vector{String}                              viscoelastic type of the element

        creation sequence:

            global shell_container = ShellArray() 

"""
mutable struct ShellArray

    numbers::Vector{Int}
    npoints::Vector{Int}
    ngauss_points::Vector{Int}
    node_orders::Vector{Vector{Int}}
    stiffness_properties::Vector{Vector{Float64}}
    ratio_infty::Vector{Any}
    mass_properties::Vector{Vector{Float64}}
    thickness::Vector{Float64}
    area::Vector{Float64}
    relative_rotations::Vector{Vector{RV3}}
    gauss_reference_frames::Vector{Vector{NodeFrame}}
    gauss_strains::Vector{Vector{Array{Float64,2}}}
    gauss_jacobians::Vector{Vector{Array{Float64,2}}}
    gauss_base_vectors::Vector{Vector{Vector{Vec3}}}
    gauss_normals::Vector{Vector{Vec3}}
    init_deformations::Vector{Vector{Array{Float64,2}}}
    shapes::Vector{Vector{Vector{Float64}}}
    shape_derivatives::Vector{Vector{Array{Float64,2}}}
    gauss_weights::Vector{Vector{Float64}}
    x_gauss::Vector{Vector{Vector{Float64}}}
    stiffness_matrix::Vector{Vector{Matrix}}
    mass_kernel::Vector{Matrix}
    mass_matrix::Vector{Matrix}
    stresses::Vector{Any}
    time_constants::Vector{Any}
    visco_strains_g::Vector{Vector{Any}}
    strains::Vector{Any}
    strains_g::Vector{Vector{Any}}
    visco_type::Vector{String}


    function ShellArray()

        numbers = []
        npoints = []
        ngauss_points = []
        node_orders = []
        stiffness_properties = []
        time_constants = []
        ratio_infty = []
        mass_properties = []
        thickness = []
        area = []
        relative_rotations = []
        gauss_reference_frames = []
        gauss_strains = []
        gauss_jacobians = []
        gauss_base_vectors = []
        gauss_normals = []
        init_deformations = []
        shapes = []
        shape_derivatives = []
        gauss_weights = []
        x_gauss = []
        stiffness_matrix = []
        mass_kernel = []
        mass_matrix = []
        stresses = []
        strains = []
        visco_strains_g = []
        strains_g = []
        visco_type = []
        
        return new(numbers, npoints, ngauss_points, node_orders, stiffness_properties,
            time_constants,ratio_infty,mass_properties, thickness,  area, relative_rotations, gauss_reference_frames, 
            gauss_strains, gauss_jacobians, gauss_base_vectors,
            gauss_normals, init_deformations, shapes,
            shape_derivatives, gauss_weights, x_gauss, stiffness_matrix, mass_kernel,
            mass_matrix,stresses,visco_strains_g,strains,strains_g,visco_type)
        
    end
end