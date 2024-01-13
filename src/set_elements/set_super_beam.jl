"""
    set_super_beam

        Function defining the topology of a straight beam between nodes nodes[1] and nodes[2]
        The node of reference 'ref_node' corresponds the floating center of mass.
        The first reference axis at ref_node must be aligned with the beam reference line.
        The other two reference axes define the cross-section orientation.

        It constructs the array super_beam_container[iel]  containing

        * the geometric and topological data
        * the elements  describing the linear behavior of the element  in the local frame: 
            mass::Float64               total mass of the element.
            Jrot::Array{Float64}        inertia tensor of the element.
            K_elast::Array{Float64}     reduced linear stiffness matrix.
            M_elast::Array{Float64}     reduced linear mass matrix associated to elastic modes.
            S::Array{Float64}           matrix allowing to compute the deformation-dependent contribution to elastic forces.    

        To get more infomation on the contruction method, see function 'super_beam_matrix_kernel.jl'.

        input: 
            nbr::Int                                element number
            ref_node::Int                           reference node at mid-length
            nodes::Vector{Int}                      end nodes of the beam element
            stiffness_properties::Vector{Float64}   stiffness properties of the element mass_properties::Vector{Float64}        mass properties of the element
            nl_correction::Bool                     nonlinear correction parameter (geometric  stiffness correction) : 'true' if omitted.

        Calling sequences: 
            set_super_beam(nbr,ref_node,nodes,stiffness_properties,mass_properties,nl_correction)
            set_super_beam(nbr,ref_node,nodes,stiffness_properties,mass_properties)

"""
function set_super_beam(nbr::Int,ref_node::Int,nodes::Vector{Int}, stiffness_properties::Vector{Float64},mass_properties::Vector{Float64},nl_correction::Bool)


    mc = Main.model_container
    nc = Main.node_container



    if mc.SuperBeams == 0
         global super_beam_container = SuperBeamArray()
    else
         super_beam_container = SetElements.super_beam_container
    end

    sbc =   SetElements.super_beam_container

    append!(sbc.numbers,nbr)

    size(stiffness_properties,1) != 6 && error("element ", nbr, " wrong dimension of stiffness properties")
    size(mass_properties,1) != 4 && error("element ", nbr, " wrong dimension of stiffness properties")
    n1=findfirst(x -> x==nodes[1], nc.node_numbers)[1]
    n2=findfirst(x -> x==nodes[2], nc.node_numbers)[1]
    n_ref=findfirst(x -> x == ref_node, nc.node_numbers)[1]
    x_1=copy(nc.init_positions[n1])
    x_2=copy(nc.init_positions[n2])
    x_0=copy(nc.init_positions[n_ref])
    psi_0=copy(nc.init_orientations[n_ref])
    psi_rel = Vector{RV3}(undef,2)
    psi_rel[1]=RV3(-psi_0, nc.init_orientations[n1])
    psi_rel[2]=RV3(-psi_0, nc.init_orientations[n2])
    """R_1 = rot(-psi_rel[1]).mat
    R_2 = rot(-psi_rel[2]).mat
    R = cat(R_1,R_1,R_2,R_2; dims = (1,2))"""

    x_mid = 0.5*(x_1 + x_2)

    length = norm((x_2-x_1).v)

    PREC=sqrt(eps(Float64))

    if norm((x_mid-x_0).v)/length > PREC
        error(" super_beam element ", nbr, " incorrect position of reference node")
    elseif norm(rot(psi_0).mat[:,1]-1.0/length*(x_2 - x_1).v) > PREC
        error(" super_beam element ", nbr, "  reference frame not aligned on beam axis")
    end

    append!(sbc.ref_node_orders,n_ref)
    push!(sbc.boundary_node_orders,[n1,n2])
    append!(sbc.length,length)
    push!(sbc.local_node_orientations,psi_rel)
    mass, Jrot,  K_elast, M_elast, S = super_beam_matrix_kernel(length, stiffness_properties, mass_properties)

    append!(sbc.mass,mass)
    push!(sbc.Jrot, Jrot)
    push!(sbc.K_elast, K_elast)
    push!(sbc.M_elast, M_elast)
    push!(sbc.S, S)
    append!(sbc.nl_correction,nl_correction)

    mc.SuperBeams  += 1

    loc_x = [copy(nc.locs[n1]);copy(nc.locs[n2]);copy(nc.locs[n_ref])]
    loc_v = collect(mc.max_v +i for i=1:18)
    mc.max_v += 18
    loc_int = Int[]
    loc_mult = Int[]
    nodes = [nodes; ref_node]
    append_element(nbr,"super_beam",nodes,loc_x,loc_int,loc_v,loc_mult)

end

function set_super_beam(nbr::Int,ref_node::Int,nodes::Vector{Int},stiffness_properties::Vector{Float64},mass_properties::Vector{Float64})

"""
    if not specified, the geometric correction is automatically included
"""
    nl_correction = true
    set_super_beam(nbr,ref_node,nodes,stiffness_properties,mass_properties,nl_correction)
end
