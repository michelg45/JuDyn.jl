"""
    set_beam
    
        function to define the topology and structural / geometric / properties of a beam element.

        constructs the  "Main.SetElements.beam_container" array  (::BeamArray type) collecting the intial data.

        Input data:

        nbr::Int                                number of the element
        connected_nodes::Vector{Int}            set of nodes of the element
        ref_orientation::RV3                    orientation of the beam reference frame.
        stiffness_properties::Vector{Float64}   bending and shear stiffnesses per unit length.
        mass_properties::Vector{Float64}        mass and inertia properties per unit length.
        constant_inertia::Bool                  if true, a beam element with contant velocity field will be generated
                                                ("beam_constant_inertia.jl" function).

        visco_type::String                      "none", "visco_QS" or "visco_dyn"
        ratio_infty::Int                        fraction of elastic stiffness at infinity 
        tau_E::Float64                          viscoeleastic time constant in bending 
        tau_S::Float64                          viscoeleastic time constant in shear


        2 instances of the calling sequence:

        for an elastic beam:

            set_beam(nbr,nodes,ref_orientation,stiffness_properties,
            mass_properties,constant_inertia)

            or

            set_beam(nbr,nodes,ref_orientation,stiffness_properties,
            mass_properties)

        for a viscoelastic beam:

            set_beam(nbr,nodes,ref_orientation,stiffness_properties,
            mass_properties,constant_inertia,tau_E,tau_S, ratio_infty, visco_type)

            or 

            set_beam(nbr,nodes,ref_orientation,stiffness_properties,
            mass_properties,tau_E,tau_S, ratio_infty, visco_type)


            if visco_type = "visco_QS", the viscoelastic stress will be computed in a quasi-static manner.
            if visco_type = "visco_dyn", the viscoelastic stress will be integrated in time using additional
            internal variables.   

"""  
function set_beam(nbr::Int,nodes::Vector{Int},ref_orientation::RV3,stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},constant_inertia::Bool)

    visco_type = "none"
    ratio_infty = [1.0]
    tau_E = [0.0]
    tau_S = [0.0]

    set_beam(nbr,nodes,ref_orientation,stiffness_properties,
    mass_properties,constant_inertia,tau_E,tau_S, ratio_infty, visco_type)
    
    return
end

function set_beam(nbr::Int,nodes::Vector{Int},ref_orientation::RV3,stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},constant_inertia::Bool,tau_E::Float64,tau_S::Float64)

    visco_type = "damped"
    ratio_infty = 1.0

    set_beam(nbr,nodes,ref_orientation,stiffness_properties,
    mass_properties,constant_inertia,tau_E,tau_S, ratio_infty, visco_type)
    
    return
end

function set_beam(nbr::Int,nodes::Vector{Int},stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64})

    nc = Main.node_container

    x_1=copy(nc.init_positions[nodes[1]])
    x_2=copy(nc.init_positions[nodes[2]])

    ref_orientation = frame_on_line(x_1,x_2)
    
    constant_inertia = false

    set_beam(nbr,nodes,ref_orientation,stiffness_properties,mass_properties,constant_inertia)
    
    return
end


function set_beam(nbr::Int,nodes::Vector{Int},ref_orientation::RV3,stiffness_properties::Vector{Float64},
        mass_properties::Vector{Float64},constant_inertia::Bool, tau_E::Vector{Float64}, tau_S::Vector{Float64}, ratio_infty::Vector{Float64}, visco_type::String)


    mc = Main.model_container
    nc = Main.node_container


    uniform_rotation = mc.uniform_rotation 
    rotation_speed = mc.rotation_speed


    if mc.Beams == 0
         global beam_container = BeamArray()
    else
         beam_container = SetElements.beam_container
    end

    bc =   SetElements.beam_container

    append!(bc.numbers,nbr)
    push!(bc.visco_type,visco_type)

    if visco_type != "none" 
        
        n_branches = size(tau_E,1)
        
        time_constants = Vector{Any}(undef,n_branches)
        visco_strains =  Vector{Any}(undef,n_branches)

        for nb = 1:n_branches  
            time_constants[nb] = zeros(6)
            visco_strains[nb] = zeros(6,2)
            time_constants[nb][1] = tau_E[nb]
            time_constants[nb][5:6] .= tau_E[nb]
            time_constants[nb][2:4] .= tau_S[nb]
        end
    else 
        time_constants = []
        visco_strains = []
         
    end


    size(stiffness_properties,1) != 6 && error("element ", nbr, " wrong dimension of stiffness properties")
    size(mass_properties,1) != 4 && error("element ", nbr, " wrong dimension of stiffness properties")
    n1=findfirst(x -> x==nodes[1], nc.node_numbers)[1]
    n2=findfirst(x -> x==nodes[2], nc.node_numbers)[1]
    x_1=copy(nc.init_positions[n1])
    x_2=copy(nc.init_positions[n2])
    length = norm2(x_2-x_1)
    PREC=sqrt(eps(Float64))
    if norm(rot(ref_orientation).mat[:,1]-1.0/length*(x_2 - x_1).v) > PREC
        println(1.0/length*(x_2 - x_1).v)
        println(norm(rot(ref_orientation).mat[:,1]-1.0/length*(x_2 - x_1).v))
        error(" beam element ", nbr, "  reference frame not aligned on beam axis")
    end
    psi_rel = Vector{RV3}(undef,2)
#    psi_rel[1]=RV3(-ref_orientation, nc.init_orientations[n1])
#    psi_rel[2]=RV3(-ref_orientation,nc.init_orientations[n2])

    psi_rel[1]=RV3(-nc.init_orientations[n1],ref_orientation)
    psi_rel[2]=RV3(-nc.init_orientations[n2],ref_orientation)

    push!(bc.node_orders,[n1,n2])
    push!(bc.length,norm((x_2-x_1).v))
    push!(bc.local_node_orientations,psi_rel)
    push!(bc.stiffness_properties,stiffness_properties)
    push!(bc.mass_properties, mass_properties)
    push!(bc.constant_inertia,constant_inertia)
    push!(bc.stresses,zeros(6))
    push!(bc.visco_type,visco_type)
    push!(bc.time_constants,time_constants)
    push!(bc.strains,zeros(6,2))
    push!(bc.visco_strains,visco_strains)
    push!(bc.ratio_infty,ratio_infty)
    



    mc.Beams  += 1

    loc_x = [copy(nc.locs[n1]);copy(nc.locs[n2])]
    if constant_inertia == true 
        loc_v = collect(mc.max_v +i for i=1:6)
        mc.max_v += 6
    else
        loc_v = collect(mc.max_v +i for i=1:12)
        mc.max_v += 12
    end

    loc_int = Int[]
    loc_mult = Int[]

    append_element(nbr,"beam",[n1,n2],loc_x,loc_int,loc_v,loc_mult)

end

function set_beam(nbr::Int,nodes::Vector{Int},ref_orientation::RV3,stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64})
    constant_inertia = false
    set_beam(nbr,nodes,ref_orientation,stiffness_properties,
    mass_properties,constant_inertia)
end

function set_beam(nbr::Int,nodes::Vector{Int},ref_orientation::RV3,stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},tau_E::Vector{Float64}, tau_S::Vector{Float64}, ratio_infty::Vector{Float64}, visco_type::String)
    constant_inertia = false
    set_beam(nbr,nodes,ref_orientation,stiffness_properties,
    mass_properties,constant_inertia,tau_E,tau_S,ratio_infty,visco_type)
end