"""
    set_shell
    
        function to define the topology and structural / geometric / properties of a shell element.

        Defines the topology of a shell  with variable  set of connected nodes and number of Gauss integration points.
        constructs the  "Main.SetElements.shell_container" array  (::ShellArray type) collecting the intial data.

        Input data:

        nbr::Int                                number of the element
        connected_nodes::Vector{Int}            set of nodes of the element
        thickness::Float64                      shell thickness
        stiffness_properties::Vector{Float64}   [E, nu]
        mass_properties::Vector{Float64}        [rho]
        ngauss_points::Int                      number of Gauss points

        visco_type::String                      "none", "maxwell" or "damped"
        ratio_infty::Int                        fraction of elastic stiffness at infinity 
        tau_b::Float64                          viscoeleastic time constant in bending 
        tau_S::Float64                          viscoeleastic time constant in shear


        2 instances of the calling sequence:

        for an elastic shell:

            set_shell(nbr,connected_nodes, thickness, stiffness_properties,mass_properties,ngauss_points)

        for a viscoelastic shell:

            set_shell(nbr,connected_nodes, thickness, stiffness_properties,mass_properties,ngauss_points,
            tau_B,tau_S, ratio_infty, visco_type)

"""  
function set_shell(nbr::Int,connected_nodes::Vector{Int}, thickness::Float64, stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},ngauss_points::Int)
    
    visco_type = "none"
    ratio_infty = 1.0
    tau_B = 0.0
    tau_S = 0.0

    set_shell(nbr,connected_nodes, thickness, stiffness_properties,
    mass_properties,ngauss_points,tau_B,tau_S, ratio_infty, visco_type)

    return
end

function set_shell(nbr::Int,connected_nodes::Vector{Int}, thickness::Float64, stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},ngauss_points::Int,tau_B::Float64,tau_S::Float64)
    
    visco_type = "damped"
    ratio_infty = 1.0

    set_shell(nbr,connected_nodes, thickness, stiffness_properties,
    mass_properties,ngauss_points,tau_B,tau_S, ratio_infty, visco_type)

    return
end

function set_shell(nbr::Int,connected_nodes::Vector{Int}, thickness::Float64, stiffness_properties::Vector{Float64},
    mass_properties::Vector{Float64},ngauss_points::Int,
    tau_B::Float64,tau_S::Float64, ratio_infty::Float64, visco_type::String)

global Nnodes = size(connected_nodes, 1)

ndim = Nnodes*6
ng = 1
visco_type == "maxwell" &&  (ng = ngauss_points) 

global H_init = Vector{NodeFrame}(undef,Nnodes)
global R_rel = Vector{RV3}(undef,Nnodes)
global H_0g = Vector{NodeFrame}(undef,ngauss_points)
global f_0g = Vector{Array{Float64,2}}(undef,ngauss_points)
global J = Vector{Array{Float64,2}}(undef,ngauss_points)
global g_0g = Vector{Vector{Vec3}}(undef,ngauss_points)
global n_0g = Vector{Vec3}(undef,ngauss_points)
global K_0g = Vector{Matrix}(undef,ngauss_points)
global M_NN = zeros(ndim,ndim)
global strains_g = Vector{Matrix}(undef,ng)
global visco_strains_g = Vector{Matrix}(undef,ng)
if visco_type == "maxwell"
    for i = 1:ngauss_points
        visco_strains_g[i] = zeros(12,2)
        strains_g[i] = zeros(12,2)
    end
end

    mc = Main.model_container
    nc = Main.node_container

    if mc.Shells == 0
         global shell_container = ShellArray()
    else
         shell_container = SetElements.shell_container
    end

    sc =   SetElements.shell_container

    append!(sc.numbers,nbr)
    push!(sc.visco_type,visco_type)
    append!(sc.npoints,Nnodes)
    append!(sc.ngauss_points,ngauss_points)
    append!(sc.thickness,thickness)

    time_constants = zeros(12)

    if visco_type != "none" 
        
        nu = stiffness_properties[2]
        tau_shear = tau_S
        tau_ext = 1.0/3.0*((1-2.0*nu)*tau_B + 2.0*(1.0+nu)*tau_S)
        time_constants[1] = tau_ext
        time_constants[2:4] .= tau_shear
        time_constants[5] = tau_ext
        time_constants[6] = 10000.0
        time_constants[7] = tau_shear
        time_constants[8] = tau_ext
        time_constants[9] = tau_shear
        time_constants[10] = tau_ext
        time_constants[11] = tau_shear
        time_constants[12] = 10000.0
         
    end

    size(stiffness_properties,1) != 2 && error("element ", nbr, " wrong dimension of stiffness properties")
    size(mass_properties,1) != 1 && error("element ", nbr, " wrong dimension of mass properties")
    rho = mass_properties[1]
    
    global node_orders = zeros(Int, Nnodes)
    global loc_x = Int[]
    global x_nodes = Vec3[]

    
    for in = 1:Nnodes
        inode = findfirst(x -> x==connected_nodes[in], nc.node_numbers)[1]
        H_init[in] = NodeFrame(copy(nc.init_positions[inode]),copy(nc.init_orientations[inode]))
        node_orders[in] = inode
        push!(x_nodes,nc.init_positions[inode])
        loc_x = [loc_x; nc.locs[inode]]
    end

    x_1 = x_nodes[1]
    x_2 = x_nodes[2]

    (Nnodes == 3 ||  Nnodes == 6) ?  x_3 = x_nodes[3] : x_3 = x_nodes[4]

    n_1 = 1.0/norm2(x_2-x_1)*(x_2-x_1)
    n_2 = 1.0/norm2(x_3-x_1)*(x_3-x_1)
    n_3 = crossp(n_1,n_2)
    n_2 = crossp(n_3,n_1)


    R = RV3(Mat3([n_1.v n_2.v n_3.v]))
 
    for in = 1:Nnodes
        R_rel[in] = RV3(H_init[in].p,-R)
        H_init[in].p  = R
    end

    loc_v = collect(mc.max_v +i for i=1:Nnodes*5)
    mc.max_v += 5*Nnodes

    x_gauss, gauss_weights = gauss_points(2,ngauss_points)

    global shapes = Vector{Vector}(undef,ngauss_points)
    global shape_derivatives = Vector{Array{Float64,2}}(undef,ngauss_points)

    global M_kernel = shell_mass_kernel(rho, thickness)

    area = 0.0
    for ing = 1:ngauss_points
        xi = x_gauss[ing]
        shapes[ing], shape_derivatives[ing] = shape_functions_2D(Nnodes,xi)
        H_0g[ing], f_0g[ing], g_0g[ing], n_0g[ing], J[ing] = frame_solve1(H_init,shapes[ing], shape_derivatives[ing], Nnodes) 
        d_area = gauss_weights[ing]*det(J[ing])
        K_0g[ing] = shell_stiffness_matrix(stiffness_properties, thickness, f_0g[ing])
        K_0g[ing] *= d_area
        for ni = 1:Nnodes
            ix = [(ni-1)*6+i for i in 1:6]
            for nk = 1:Nnodes
                kx = [(nk-1)*6+k for k in 1:6]
                M_NN[ix,kx] += shapes[ing][ni]*shapes[ing][nk]*d_area*M_kernel
            end
        end
        area += d_area
    end


    append!(sc.area,area)
    push!(sc.node_orders,node_orders)
    push!(sc.stiffness_properties,stiffness_properties)
    push!(sc.mass_properties, mass_properties)
    push!(sc.relative_rotations,R_rel) 
    push!(sc.gauss_reference_frames,H_0g) 
    push!(sc.gauss_strains,copy(f_0g)) 
    push!(sc.gauss_jacobians,J) 
    push!(sc.init_deformations,f_0g) 
    push!(sc.gauss_base_vectors,g_0g) 
    push!(sc.gauss_normals,n_0g) 
    push!(sc.x_gauss,x_gauss)
    push!(sc.gauss_weights,gauss_weights)
    push!(sc.shapes,shapes)
    push!(sc.shape_derivatives,shape_derivatives)
    push!(sc.stiffness_matrix,K_0g)
    push!(sc.mass_kernel,M_kernel)
    push!(sc.mass_matrix,M_NN)
    push!(sc.stresses,zeros(12))
    push!(sc.time_constants,time_constants)
    push!(sc.strains,zeros(12,2))
    push!(sc.strains_g,strains_g)
    push!(sc.visco_strains_g,visco_strains_g)
    push!(sc.ratio_infty,ratio_infty)
    push!(sc.visco_type,visco_type)



    mc.Shells  += 1

    loc_int = Int[] 
    loc_mult = Int[]

  append_element(nbr,"shell",node_orders,loc_x,loc_int,loc_v,loc_mult)

end

