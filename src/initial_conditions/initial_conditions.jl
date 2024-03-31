"""
    initial_conditions
"""
function  initial_conditions()

    Ndofs=copy(Main.model_container.Ndofs)

    global  x_0=Vector{Float64}(zeros(Ndofs))
    global  xdot_0=Vector{Float64}(zeros(Ndofs))
    global  Dx_0 = Vector{Float64}(zeros(Ndofs))  


    get_initial_position()

end

function  initial_conditions(init_file::String,sol_type::String)

    mc = Main.model_container

    Ndofs=copy(mc.Ndofs)

    global  x_0=Vector{Float64}(zeros(Ndofs))
    global  xdot_0= Vector{Float64}(zeros(Ndofs))
    global  Dx_0 = Vector{Float64}(zeros(Ndofs))  
    x_1 = Vector{Float64}(zeros(Ndofs))

    h5_file = init_file*".h5"
    jld_file = init_file*".jld"

    println("Initial position read from file ", h5_file)


    if sol_type == "static"
        times, y, ydot, p, str_energy, pot_energy, ext_work, node_positions, nitmax, model_data = read_results_static(h5_file)
    elseif sol_type == "dynamic" 
       times, y, ydot, p, str_energy, pot_energy, kin_energy, ext_work, node_positions,  nitmax, model_data = read_results_dynamic(h5_file)
    end
   
    Nodes = mc.Nodes
    Ndofs_x = mc.Ndofs_x

    Nodes  != model_data[1] && (error("incompatible initial model"))

    istep = size(times,1) 

    x_1[1:Ndofs_x] = y[istep,1:Ndofs_x]
    x_0[1:Ndofs_x] = y[istep-1,1:Ndofs_x]
    

   
    nc = Main.node_container
    
    psi_0 = zeros(3)
    psi_1 = zeros(3)
    Dx = zeros(3)    


    for inode = 1:Nodes
        psi_0 .= 0.0
        psi_1 .= 0.0
        Dx .= 0.0 
        inv_x = nc.inv_locs[inode][1:3]
        inv_psi = nc.inv_locs[inode][4:6]
        for ic = 1:3
            inv_xc = inv_x[ic]
            inv_xc > 0 && (Dx[ic] = x_1[inv_xc] -x_0[inv_xc])
            inv_psic = inv_psi[ic] 
            if inv_psic > 0
                psi_0[ic] = x_0[inv_psic]
                psi_1[ic] = x_1[inv_psic]
            end
        end

        Dpsi = RV3(RV3(psi_1),RV3(-psi_0))
        for ic = 1:3
            inv_x[ic] > 0 && (Dx_0[inv_x[ic]] = Dx[ic])
            inv_psi[ic] > 0 && (Dx_0[inv_psi[ic]] = Dpsi[ic])
        end
    end


    mc.initial_shape = true
    mc.init_file = init_file

end