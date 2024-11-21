

"""
dynamic_solve

    Function performing theta time integration of a system describd in first-order form
    using the generalized - alpha method.
    Its parameters are read from the file JSON_file.

        T = period of analysis
        k = scaling factor
        PREC = precision parameter
        NitMax = maximum number of iterations
        Npas = number of time steps
        rho_inf = spectral radius at infinity
        elimv = true: elimination of velocities before solution of correction equation.
        verbose = true: print at end of each time step.
        problem_name = name definition of the HDF5 file.
        save = true: output to HDF5 file.
        energy_balance = true: re-evaluation of energies at the end of the time step.

"""
function dynamic_solve(JSON_file)

    #
    # timer initialization
    #

    my_timer()  =  (Dates.hour(now())*3600+Dates.minute(now())*60+Dates.second(now()))*1000+Dates.millisecond(now())   

     t_start = my_timer()

    #
    # Reading solution and output parameters from JSON file
    #


    open(JSON_file,"r") do f
           global dicts
           dicts = read(f,String)
    end

    dict =  JSON.parse(dicts)



    T = dict["Time_Integrator"]["parameters"]["T"]
    k = dict["Time_Integrator"]["parameters"]["k"]
    PREC = dict["Time_Integrator"]["parameters"]["PREC"]
    NitMax = dict["Time_Integrator"]["parameters"]["NitMax"]
    Npas = dict["Time_Integrator"]["parameters"]["Npas"]
    rho_inf = dict["Time_Integrator"]["parameters"]["rho_inf"]
    verbose = dict["Time_Integrator"]["parameters"]["verbose"]

    is_eigvals = haskey(dict,"Eigenvalue_Analysis")
    if is_eigvals == false
        eigvals = false
    else
        eigvals = dict["Eigenvalue_Analysis"]["eigvals"]
        eigvals == true ?  (eig_freq = dict["Eigenvalue_Analysis"]["eig_freq"];
        max_vals = dict["Eigenvalue_Analysis"]["max_vals"]) : (max_vals = 1; eig_freq = 0)
    end

    
    problem_name = dict["Output"]["name"]
    save = dict["Output"]["parameters"]["save_h5"]
    save_freq = dict["Output"]["parameters"]["save_freq"]
    energy_balance = dict["Output"]["parameters"]["energy_balance"]

    

    # external_force = Main.InputFunctions.input_functions[dict["Input"]["external_force"]]

    #
    # computing solution parameters
    #


    global h = T/(Npas-1)

    global delta_m=0.5*(3.0*rho_inf-1.0)/(rho_inf+1.0)
    global delta_f=rho_inf/(rho_inf+1.0)
    global theta=0.5 + delta_f- delta_m
    global theta_p=(1.0-delta_m)/(h*theta*(1.0-delta_f))

    #
    # output files. mat_file can be crated only if save = true
    #

    sh5 = ".h5"

    h5_file = problem_name*sh5

    isfile(h5_file) && rm(h5_file)


    #
    # Reading model topology  from the diffrent containers
    #

    element_container = Main.element_container
    model_container = Main.model_container
    Nnodes = model_container.Nodes
    Ndofs = model_container.Ndofs
    Ndofs_x= model_container.Ndofs_x
    Ndofs_int= model_container.Ndofs_int
    Ndofs_v= model_container.Ndofs_v
    Ndofs_mult= model_container.Ndofs_mult
    Ndofs_q = model_container.Ndofs_q
    Nel = model_container.Elements
    matrix_update = model_container.matrix_update
    initial_shape =  model_container.initial_shape
    element_numbers = element_container.element_numbers
    element_types = element_container.element_types



    #
    # Initialization ov vectors and arrays
    #

    global p = Vector{Float64}(zeros(Ndofs))
    global res = Vector{Float64}(zeros(Ndofs))
    global y_n = Vector{Float64}(zeros(Ndofs))
    global ydot_n = Vector{Float64}(zeros(Ndofs))
    global ydot_np1 = Vector{Float64}(zeros(Ndofs))
    global Dy = Vector{Float64}(zeros(Ndofs))
    global dy = Vector{Float64}(zeros(Ndofs))
    global wy = Vector{Float64}(zeros(Ndofs))

    y_n[:] = Main.InitialConditions.x_0[:]
    Dy[:] = Main.InitialConditions.Dx_0[:]

    #
    # Index vectors of sparse system matrix
    #

    global sps = sparse_matrix_S

    I_qq = sps.I_index_qq
    J_qq = sps.J_index_qq
    I_vq = sps.I_index_vq
    J_vq = sps.J_index_vq

    nz_qq = size(I_qq,1)
    nz_vq = size(I_vq,1)

#
#   pre-factorization of system matrix
#
    global F = lu(sparse(I_qq,J_qq,rand(nz_qq)))
    global S_vq = sparse(I_vq,J_vq,sps.S_vq)

    #
    # Ordering of degrees of freedom for the global system
    # _x :      inematic variables at node frames
    # _int:     internal variables to elements
    # _mult:    multipliers
    # _v:       velocities (defined at element level)
    #

    iloc_x=[i for i in 1:Ndofs_x]
    Ndofs_int > 0 ? iloc_int=[i for i in Ndofs_x+1:Ndofs_x+Ndofs_int] :  iloc_int=Int[]
    Ndofs_mult > 0 ? iloc_mult=[i for i in Ndofs_x+Ndofs_int+1:Ndofs_x+Ndofs_int+Ndofs_mult] :  iloc_mult=Int[]
    iloc_v=[i for i in Ndofs_q+1:Ndofs_q+Ndofs_v]

    iloc_xiv = [iloc_x; iloc_int; iloc_v]
    iloc_iv = [iloc_int; iloc_v]
    iloc_q = [iloc_x;iloc_int;iloc_mult]
    iloc_xi = [iloc_x;iloc_int]

    
    time_elements = 0.0

    init_frames(y_n, ydot_n,initial_shape) 

    

    itime = 1
    
    global niter = 0

    global ext_work = 0.0
    global sum_ext_work = 0.0

    # global inc_ext_work = 0.0

    global itime_vals_saved = 0

    global itime_saved = 0

    res[:] .= 0.0

    kin_energy, pot_energy, str_energy = dynamic_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,theta_p,ext_work,itime,h,niter)

    increment_frames(Dy,y_n,Nnodes)

    Dy[:] .= 0.0

    if save == true
        file, dsets = create_h5_file_dynamic(h5_file,Npas,save_freq,eigvals,eig_freq,max_vals)

        vals = zeros(max_vals)
        times = (itime -1)*h

        itime_saved, itime_vals_saved = record_on_h5_file_dynamic(dsets,itime,itime_saved,itime_vals_saved,times,niter,y_n,ydot_np1,p,
            pot_energy,kin_energy,str_energy,ext_work,Npas,save,save_freq,eigvals,eig_freq, vals)


    end

     

    #
    # pseudo-velocity wy is initially set to 0.
    #

    wy[iloc_xiv]=(1.0 - delta_f)/(1.0 - delta_m)*ydot_n[iloc_xiv] 

        println("start solve")
        println("number of time steps: ", Npas)
        println("spectral radius: ", rho_inf)

        nit_tot = 0

    for itime = 2:Npas

        times = (itime-1)*h


        #
        # solution increment Dy and solution time derivative ydot_np1 at time n+1 are initially set to 0.
        #

        ydot_np1[:] .= 0.0
        # ydot_np1[iloc_v] .= 0.0
        # ydot_np1[iloc_xi] = ydot_n[iloc_xi]

        # Dy .= 0.0

        #
        # computation of Wy and Dy predictors
        #

        Dy[iloc_xiv]=h*(1.0-theta)*wy[iloc_xiv]
        if Ndofs_mult > 0
            Dy[iloc_mult] .= 0.
        end
        wy[iloc_xiv]=1.0/(1.0-delta_m)*(delta_f*ydot_n[iloc_xiv] - delta_m*wy[iloc_xiv])
        
		Dy[iloc_xiv] += h*theta*wy[iloc_xiv]
        
        #
        #
        # Arbitrary inittialization of test and TOL to force execution of first iteration
        #
        


        global test = 1.e7
        global niter = 0
        global TOL = PREC
        tol_sparse = eps(Float64)*10.0

        #
        # Nodes current position and orientation is incremented by Dy at nodes
        #

        # increment_frames(Dy,Nnodes)

        #
        # iteration on equilibrium, velocity and contraint equations
        #

        # global p_init = zeros(Ndofs)
        # global p_av_init = 0.5*p

        while (test > TOL && niter < NitMax)

            #
            # iteration matrix, residual vector and load vector are reinitialized.
            # iteration counter is set to 0
            #

            niter += 1

            res[:] .= 0.0

            ext_work = 0.0

            #
            # loop on elements
            #

            t_element_start = my_timer()

            global update = false 
            (matrix_update == true ||  niter == 1) && (update = true)

            if update == true
                kin_energy, pot_energy, str_energy = dynamic_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,theta_p,ext_work,itime,h,niter)
            else
                kin_energy, pot_energy, str_energy = dynamic_element_forces(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter)
            end

            t_element_stop = my_timer()

            time_elements += t_element_stop - t_element_start

            #
            # convergence tolerance TOL is computed at first iteration
            #

            if niter == 1
                TOL =  norm(p) + norm(res[iloc_x] - p[iloc_x]) + norm(res[iloc_v]) + norm(res[iloc_mult])
                Ndofs_int > 0  && (TOL += norm(res[iloc_int]))
                TOL = (TOL+1.0)*PREC
            end

            # println("itime, niter, TOL, test ", itime, " ", niter, " ",TOL, " ", test)

            #
            # add load to residual vector
            #

            test=norm(res[iloc_x]) + norm(res[iloc_v])

            Ndofs_int > 0 && (test += norm(res[iloc_int]))
            Ndofs_mult > 0 && (test += norm(res[iloc_mult]))

 
            update == true  && (lu!(F,sparse(I_qq,J_qq,sps.S_qq));S_vq = sparse(I_vq,J_vq,sps.S_vq))

            dy[iloc_q] = F\res[iloc_q]

            dy[iloc_v] = S_vq[:,sps.subset_q]*dy[sps.subset_q] - res[iloc_v]
            Dy += dy
            ydot_np1[iloc_xiv] += theta_p*dy[iloc_xiv]

            # correct_frames(dy,theta_p,Nnodes)

        end

        sum_ext_work += ext_work
        
        if energy_balance == true
        
            kin_energy, pot_energy, str_energy = dynamic_element_forces(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter)
        end

        if eigvals == true &&  (mod(itime,eig_freq) == 0)

            vals = eig_solve(y_n,Dy,ydot_np1,max_vals,itime,h,niter)

        end

        #
        # Solution update for pseudo-velocity wy, velocities and multipliers
        #

        wy[iloc_xiv] += (1.0-delta_f)/(1.0-delta_m)*ydot_np1[iloc_xiv]


        y_n[iloc_iv]  += Dy[iloc_iv]

        Ndofs_mult > 0 &&  (y_n[iloc_mult] = Dy[iloc_mult])


        #
        # solution update at nodal frames
        #

        increment_frames(Dy,y_n,Nnodes)

        #
        # shift in time solution time derivative
        #

        ydot_n[:] = ydot_np1[:]

        #
        # add solution step to solution database
        #


        if save == true 

            itime_saved, itime_vals_saved = record_on_h5_file_dynamic(dsets,itime,itime_saved,itime_vals_saved,times,niter,y_n,ydot_np1,p,
            pot_energy,kin_energy,str_energy,sum_ext_work,Npas,save,save_freq,eigvals,eig_freq, vals)

            

        end

        verbose == true && (println("itime, niter, TOL, test ", itime, " ", niter, " ",TOL, " ", test))

        niter == NitMax && (println("NitMax: ", NitMax, " reached at step ", itime, " time ",time," TOL ",TOL, " test ",test))

        nit_tot += niter

    end

    t_stop = my_timer()

    println("elapsed time: ", t_stop - t_start, " milliseconds")
    println("time in elements: ", time_elements , " milliseconds")


    nit_mean = nit_tot/Npas

    println("average number of iterations: ", nit_mean)

    save == true && (close(file))


    return h5_file
end
