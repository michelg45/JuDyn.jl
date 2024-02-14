function static_solve(JSON_file::String,uniform_rotation::Bool)

    """
        the function "static_solve" computes the static nonlinear solution
        of the system in several increments.

        the time variation of the forcing functions should be set to "linear"
        and the time internal T set to 1.0.

        Its parameters are read from the file JSON_file.

        T = time internal of the analysis
        k = scaling factor
        PREC = precision parameter
        NitMax = maximum number of iterations
        Npas = number of time steps
        elimv = true: elimination of velocities before solution of correction equation.
        verbose = true: print at end of each time step.

        problem_name = name definition of the HDF5 file.
        save = true: output to HDF5 file.
        save_mat = true: output to MAT file.
        energy_balance = true: re-evaluation of energies at the end of the time step.


    """



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
    save_mat = dict["Output"]["parameters"]["save_mat"]
    energy_balance = dict["Output"]["parameters"]["energy_balance"]

    # external_force = Main.InputFunctions.input_functions[dict["Input"]["external_force"]]

    #
    # computing solution parameters
    #


    h = T/(Npas-1)


    theta_p= 0.0

    #
    # output files. mat_file can be crated only if save = true
    #

    sh5 = ".h5"
    smat = ".mat"

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
    Ndofs_q = model_container.Ndofs_q
    uniform_rotation == true && (Ndofs_v = model_container.Ndofs_v)
    Ndofs_int= model_container.Ndofs_int
    Ndofs_mult= model_container.Ndofs_mult
    Nel = model_container.Elements
    matrix_update = model_container.matrix_update
    initial_shape = model_container.initial_shape
    element_numbers = element_container.element_numbers
    element_types = element_container.element_types

    #
    # Initialization ov vectors and arrays
    #


    global p = Vector{Float64}(zeros(Ndofs))
    global res = Vector{Float64}(zeros(Ndofs))
    global y_n = Vector{Float64}(zeros(Ndofs))
    global Dy = Vector{Float64}(zeros(Ndofs))
    global dy = Vector{Float64}(zeros(Ndofs))
    global ydot_n = Vector{Float64}(zeros(Ndofs))
    global ydot_np1 = Vector{Float64}(zeros(Ndofs))

    y_n[:] = Main.InitialConditions.x_0[:]
    Dy[:] = Main.InitialConditions.Dx_0[:]

 
    sps = sparse_matrix_S

    I_qq = sps.I_index_qq
    J_qq = sps.J_index_qq


    nz_qq = size(I_qq,1)

    uniform_rotation == true && (I_vq = sps.I_index_vq; J_vq = sps.J_index_vq)
#
#   pre-factorization of system matrix
#
global F = lu(sparse(I_qq,J_qq,rand(nz_qq)))

    #
    # Ordering of degrees of freedom for the global system
    # _x :      inematic variables at node frames
    # _int:     internal variables to elements
    # _mult:    multipliers
    # _v:       velocities (defined at element level)
    #

    iloc_x=[i for i in 1:Ndofs_x]
    Ndofs_int > 0 ? iloc_int=[i for i in Ndofs_x+1:Ndofs_x+Ndofs_int] :  iloc_int=[]
    Ndofs_mult > 0 ? iloc_mult=[i for i in Ndofs_x+Ndofs_int+1:Ndofs_x+Ndofs_int+Ndofs_mult] :  iloc_mult=[]
    iloc_q = [iloc_x;iloc_int;iloc_mult]
    uniform_rotation == true && (iloc_v=[i for i in Ndofs_q+1:Ndofs_q+Ndofs_v]; 
        iloc_xiv = [iloc_x; iloc_int; iloc_v]; iloc_iv = [iloc_int; iloc_v])

    iloc_xim = [iloc_x;iloc_int;iloc_mult]
    iloc_xi = [iloc_x;iloc_int]

    init_frames(y_n,ydot_n,initial_shape)

    #
    # computation of initial energies
    #


    #
    # if save == true build HDF5 solution database
    #

    lambda_n = zeros(1)
    bounds = zeros(1)
    
    global itime_vals_saved = 1

    itime = 1

    global ext_work = 0.0

    niter = 1

    pot_energy, str_energy = static_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_n,res,p,ext_work,itime,h,niter,theta_p)
    
    increment_frames(Dy,y_n,Nnodes)

    if save == true

        file, dsets = create_h5_file_static(h5_file,Npas,eigvals,eig_freq,max_vals)

        

        itime = 1

        vals = zeros(max_vals)
        times = (itime -1)*h
        niter = 0
        pot_energy = 0.0
        str_energy = 0.0
        ext_work = 0.0



        record_on_h5_file_static(dsets,itime,itime_vals_saved,times,niter,y_n,ydot_np1,p,
            pot_energy,str_energy,ext_work,eigvals,eig_freq, vals,lambda_n,bounds)


    end

    #
    # pseudo-velocity wy is initially set to 0.
    #
    #
    # timer initialization
    #

        t_start = (Dates.hour(now())*3600+Dates.minute(now())*60+Dates.second(now()))*1000+Dates.millisecond(now())

        println("start static solve")

        Npas == 1 && (Npas = 2)

        nit_tot = 0

        println("Npas, NitMax, eigvals ", Npas, NitMax, eigvals)

    for itime = 2:Npas

        times = (itime-1)*h


        Dy .= 0.

        test = 1.e7
        global niter = 0
        TOL = PREC

        
        while (test > TOL && niter < NitMax)

            global update = false 
            (matrix_update == true ||  niter == 1) && (update = true)    

            #
            # iteration matrix, residual vector and load vector are reinitialized.
            # iteration counter is set to 0
            #

            niter +=1

            p  .= 0.
            res .= 0.

            sps.S_qq .= 0.0

            uniform_rotation == true && (sps.S_vq .= 0.0)

            ext_work = 0.

            global update = false 
            (matrix_update == true ||  niter <= 2) && (update = true)

            #
            # loop on elements
            #
            if update == true
                pot_energy, str_energy = static_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter,theta_p)
            else
                pot_energy, str_energy = static_element_forces(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter)
            end

    
            # convergence tolerance TOL is computed at first iteration
            #

            if niter == 1

                TOL =  norm(p)+ norm(res[iloc_x]) + norm(res[iloc_mult])
                uniform_rotation == true && (TOL += norm(res[iloc_v]))
                Ndofs_int > 0 && (TOL += norm(res[iloc_int]))
                TOL = (TOL+1.0)*PREC
            end

            #
            # add load to residual vector
            #
            
            test=norm(res[iloc_x])
            Ndofs_int > 0 && (test += norm(res[iloc_int]))
            Ndofs_mult > 0 && (test += norm(res[iloc_mult]))
            uniform_rotation == true && (test += norm(res[iloc_v]))

           
            update == true  && (lu!(F,sparse(I_qq,J_qq,sps.S_qq)))
            
            (update == true && uniform_rotation == true) 

            

            dy[iloc_q] = F\res[iloc_q]

            

            (update == true && uniform_rotation == true)  && (S_vq = sparse(I_vq,J_vq,sps.S_vq); 
                dy[iloc_v] = S_vq[:,sps.subset_q]*dy[sps.subset_q] - res[iloc_v])

            Dy += dy

            # correct_frames(dy,theta_p,Nnodes)

        end

        verbose == true && (println("convergence ", "itime ", itime, " niter ", niter, " test ", test, " TOL ", TOL))
  


        if energy_balance == true

            ext_work = 0.

            pot_energy, str_energy = static_element_forces(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter)

        end

        if eigvals == true &&  (mod(itime,eig_freq) == 0) 

            vals = eig_solve(y_n,Dy,ydot_np1,max_vals,itime,h,niter)
            itime_vals_saved += 1

        end

     
        Ndofs_int > 0 && ( y_n[iloc_int]  = Dy[iloc_int])
        Ndofs_mult > 0 && (y_n[iloc_mult] = Dy[iloc_mult])

        uniform_rotation == true && (y_n[iloc_v]  += Dy[iloc_v])
      
        increment_frames(Dy,y_n,Nnodes)
      
        if save == true

            record_on_h5_file_static(dsets,itime,itime_vals_saved,times,niter,y_n,ydot_np1,p,
            pot_energy,str_energy,ext_work,eigvals,eig_freq, vals,lambda_n,bounds)
            
        end

        verbose == true && (println("itime, niter, TOL, test ", itime, " ", niter, " ",TOL, " ", test))
        niter == NitMax && (println("NitMax: ", NitMax, " reached at step ", itime, " time ",times," TOL ",TOL, " test ",test))
        nit_tot += niter
    end

    t_stop = (Dates.hour(now())*3600+Dates.minute(now())*60+Dates.second(now()))*1000+Dates.millisecond(now())
    println("elapsed time: ", t_stop - t_start, " milliseconds")


    nit_mean = nit_tot/Npas

    println("average number of iterations: ", nit_mean)

    save == true && (close(file))

    return h5_file
end
