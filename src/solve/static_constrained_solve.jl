"""
    static_constrained_solve

        function which computes in several time increments on a [0.0, T] time interval 
        the static nonlinear solution subject of inequality constraints.

        The inequality constraints are defined as elements using the 'set_inequality' function. 
        They are solved in an inner loop to the equilibrium loop using the 'Ipopt' algorithm ('Jump' library).
            
        The time interval on which the forcing functions (e.g.-. 'set_node_force') and the inequality
        constraints must  match the time internal [0.0, T] of the quasi-static solution.
        

        the time variation of the forcing functions should be set to "linear"
        and the time internal T set to 1.0.

        The parameters of the solution algorithm are read from the file JSON_file.

            T = time internal of the analysis
            k = scaling factor
            PREC = precision parameter
            PREC_OPT = precision parameter on the nonlinear constraints  ('Ipopt').
            rho_inf = spectral radius if the the integration algorithm.
            NitMax = maximum number of iterations
            Npas = number of time steps
            verbose = true: print at end of each time step.
            verbose_constraints = true: print of 'Ipopt' summary model.

            problem_name = name definition of the HDF5 file.
            save = true: output to HDF5 file.
            energy_balance = true: re-evaluation of energies at the end of the time step.

        Output results are saved on "name.h5" file using the "record_on_h5_file_static.jl" function.

"""
function static_constrained_solve(JSON_file::String,uniform_rotation::Bool)


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
    rho_inf = dict["Time_Integrator"]["parameters"]["rho_inf"]
    haskey(dict["Time_Integrator"]["parameters"],"PREC_OPT") == true ? (PREC_OPT = dict["Time_Integrator"]["parameters"]["PREC_OPT"]) : PREC_OPT = 1.0e-12
    NitMax = dict["Time_Integrator"]["parameters"]["NitMax"]
    Npas = dict["Time_Integrator"]["parameters"]["Npas"]
    verbose = dict["Time_Integrator"]["parameters"]["verbose"]
    haskey(dict["Time_Integrator"]["parameters"],"verbose_constraints") == true ? (verbose_constraints = dict["Time_Integrator"]["parameters"]["verbose_constraints"]) : verbose_constraints = false

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
    energy_balance = dict["Output"]["parameters"]["energy_balance"]

    # external_force = Main.InputFunctions.input_functions[dict["Input"]["external_force"]]

    #
    # computing solution parameters
    #


    h = T/(Npas-1)

    delta_m=0.5*(3.0*rho_inf-1.0)/(rho_inf+1.0)
    delta_f=rho_inf/(rho_inf+1.0)
    theta=0.5 + delta_f- delta_m
    global theta_p=(1.0-delta_m)/(h*theta*(1.0-delta_f))


   

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
    Ndofs_v= model_container.Ndofs_v
    uniform_rotation == true && (Ndofs_v = model_container.Ndofs_v)
    Ndofs_int= model_container.Ndofs_int
    Ndofs_mult= model_container.Ndofs_mult
    Nel = model_container.Elements
    matrix_update = model_container.matrix_update
    element_numbers = element_container.element_numbers
    element_types = element_container.element_types
    Nbounds = model_container.Inequalities
    Nlambda = Nbounds
    initial_shape = model_container.initial_shape
    init_file = model_container.init_file
    mc = model_container
    mc.Shells > 0 && (sc = Main.SetElements.shell_container)
    mc.Beams > 0 && (bc = Main.SetElements.beam_container)
    
    Ndofs_v > 0  && (error("static constrained solution: Ndofs_v must be 0 "))

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
    global wy = Vector{Float64}(zeros(Ndofs))
    


    y_n[:] = Main.InitialConditions.x_0[:]
    Dy[:] = Main.InitialConditions.Dx_0[:]

   

    if Nbounds > 0

        idx_upper_bounds, idx_lower_bounds, lambda_n, bounds = init_element_constraints(Nel)

        println("lambda_n ", lambda_n)
        println("bounds ", bounds)

        n_upper_bounds = size(idx_upper_bounds,1)
        n_lower_bounds = size(idx_lower_bounds,1)

    else
        lambda_n = zeros(1)
        bounds = zeros(1)

    end

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

    cf = Main.Frames.current_frames

        #
    # computation of initial energies
    #


    #
    # if save == true build HDF5 solution database
    #

    itime = 1
    niter = 1

    res[:] .= 0.0

    global ext_work = 0.0

    global itime_vals_saved = 1

    pot_energy, str_energy = static_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_n,res,p,ext_work,itime,h,niter,theta_p)

    increment_frames(Dy,y_n,Nnodes)

    Dy[:] .= 0.0



    if save == true
        file, dsets = create_h5_file_static(h5_file,Npas,eigvals,eig_freq,max_vals)
    
        vals = zeros(max_vals)
        times = (itime -1)*h
        niter = 0

        record_on_h5_file_static(dsets,itime,itime_vals_saved,times,niter,y_n,ydot_n,p,
            pot_energy,str_energy,ext_work,eigvals,eig_freq, vals,lambda_n,bounds)

    end

    #
    # pseudo-velocity wy is initially set to 0.
    #
    #
    # timer initialization
    #

        t_start = (Dates.hour(now())*3600+Dates.minute(now())*60+Dates.second(now()))*1000+Dates.millisecond(now())

        println("start constrained static solve")

        Npas == 1 && (Npas = 2)

        nit_tot = 0

        println("Npas, NitMax, eigvals ", Npas, NitMax, eigvals)

        # Nbounds > 0 && (global lambda_n = zeros(Nlambda))

    Ndofs_int > 0  && (wy[iloc_int] = ydot_n[iloc_int])
    
    

    for itime = 2:Npas

        times = (itime-1)*h

        ydot_np1[iloc_int] .= 0.0

            
        test = 1.e7

        niter = 0

        TOL = PREC

        PRC = 1.e-8

        res .= 0.0

        if Ndofs_int  > 0

            Dy[iloc_int]=h*(1.0-theta)*wy[iloc_int]
            wy[iloc_int]=1.0/(1.0-delta_m)*(delta_f*ydot_n[iloc_int] - delta_m*wy[iloc_int])
		    Dy[iloc_int] += h*theta*wy[iloc_int] 
            
        end
               
        while (test > TOL && niter < NitMax || niter < 2)

            #
            # iteration matrix, residual vector and load vector are reinitialized.
            # iteration counter is set to 0
            #
           

            niter +=1

            p  .= 0.
            res .= 0.

            ext_work = 0.

            global update = false 
            (matrix_update == true ||  niter == 1) && (update = true)

            #
            # loop on elements
            #
            if update == true
                pot_energy, str_energy = static_element_system(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter,theta_p)
            else
                pot_energy, str_energy = static_element_forces(Nel,element_numbers,element_types, y_n,Dy,ydot_np1,res,p,ext_work,itime,h,niter)
            end

            if Nbounds > 0 
                B, x_n, bounds  = element_constraints(Ndofs,Nel,Nbounds,bounds,y_n, Dy,itime,h)
            end   

    
            # convergence tolerance TOL is computed at first iteration
            #

            if niter == 1
                TOL =  (norm(p)+ norm(res[iloc_x] - p[iloc_x]) + norm(res[iloc_mult]))/Ndofs
                uniform_rotation == true && (TOL += norm(res[iloc_v]))
                Ndofs_int > 0 && (TOL += norm(res[iloc_int]))
                TOL = (TOL+1.0)*PREC
            end

            #
            # add load to residual vector
            #
            if Nbounds > 0  && niter > 1 (test=norm((res+B*lambda_n)[iloc_x])/Ndofs_x + norm(dy)/Ndofs_x +norm(lambda_n.*(x_n-bounds))/Nbounds)
            # if Nbounds > 0  && niter > 1 (test=norm((res+B*lambda_n)[iloc_x]))    
            else 
                test=(norm(res[iloc_x]))/Ndofs 
            end

            Ndofs_int > 0 && (test += norm(res[iloc_int]))
            Ndofs_mult > 0 && (test += norm(res[iloc_mult]))
            uniform_rotation == true && (test += norm(res[iloc_v]))

            
            update == true  && (lu!(F,sparse(I_qq,J_qq,sps.S_qq)))

            if Nbounds > 0 
                rhs = [res[iloc_q] B]
            else 
                rhs = res[iloc_q]
            end

             sol = F\rhs

             if Nbounds > 0 

                G = sol[:, 2:(Nlambda+1)]
                A = B'*G
                u = B'*sol[:, 1]

                normG = zeros(Nlambda)
                for i = 1:Nlambda
                    normG[i] = norm(G[:,i])
                end
                
                norm(lambda_n) == 0.0 && (lambda_n = A\(bounds-u - x_n))

                model = Model(Ipopt.Optimizer)

                verbose_constraints == false && set_silent(model)
                @variable(model, lambda[i=1:Nlambda],start=lambda_n[i])
                n_lower_bounds > 0 && ( @constraint(model,  (u +  A*lambda + x_n - bounds)[idx_lower_bounds]  .>= 0.0))
                n_lower_bounds > 0 && (@constraint(model,lambda[idx_lower_bounds] .>= 0.0))
                n_upper_bounds > 0 && (@constraint(model,  (u +  A*lambda + x_n - bounds)[idx_upper_bounds] .<= 0.0 ))
                n_upper_bounds > 0 && (@constraint(model,lambda[idx_upper_bounds] .<= 0.0))
                # @NLconstraint(model, [i=1:Nlambda], (lambda[i]*(u +  A*lambda + x_n - bounds)[i]) == 0.0 )
                @NLconstraint(model, [i=1:Nlambda], -PREC_OPT <= (lambda[i]*(u +  A*lambda + x_n - bounds)[i]) <= PREC_OPT)
                optimize!(model)
                verbose_constraints == true && solution_summary(model)


                lambda_n = value.(lambda)

"""                println("itime ",itime, " niter ", niter)
                println(" lambda_n ", lambda_n)
                println("A ", A)
                println("bounds ", bounds)
                println("u ", u)
                println("x_n ", x_n)"""

            end
            
            Nbounds > 0 ?  dy = sol[:, 1] + G*lambda_n :  dy = sol[:, 1]
         
 
            Nbounds > 0 && niter == 1 && (TOL +=norm(B*lambda_n)*PREC; PRC = PREC*norm(bounds)*norm(lambda_n))


            uniform_rotation == true && (S_vq = sparse(I_vq,J_vq,sps.S_vq); 
                dy[iloc_v] = S_vq[:,sps.subset_q]*dy[sps.subset_q] - res[iloc_v])

            Dy += dy

            Ndofs_int > 0 && (ydot_np1[iloc_int] += theta_p*dy[iloc_int])

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

        if Ndofs_int > 0
            wy[iloc_int] += (1.0-delta_f)/(1.0-delta_m)*ydot_np1[iloc_int]
            y_n[iloc_int]  += Dy[iloc_int]
            ydot_n[iloc_int] = ydot_np1[iloc_int]
        end

        Ndofs_mult > 0 && (y_n[iloc_mult] = Dy[iloc_mult])

        uniform_rotation == true && (y_n[iloc_v]  += Dy[iloc_v])
      
        increment_frames(Dy,y_n,Nnodes)


        if save == true

            record_on_h5_file_static(dsets,itime,itime_vals_saved,times,niter,y_n,ydot_np1,p,
            pot_energy,str_energy,ext_work,eigvals,eig_freq, vals,lambda_n,bounds)

        end

     #   verbose == true && (println("itime, niter, TOL, test ", itime, " ", niter, " ",TOL, " ", test))
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
