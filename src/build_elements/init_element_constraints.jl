
"""
    init_element_constraints
"""
function init_element_constraints(Nel::Int64)

    ec =  Main.element_container
    mc = Main.model_container
    inc = Main.SetElements.inequality_container

    init_file = mc.init_file

    global idx_upper_bounds = Vector{Int}[]
    global idx_lower_bounds = Vector{Int}[]

    Nbounds = mc.Inequalities
    init_file = mc.init_file
    
    nbound =0
   
    for iel = 1:Nel
        el_type = ec.element_types[iel]
        if el_type == "inequality" 
            nbr = ec.element_numbers[iel]
            iel2 = findfirst(x -> x == nbr, inc.number)[1]
            bound_type = inc.bound_type[iel2]
            nbound += 1
            bound_type == "upper" && (idx_upper_bounds = [idx_upper_bounds; nbound])
            bound_type == "lower" && (idx_lower_bounds = [idx_lower_bounds; nbound])
        end

    end

    if init_file == ""
        global lambda_n = zeros(Nbounds)
        global bounds = zeros(Nbounds)
    else
        h5_file = init_file*".h5"
        Lambda, Bounds = read_results_constraints(h5_file)
        step_n = size(Lambda,1)
        global lambda_n = Lambda[step_n,:]
        global bounds =  Bounds[step_n,:]
    end


    return  idx_upper_bounds, idx_lower_bounds, lambda_n, bounds

end