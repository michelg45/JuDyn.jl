"""
    solve

    function activating the different solution methods available in the JUDYN library.

    3 solution methods are currently available:

        "static" :              step-by-step Newton-Raphson iteration.
        "dynamic" :             dynamic response using the first-order generalized-alpha method.
        "static_constrained" :  step-by-step Newton-Raphson iteration under inequality constraints. 

        Input:
            sol_type::String        "static", "dynamic" or "static_constrained".
            JSON_file::String       JSON file specifying solution parameters.
            uniform_rotation::Bool  uniform_rotation = true allows computing the static or dynamic 
                                    response in a uniformly rotating reference frame.

        Calling sequences:

            h5_file = solve(JSON_file,sol_type,uniform_rotation)
            h5_file = solve(JSON_file) => sol_type => "dynamic", uniform_rotation = false
            h5_file = solve(JSON_file,sol_type) => sol_type = uniform_rotation = false

"""
function solve(JSON_file::String)
    global sparse_matrix_S =  sparse_matrix_prepro()

    println("solution method: generalized-alpha")
    # return h5_file, mat_file = dynamic_solve(JSON_file)
    return h5_file = dynamic_solve(JSON_file)
end

function solve(JSON_file::String,sol_type::String,uniform_rotation::Bool)
    global sparse_matrix_S =  sparse_matrix_prepro()
    if sol_type == "static"
        println("solution method: static")
        return h5_file = static_solve(JSON_file,uniform_rotation)
    elseif sol_type  == "dynamic"
        println("solution method: generalized-alpha")
        return h5_file = dynamic_solve(JSON_file)
    elseif sol_type  == "static_constrained"
        return h5_file = static_constrained_solve(JSON_file,uniform_rotation)    
    else
        error("unknown solution method ")
    end
end

function solve(JSON_file::String,sol_type::String)
    uniform_rotation = false
    solve(JSON_file,sol_type,uniform_rotation)
end