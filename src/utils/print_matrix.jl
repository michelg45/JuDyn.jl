macro Name(arg)
    string(arg)
end


function print_matrix(S::Array{Float64,2})
    println()
    a = @Name(S)
    println(a)
    println()
    show(stdout, "text/plain",S)
    println()
    return
end
