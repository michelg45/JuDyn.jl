"""
    find_positive_integers

        Function finding the list and number of positive integers in a vector of integers

        Input : v::Array{Int}
        Output : (n, w) number and list of positive integers

        Calling sequence: (n,w) = find_positive_integers(v)

"""
function find_positive_integers(v::Array{Int})


    cartindx = findall(x -> x > 0, v)
    n = size(cartindx,1)
    index = [cartindx[i] for i in 1:n]
    w = v[index]
    return n, w
    end

function end_BC()
    return
end