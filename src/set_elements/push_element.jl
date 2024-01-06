function push_element(S::Array,res::Vector,p::Vector,loc::Array{Int,1},S_el::Array,res_el::Vector,p_el::Vector)

    loc_el = nonnuls((loc .> 0).*[i for i in 1:size(loc,1)])
    loc = nonnuls(loc)

    S[loc,loc] += S_el[loc_el,loc_el]
    res[loc]   += res_el[loc_el]
    p[loc]     += p_el[loc_el]

    return
end

function push_element(S::Array,res::Vector,loc::Array{Int,1},S_el::Array,res_el::Vector)

    loc_el = nonnuls((loc .> 0).*[i for i in 1:size(loc,1)])
    loc = nonnuls(loc)

    S[loc,loc]   += S_el[loc_el,loc_el]
    res[loc]     += res_el[loc_el]

    return
end

function push_element(S::Array,loc1::Array{Int,1},loc2::Array{Int,1},S_el::Array{Float64,2})

    loc_el_1 = nonnuls((loc1 .> 0).*[i for i in 1:size(loc1,1)])
    loc_el_2 = nonnuls((loc2 .> 0).*[i for i in 1:size(loc2,1)])
    loc1 = nonnuls(loc1)
    loc2 = nonnuls(loc1)

    S[loc1,loc2] += S_el[loc_el_1,loc_el_2]

    return
end

