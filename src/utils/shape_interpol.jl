function shape_interpol(node_positions::Array,t::Float64, T::Float64,Npas::Int)
    i1 = lower_lim(t,T,Npas)
    h = T/(Npas-1)
    t1 = (i1 - 1)*h
    if t < T
        inter_shape = node_positions[i1,:,:] + 
            (node_positions[i1+1,:,:] - node_positions[i1,:,:])*(t-t1)/h
    else
        inter_shape = node_positions[Npas,:,:]
    end
    return inter_shape
end