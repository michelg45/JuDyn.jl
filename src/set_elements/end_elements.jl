function end_elements()
    ec = Main.element_container
    mc = Main.model_container
    shift_int = copy(mc.max_x)
    shift_mult = shift_int + mc.max_int
    shift_v =  shift_mult + mc.max_mult
    nel = mc.Elements
    for i = 1:nel
        if ec.n_int[i] > 0
            ec.loc_int[i] =  ec.loc_int[i]  .+ shift_int
            ec.inv_loc_int[i] =  ec.inv_loc_int[i]  .+ shift_int
        end

        if ec.n_mult[i] > 0
            ec.loc_mult[i]  =  ec.loc_mult[i]  .+ shift_mult
            ec.inv_loc_mult[i]  =  ec.inv_loc_mult[i]  .+ shift_mult
        end
    end
    mc.end_of_elements=true
    temp = reduce(vcat, ec.loc_int)
    size(temp,1) > 0 ? mc.max_v = findmax(temp)[1] : mc.max_v = 0
    temp = reduce(vcat, ec.loc_v)
    size(temp,1) > 0 ? mc.max_v = findmax(temp)[1] : mc.max_v = 0
    temp = reduce(vcat, ec.loc_mult)
    size(temp,1) > 0 ? mc.max_mult = findmax(temp)[1] : mc.max_mult = 0

    return
end