"""
    assemble
"""
function assemble(sol_type::String,uniform_rotation::Bool)

    mc=Main.model_container

    if sol_type != "static" && sol_type != "dynamic" && sol_type != "static_constrained"
        error(" assemble : incorrect solution type ", sol_type)
    end

    struc_loc_q, struc_loc_v = get_struc_loc(sol_type,uniform_rotation)

    mc.struc_loc_q=[mc.struc_loc_q;struc_loc_q]
    mc.struc_loc_v=[mc.struc_loc_v;struc_loc_v]

    mc.end_of_assembly=true

    return

end

function assemble()
    sol_type = "dynamic"
    uniform_rotation = false
    assemble(sol_type,uniform_rotation)
end

function assemble(sol_type::String)
    uniform_rotation = false
    assemble(sol_type,uniform_rotation)
end