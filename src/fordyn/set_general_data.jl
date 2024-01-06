"""
    set_general_data

        function allowing to set the following data / parameters for the model: Name::String, gravity::Vec3, 
        uniform_rotation::Bool, rotation_speed::Vec3, matrix_update::Bool.

        the function "set_general_data" has different instances:

            - set_general_data(Name,gravity,uniform_rotation,rotation_speed,matrix_update)
                all parameters need to be defined.
            - set_general_data(Name,gravity,uniform_rotation,rotation_speed)
                all parameters except "matrix_update" need to be defined.
            - set_general_data(Name,gravity) 
                "uniform_rotation" is set to false, "rotation_speed"  to Vec3(). 
            -  set_general_data(Name,gravity,matrix_update) 
                idem.
"""
function set_general_data(Name::String,gravity::Vec3,uniform_rotation::Bool,rotation_speed::Vec3)
    mc = Main.model_container
    mc.Name = Name
    mc.gravity = gravity
    mc.uniform_rotation = uniform_rotation
    mc.rotation_speed = rotation_speed
    return
end

function set_general_data(Name::String,gravity::Vec3)
    mc = JuDyn.FORDYN.model_container
    mc.Name = Name
    mc.gravity = gravity
    mc.uniform_rotation = false
    mc.rotation_speed = Vec3(0.0,0.0,0.0)
    return
end

function set_general_data(Name::String,gravity::Vec3,matrix_update::Bool)
    mc = Main.model_container
    mc.Name = Name
    mc.gravity = gravity
    mc.matrix_update = matrix_update
    mc.uniform_rotation = false
    mc.rotation_speed = Vec3(0.0,0.0,0.0)
    return
end

function set_general_data(Name::String,gravity::Vec3,uniform_rotation::Bool,rotation_speed::Vec3,matrix_update::Bool)
    mc = Main.model_container
    mc.Name = Name
    mc.gravity = gravity
    mc.uniform_rotation = uniform_rotation
    mc.rotation_speed = rotation_speed
    mc.matrix_update = matrix_update
    return
end