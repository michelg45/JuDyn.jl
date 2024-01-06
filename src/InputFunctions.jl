__precompile__()

"""
    InputFunctions

        This mModule is currently a mess.
        It contains the definition of different input functions to define
        - external loading on a specified dof (called by "Solve.jl")
        - internal loading, thus called by elements (e.g. ground hinge)

        The input functions are transmitted through a dictionary ("Dict" structure)

"""
module InputFunctions
using ..Nodes
input_functions = Dict()

include("input_functions/step.jl")
input_functions["step"] = step
export step

include("input_functions/inverse_step.jl")
input_functions["inverse_step"] = inverse_step
export inverse_step

include("input_functions/double_step.jl")
input_functions["double_step"] = double_step
export double_step

include("input_functions/triangle_impulse.jl")
input_functions["triangle_impulse"] = triangle_impulse
export triangle_impulse

include("input_functions/linear.jl")
input_functions["linear"] = linear
export linear

include("input_functions/linear_velocity.jl")
input_functions["linear_velocity"] = linear_velocity
export linear_velocity

include("input_functions/crank_motion.jl")
input_functions["crank_motion"] = crank_motion
export crank_motion

include("input_functions/null_torque.jl")
input_functions["null_torque"] = null_torque
export null_torque

include("input_functions/haug_rotation_profile.jl")
input_functions["haug_rotation_profile"] = haug_rotation_profile
export haug_rotation_profile

include("input_functions/shaft_velocity_profile.jl")
input_functions["shaft_velocity_profile"] = shaft_velocity_profile
export shaft_velocity_profile

include("input_functions/shaft_velocity_profile2.jl")
input_functions["shaft_velocity_profile2"] = shaft_velocity_profile2
export shaft_velocity_profile2

include("input_functions/piecewise.jl")
input_functions["piecewise"] = piecewise
export piecewise

export input_functions

end # module
