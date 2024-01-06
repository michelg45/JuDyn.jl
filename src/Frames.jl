__precompile__()

"""
    Frames

        Module for  management of nodal frame configuration in array "current_frames".
        It performs the following operations through call to functions:

        init_frames:        frame intialization at the beginning of the solution phase.
        increment_frames :  frame incrementation at the beginning of a time step

"""
module Frames

using LinearAlgebra
using VMLS
using ..MyAlgebra

dir = "frames/"


include(dir*"CurrentFrame.jl")
include(dir*"init_frames.jl")
include(dir*"increment_frames.jl")






export CurrentFrame
export init_frames
export increment_frames


end # module