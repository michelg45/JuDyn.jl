"""
    PrismaticJointArray

Data structure for the `SetElements.prismatic_joint_container` array, created by the `set_prismatic_joint` function. It contains the following data for each prismatic joint element of the element set:

|   |  |
|:-----------------------------------|:---------------------------------------------|
| number::Vector{Int} | element numbers |
| node_orders::Vector{Vector{Int}} | set of  node pairs connected by the prismatic joint element. |
| l_0::Vector{Float64} | vector of initial lengths |
| RV_rel::Vector{RV3} | vector of relative orientations between end nodes |
| RV_line::Vector{RV3} | vector of frames with ``x`` axis aligned on prismatic joint axis. |
| scale_factor::Vector{Float64} | set of constraint scale factors. |
| mode::Array{String,1} | set of strings  take thing values "free", "force", "driven", "spring". |
| time_function::Vector{String}| time function  strings from `ÌnputFunctions` | d
| params::Vector{Vector{Float64}} | sets of parameters defining either the time functions  of "force" and "driven" modes  or  the viscoelastic spring parameters ("spring" mode). | 

"""
mutable struct PrismaticJointArray

    number::Vector{Int}
    node_orders::Vector{Vector{Int}}
    l_0::Vector{Float64}
    RV_rel::Vector{RV3}
    RV_0::Vector{RV3}
    RV_line::Vector{RV3}
    scale_factor::Vector{Float64}
    mode::Array{String,1}
    time_function::Vector{String}
    params::Vector{Vector{Float64}}

    function PrismaticJointArray()

        number = []
        node_orders = []
        l_0 = []
        RV_rel = []
        RV_0 = []
        RV_line = []
        scale_factor = []
        mode = []
        time_function = []
        params = []

        return new(number, node_orders,l_0,RV_rel,RV_0,RV_line,scale_factor, mode, time_function, params)

    end


end

prism_container = PrismaticJointArray()