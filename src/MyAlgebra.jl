__precompile__()

"""
    MyAlgebra

"MyAlgebra.jl" is a module  introducing different algebraic user types and the associated operations.

User Types are noted with capital letters (_Vec3_, _Mat3_, _RV3_, _Quat_, _NodeFrame_)  and ordinary functions start with lower case letters.  Constructors provide type conversion.

The user types are:
>
>  _Vec3_ :    vector with 3 components (e.g. position, velocity, Euler angles). 
> 
> _Mat3_ :    3 x 3 matrix (e.g. rotation and tangent operators).
>
> _RV3_ :  rotation vector `` \\mathbf{n}  \\phi ``.
>
> _Quat_ :  quaternion  `` \\mathbf{q} = [q_0,  q_1, q_2, q_3] `` .
>
> _NodeFrame_ :   compound structural type to describe a nodal frame `` \\mathbf{H} (\\mathbf{x}, \\mathbf{p}) `` ,  with translation and rotation parts  _(x::Vec3, p::RV3)_.   
>


|           |                                                                     |
|:----------|:--------------------------------------------------------------------|
|  _Vec3_ |    vector with 3 components (e.g. position, velocity, Euler angles). |
|  _Mat3_ |    3 x 3 matrix (e.g. rotation and tangent operators). |
|  _RV3_ |  rotation vector `` \\mathbf{n}  \\phi ``. |
| _Quat_ |  quaternion  `` \\mathbf{q} = [q_0,  q_1, q_2, q_3] `` . |
|  _NodeFrame_ |   compound structural type to describe a nodal frame `` \\mathbf{H} (\\mathbf{x}, \\mathbf{p}) `` ,  with translation and rotation parts  _(x::Vec3, p::RV3)_.  | 

"""
module MyAlgebra

using LinearAlgebra
using Printf

"""
    PREC
        constant : PREC = sqrt(eps(Float64))
"""
const PREC=sqrt(eps(Float64))

"""
    PRC
        constant : PRC = eps(Float64)
"""
const PRC = eps(Float64)

Base.show(io::IO, f::Float64) = @printf(io, "%.6e", f)


# 

dir = "my_algebra/"

include(dir*"vec3.jl")
include(dir*"vec3_base.jl")
export Vec3
export dotp,norm2



include(dir*"mat3.jl")
include(dir*"mat3_base.jl")
export Mat3, I3

include(dir*"rv3.jl")
include(dir*"rv3_base.jl")
export RV3
export tilde,crossp

include(dir*"quat.jl")
include(dir*"quat_base.jl")
export Quat

include(dir*"rp3.jl")
include(dir*"rp3_base.jl")
export RP3

include(dir*"rot.jl")
export rot,rot2

include(dir*"invrot.jl")
export invrot

include(dir*"rv3_comp_rule.jl")

include(dir*"tang.jl")
export tang

include(dir*"dtang.jl")
export Dtang

include(dir*"invtang.jl")
export invtang

include(dir*"dinvtang.jl")
export Dinvtang,DinvtangT

include(dir*"euler_to_rv.jl")
export euler_to_RV

include(dir*"NodeFrame.jl")
include(dir*"exp_SE3.jl")
include(dir*"log_SE3.jl")
include(dir*"invT_SE3.jl")
include(dir*"DinvT_SE3.jl")
include(dir*"frame_solve.jl")
include(dir*"frame_solve1.jl")
include(dir*"compute_strains.jl")
include(dir*"Adj.jl")
include(dir*"tilde.jl")

export NodeFrame,exp_SE3,log_SE3,invT_SE3,DinvT_SE3, frame_solve,frame_solve1,compute_strains, Adj,tilde

include(dir*"shape_functions_1D.jl")
include(dir*"shape_functions_2D.jl")
export shape_functions_1D, shape_functions_2D
include(dir*"gauss_points.jl")
export gauss_points
include(dir*"eye.jl")
export eye
include(dir*"diagonal.jl")
export diagonal

"""include(dir*"levicivita.jl")
export levicivita
include(dir*"alpha.jl")
include(dir*"beta.jl")
include(dir*"gamma.jl")
include(dir*"delta.jl")
include(dir*"bernoulli.jl")
export bernoulli,bernoulli2
export alpha,beta,gamma,delta"""



end
