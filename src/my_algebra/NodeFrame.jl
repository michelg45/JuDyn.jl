"""
    NodeFrame

        NodeFrame is a structural type to describe a nodal frame. 
        Its translation part is described by a vector p::Vec3 and
        is rotation part by a rotation vector p::RV3.

        mutable struct NodeFrame
            x::Vec3
            p::RV3
        end

```@example

        a1 = Vec3(0.0,1.0,0.0)
        Vec3([0.000000e+00, 1.000000e+00, 0.000000e+00])
        phi1 = RV3(0.0, 2.0, -1.0)
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00])
        H1 = NodeFrame(a1,phi1)
        NodeFrame(Vec3([0.000000e+00, 1.000000e+00, 0.000000e+00]), 
        RV3([0.000000e+00, 2.000000e+00, -1.000000e+00]))

```

"""
mutable struct NodeFrame
    x::Vec3
    p::RV3
end

function Base.:*(h1::NodeFrame, h2::NodeFrame)
    x = h1.x + rot(h1.p, h2.x)
    p = RV3(h1.p, h2.p)
    return NodeFrame(x,p)
end

function Base.:-(h::NodeFrame)
    x = -rot(-h.p, h.x)
    p = -h.p
    return NodeFrame(x,p)
end

NodeFrame() = NodeFrame(Vec3(), RV3())

