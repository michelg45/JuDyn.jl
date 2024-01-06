"""
    NodeFrame

    NodeFframe is a structural type to describe a nodal frame. Its translation part is descrived by a ::Vec3 type
    and is rotation part by a ::RV3 type.

        mutable struct NodeFrame
            x::Vec3
            p::RV3
        end

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

