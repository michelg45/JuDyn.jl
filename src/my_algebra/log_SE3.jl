function log_SE3(h::NodeFrame)
    return log_SE3 = [invtang(-h.p, h.x).v; (h.p).v]
end