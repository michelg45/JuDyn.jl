function log_SO3(h::NodeFrame)
    return log_SO3 = [(h.x).v; (h.p).v]
end