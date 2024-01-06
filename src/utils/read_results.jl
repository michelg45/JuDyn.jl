
    using HDF5


    function read_results_dynamic(h5_file)

        """
            "read_results" is a function reading the full content of the HDF5 file
            collecting the results of the analysis.

        """

    file=h5open(h5_file,"r")
    gtimes = file["times"]
    gp=file["p"]
    gy=file["y"]
    gydot=file["ydot"]
    gstr = file["str_energy"]
    gpot = file["pot_energy"]
    gkin = file["kin_energy"]
    gwork = file["ext_work"]
    gnitmax = file["nitmax"]
    gmdata = file["model_data"]
    g_node_pos = file["node_positions"]
    times = read(gtimes, "times")
    p = read(gp, "p")
    y = read(gy, "y")
    ydot = read(gydot, "ydot")
    str_energy = read(gstr, "str_energy")
    pot_energy = read(gpot, "pot_energy")
    kin_energy = read(gkin, "kin_energy")
    ext_work = read(gwork, "ext_work")
    nitmax = read(gnitmax, "nitmax")
    model_data = read(gmdata, "model_data")
    node_positions= read(g_node_pos,"node_positions")

    close(file)

    return times, y, ydot, p, str_energy, pot_energy, kin_energy, ext_work, nitmax, node_positions, model_data
    end

    function read_results_static(h5_file)

        """
            "read_results" is a function reading the full content of the HDF5 file
            collecting the results of the analysis.

        """

    file=h5open(h5_file,"r")
    gtimes = file["times"]
    gp=file["p"]
    gy=file["y"]
    gydot=file["ydot"]
    gstr = file["str_energy"]
    gpot = file["pot_energy"]
    gwork = file["ext_work"]
    gnitmax = file["nitmax"]
    gmdata = file["model_data"]
    g_node_pos = file["node_positions"]
    times = read(gtimes, "times")
    p = read(gp, "p")
    y = read(gy, "y")
    ydot = read(gydot, "ydot")
    str_energy = read(gstr, "str_energy")
    pot_energy = read(gpot, "pot_energy")
    ext_work = read(gwork, "ext_work")
    nitmax = read(gnitmax, "nitmax")
    node_positions= read(g_node_pos,"node_positions")
    model_data = read(gmdata, "model_data")
    close(file)

    return times, y, ydot, p, str_energy, pot_energy,  ext_work, nitmax, node_positions, model_data


    end

    function read_results(h5_file,sol_type)

        sol_type == "static" ?  read_results_static(h5_file) : read_results_dynamic(h5_file)

    end

    function read_results(h5_file)

        read_results_dynamic(h5_file)

    end

    
    function read_results_constraints(h5_file)

    file=h5open(h5_file,"r")
    g_mult = file["multipliers"]
    g_bounds = file["bounds"]

    lambda = read(g_mult, "multipliers")
    bounds = read(g_bounds, "bounds")

    close(file)

    return lambda, bounds


    end