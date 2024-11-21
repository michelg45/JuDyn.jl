function test_top_with_frame_link()
    # ===================================================================
    # introduction og general model data 
    # ===================================================================
    gravity = Vec3(0.,0.,-9.81)
    name = "top_with_frame_link"
    update = true
    set_general_data(name,gravity,update)
    # =====================================================================
    # frame orientation (rotation vector)  computed from Euler angles.
    # =====================================================================
    psi = 0.0
    theta = pi/9.0
    phi = 0.0
    Rot2 = euler_to_RV(phi,theta,psi)
    # ========================================================================
    # initial material angular velocity computed  from Euler angle derivatives.
    # ========================================================================
    phidot = 50.0
    psidot = -10.0
    thetadot = 0.0
    axe1 = Vec3(1.0,0.0,0.0)
    axe3 = Vec3(0.0,0.0,1.0)
    # ======================================================================
    # node input
    # ======================================================================
    X_2 = Vec3(0.,0.,1.3)
    x_1 = Vec3()
    x_2 = rot(Rot2,X_2)
    Rot1 = RV3()
    set_node(10,x_1,Rot1,"frame")
    set_node(20,x_2,Rot2,"frame")
    # =====================================================================
    # closing node input
    # =====================================================================
    end_nodes()
    # =====================================================================
    # element input
    # =====================================================================
    mass = 5.0
    inertia = Vec3(0.8,0.8,1.8)
    set_rigid_body(1,20,mass,inertia)
    set_frame_link(2,10,20)
    # =====================================================================
    # end of element input
    # =====================================================================
    end_elements()
    # =====================================================================
    # set boundary conditions
    # =====================================================================
    set_node_BC(10,"pinned")
    end_BC()
    # =====================================================================
    # model assembly
    # =====================================================================
    sol_type = "dynamic"
    assemble(sol_type)
    # =====================================================================
    # initial conditions input
    # =====================================================================
    Omega=axe3*phidot+rot(-phi,3)*(axe1*thetadot+rot(-theta,1)*axe3*psidot)
    v_2=rot(Rot2,crossp(Omega,X_2))
    println(v_2,Omega)
    initial_conditions()
    set_initial_velocity(20,v_2,Omega)
    end_initial_conditions()
    # =====================================================================
    # response computation : dynamic response is using generalized-α is the 
    # standard. Solution parameters are read from JSON_file.
    # ======================================================================
    input_dir = "test/top_tests/"
    JSON_file = input_dir*"top.json"
    h5_file = solve(JSON_file,sol_type)

    global times, y, ydot, p, str_energy, pot_energy, kin_energy, ext_work, nitmax  = read_results_dynamic(h5_file)

    min_pot_energy = findmin(pot_energy)[1]
    min_kin_energy = findmin(kin_energy)[1]
    max_pot_energy = findmax(pot_energy)[1]
    max_kin_energy = findmax(kin_energy)[1]

    rm(h5_file)

    return  min_pot_energy, max_pot_energy, min_kin_energy, max_kin_energy
end
