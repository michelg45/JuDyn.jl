"""
    superelement
"""
function superelement(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector,theta_p::Float64,
    alpha_stiff::Float64,matrix::Bool)

    nc = Main.node_container
    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    sec = Main.SuperElements.superelement_container
    sem = Main.SE_matrices


    iel = findfirst(x -> x == nbr, sec.numbers)[1]

    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed)
    grav = mc.gravity

    #
    # Getting dimensions of the superelement from superelement matrix file "SE_matrices"
    # Allocation of iteration matrix, residual and load vectors
    #

    imx = sec.matrix_sets[iel]

    Nrig = sem[imx].Nrig
    N_B = sem[imx].N_B
    N_I = sem[imx].N_I
    K = sem[imx].K
    M = sem[imx].M

    iq_B = [i for i in 1:N_B]
    K_BB = K[iq_B,iq_B]
    iq = [i for i in 1:N_B+N_I]

    N_I > 0 && (init_qI = N_B;
        iq_I = [i for i in 1:N_I]  .+ init_qI;
        K_II = K[iq_I,iq_I];
        K_BI = K[iq_B,iq_I])

    rotation == true && (
        S_rot = sem[imx].S_rot;
        A_rot = sem[imx].A_rot)

    ndim = (N_B + N_I + Nrig)*2


    S_el = Array{Float64,2}(zeros(ndim,ndim))
    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    DqDx = Array{Float64,2}(zeros(N_B,N_B+6))
    DqDw = Array{Float64,2}(zeros(N_B,N_B+6))
    DqdotDx = Array{Float64,2}(zeros(N_B,N_B+6))

    #
    # Getting detailed characteristics  of the superelement from "superelement_container"
    #

    node_0 = sec.ref_node_orders[iel]
    bnodes = sec.boundary_node_orders[iel]
    Nbnodes = sec.boundary_node_numbers[iel]
    Nmodes = sec.internal_mode_numbers[iel]
    int_dofs = sec.internal_mode_numbers[iel]
    node_components = sec.boundary_node_components[iel]
    X_B = sec.local_node_coordinates[iel]
    psi_rel_B = sec.local_node_orientations[iel]
    name = sec.superelement_names[iel]
    Mrig = sem[imx].Mrig
    mass = Mat3(Mrig[1:3,1:3])
    J = Mat3(Mrig[4:6,4:6])
   
    nl_corr = sec.nl_correction[iel]

    #
    # Definition of arrays for selective treatment of nodes
    #

    ix_B_tr = Array{Int,2}(zeros(Nbnodes,3))
    ix_B_rot = Array{Int,2}(zeros(Nbnodes,3))
    irotB = Vector{Bool}(undef,Nbnodes)
    R0txB = Array{Float64,2}(zeros(3,Nbnodes))

    #
    # getting localization of the superelement from "element_container"
    #

    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]
    inv_loc_int = ec.inv_loc_int[iel]
    inv_loc_q = [inv_loc_x;inv_loc_int]
    inv_loc_v = ec.inv_loc_v[iel]
    inv_loc =[inv_loc_q;inv_loc_v.+ mc.Ndofs_q]


    #
    # pull current solution at element level
    #

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)


    #
    # localization of degrees of freedom
    #

    ix_0 = [i for i in 1:3]
    itheta_0 = ix_0 .+ 3
    init_v  = N_B + N_I + 6
    iv_0 = ix_0 .+ init_v
    iomega_0 = itheta_0 .+ init_v

    
    ix_B = iq_B .+6
    ix = iq .+6
    ix_0B = [ix_0; itheta_0; ix_B]
    iv = ix .+ init_v
    iv_B = ix_B .+ init_v
    N_I > 0 &&  (ix_I = iq_I .+ 6; iv_I = ix_I .+ init_v)



    #
    # Rigib body contribution of the reference node
    #

    v_0 = Vec3(s[iv_0]+Ds[iv_0])
    vdot_0 = Vec3(sdot[iv_0])
    Omega_0= Vec3(s[iomega_0]+Ds[iomega_0])
    Omegadot_0= Vec3(sdot[iomega_0])

    Dx_0 = Vec3(Ds[ix_0])
    Dpsi_0 = RV3(Ds[itheta_0])
    x_0 = cf[nodes[3]].x + Dx_0
    psi_0 = RV3(cf[nodes[3]].psi,Dpsi_0)
    xdot_0 = Vec3(sdot[ix_0])
    R_0 = rot(psi_0)
    R_0T =  transpose(R_0)
    T_0 = tang(Dpsi_0)
    psidot_0 = Vec3(sdot[itheta_0])
    W_0 = T_0*psidot_0
    DW_0 = Dtang(psidot_0,Dpsi_0)
 
  

    W0R0T = tilde(W_0)*R_0T

    p_el[ix_0] = (mass*grav).v

if rotation == true 

    ROmega_d = rot(-psi_0,Omega_d)
    tOmega_d =  tilde(Omega_d).mat
    res_el[ix_0] =  (mass*(grav - vdot_0 - crossp(Omega_d,v_0))).v
    res_el[iv_0]= (v_0 - xdot_0 - crossp(Omega_d,x_0)).v 
    res_el[itheta_0] = -(J*Omegadot_0 + crossp(W_0+ROmega_d, J*Omega_0)).v
    res_el[iomega_0]= (Omega_0 - W_0 - ROmega_d).v
    S_el[ix_0,iv_0]= mass.mat*(theta_p*eye(3) + alpha_stiff*tOmega_d)
    S_el[iv_0,ix_0]= theta_p*eye(3) + alpha_stiff*tOmega_d
    A_0= (theta_p*I3 + alpha_stiff*tilde(ROmega_d))*T_0 + alpha_stiff*DW_0
    S_el[itheta_0,iomega_0]= (theta_p*J+alpha_stiff*tilde(W_0+ROmega_d)*J).mat
    
else 
    res_el[ix_0] =  (mass*(grav -vdot_0)).v
    res_el[itheta_0] = -(J*Omegadot_0 + crossp(W_0,J*Omega_0)).v
    res_el[iv_0]= (v_0 - xdot_0).v
    res_el[iomega_0]= (Omega_0 -W_0).v
    A_0= theta_p*T_0 + alpha_stiff*DW_0
    S_el[itheta_0,iomega_0]= (theta_p*J+alpha_stiff*tilde(W_0)*J).mat
    S_el[ix_0,iv_0]=theta_p*mass.mat
    S_el[iv_0,ix_0]= theta_p*eye(3)

end
    CM_kin_el = 1.0/2.0*(dotp(mass*v_0,v_0)+dotp(Omega_0,J*Omega_0))


    S_el[iv_0,iv_0]=  - alpha_stiff*eye(3)
    S_el[itheta_0,itheta_0]= -(tilde(J*Omega_0)*A_0).mat
    S_el[iomega_0,iomega_0]= - alpha_stiff*eye(3)
    S_el[iomega_0,itheta_0]  =  A_0.mat

    id_0 = [ix_0; itheta_0]
    idv_0 = [iv_0; iomega_0]

    index0 = [id_0; idv_0]


    #
    # Preparing contribution of internal modes
    #

    qdot_I = sdot[ix_I]
    q_I = s[ix_I] + Ds[ix_I]

    v = s[iv]+Ds[iv]
    vdot = sdot[iv]

    #
    # Preparing contribution of boundary  modes
    #

    init_xP = 6



    q_B    = Vector{Float64}[]
    qdot_B    = Vector{Float64}[]

    for iP = 1:Nbnodes

        X_P = X_B[iP]




        ncomp = node_components[iP]



        if ncomp == 6  
            psi_rel_P = psi_rel_B[iP]
            R_rel_P = rot(psi_rel_P)
            irotB[iP] = true 
        else
             irotB[iP] = false
        end

        node_P = bnodes[iP]

        ix_P  = [i for i in 1:ncomp] .+ init_xP
        iq_P  = ix_P .- 6
        ix_0P  = [ix_0; itheta_0; ix_P]
        ix_P_tr  = ix_P[1:3]
        ix_B_tr[iP,:] = ix_P_tr

        xdot_P = Vec3(sdot[ix_P_tr])

        Bp = R_0T

        x_P = cf[node_P].x + Vec3(Ds[ix_P_tr])

        Xpu = rot(-psi_0, (x_P - x_0))
        tXpu = tilde(Xpu)
        tXpuT = tXpu*T_0
        R0txB[:,iP] = Xpu.v

        u_P = Xpu - X_P
        udot_P= rot(-psi_0,(xdot_P - xdot_0)) + crossp(Xpu,W_0)
        coef = tilde(udot_P)*T_0 + tilde(W_0)*tXpuT + tXpu*DW_0

        if irotB[iP] == true

            itheta_P = ix_P[4:6]
            Dpsi_P = RV3(Ds[itheta_P])
            psidot_P = Vec3(sdot[itheta_P])
            psi_P = RV3(cf[node_P].psi,Dpsi_P) 
            T_P = tang(Dpsi_P)
            W_P = T_P*psidot_P
            DW_P = Dtang(psidot_P, psi_P)

            R_P = rot(psi_P)

            phi_P = RV3(-psi_0, RV3(psi_P,-psi_rel_P))
            phidot_P = R_rel_P *W_P - W_0
            Tm_phi = invtang(phi_P)
            Tm_phiT = transpose(Tm_phi)
            ix_B_rot[iP,:] = copy(itheta_P)

            if nl_corr == true

                q_P = [(Tm_phiT*u_P).v; phi_P.v]
                qdot_P = [udot_P.v; phidot_P.v]

                DinvT = DinvtangT(u_P,phi_P)

                Ap = Tm_phiT*tXpu- DinvT
                Bp  = (Tm_phiT*Bp).mat

                DqDx[iq_P,ix_0P] =    [-Bp  (Ap*T_0).mat  Bp  (DinvT*R_rel_P*T_P).mat; zeros(3,3) -(Tm_phiT*T_0).mat zeros(3,3) (Tm_phi*R_rel_P *T_P).mat]
                DqDw[iq_P,ix_0P] =    [-Bp  Ap.mat Bp (DinvT*R_rel_P).mat; zeros(3,3) -Tm_phiT.mat zeros(3,3) (Tm_phi*R_rel_P).mat]
                DqdotDx[iq_P,ix_0P] = [W0R0T.mat  coef.mat  -W0R0T.mat  zeros(3,3); zeros(3,3) -DW_0.mat  zeros(3,3)  (R_rel_P*DW_P).mat]

            else
                q_P = [u_P.v; phi_P.v]
                qdot_P = [udot_P.v; phidot_P.v]
                Ap = tXpu
                Bp  = Bp.mat
                DqDx[iq_P,ix_0P] =    [-Bp  (Ap*T_0).mat  Bp  zeros(3,3); zeros(3,3) -(Tm_phiT*T_0).mat zeros(3,3) (Tm_phi*R_rel_P *T_P).mat]
                DqDw[iq_P,ix_0P] =    [-Bp  Ap.mat Bp zeros(3,3); zeros(3,3) -Tm_phiT.mat zeros(3,3) (Tm_phi*R_rel_P).mat]
                DqdotDx[iq_P,ix_0P] = [W0R0T.mat  coef.mat  -W0R0T.mat  zeros(3,3); zeros(3,3) -DW_0.mat  zeros(3,3)  (R_rel_P*DW_P).mat]
            end

        else
            q_P = copy(u_P.v)
            qdot_P =  copy(udot_P.v)


            DqDx[iq_P,ix_0P] =    [-Bp.mat  tXpuT.mat  Bp.mat  ]
            DqdotDx[iq_P,ix_0P] = [W0R0T.mat  coef.mat  -W0R0T.mat ]
            DqDw[iq_P,ix_0P] = DqDx[iq_P,ix_0P]
        end

        q_B = [q_B; q_P]
        qdot_B = [qdot_B; qdot_P]

        init_xP += ncomp

    end

    qdot = [qdot_B; qdot_I]
    q = [q_B; q_I]


# Contribution of boundary nodes in local coordinates


    


    res_el[iv] = v - qdot
    Fel = K*q
    mv = M*v
    elast_kin_el = 1.0/2.0*transpose(v)*mv
    res_el[ix] = - (Fel + M*vdot)
    S_el[iv,iv] = -alpha_stiff*eye(N_B+N_I)
    S_el[ix,iv] = theta_p*M



    Om = W_0 + Omega_d
    A_omega = Om.v'*A_rot
    G = - (Om[1]*S_rot[1]+Om[2]*S_rot[2]+Om[3]*S_rot[3])
    res_el[iv] -= A_omega*q
    res_el[ix] -= G*v
    S_el[ix,iv] += alpha_stiff*G
    S_el[iv,iq] += alpha_stiff*A_omega
    
    
    S_el[iv_B,ix_0B] = theta_p*DqDx + alpha_stiff*DqdotDx
    S_el[iv_I,ix_I] = theta_p*eye(N_I)
    S_el[ix_I,ix_I] = alpha_stiff*K_II
    S_el[ix_B,ix_I] = alpha_stiff*K_BI
    S_el[ix_I,ix_0B] =  alpha_stiff*transpose(K_BI)*DqDx
    S_el[ix_B,ix_0B] = alpha_stiff*K_BB*DqDx
    str_el = 1.0/2.0*transpose(q)*Fel



#
#   Contribution of elatic modes to RB motion
#


if rotation == true
    s_qv = Vec3([q'*S_rot[i]*v  for i in 1:3]) 
    s_qvdot = [q'*S_rot[i]*vdot  for i in 1:3] 
    s_vqdot = [v'*S_rot[i]*qdot  for i in 1:3]
    res_el[itheta_0] -= (crossp(Om,s_qv)).v - s_vqdot + s_qvdot
end

#
# transformation of boundary node contribution to global coordinates
#

    res_el[id_0] += transpose(DqDw[:,id_0])*res_el[ix_B]
    S_el[id_0,:] += transpose(DqDw[:,id_0])*S_el[ix_B,:]
    res_el[ix_B] = transpose(DqDw[:,ix_B])*res_el[ix_B]
    S_el[ix_B,:] = transpose(DqDw[:,ix_B])*S_el[ix_B,:]

    res_el[ix_0] += p_el[ix_0]

    pot_el = - transpose(p_el[ix_0])*x_0.v

    kin_el = CM_kin_el + elast_kin_el

    if matrix == false
        push_element_sparse(res,p,iel,inv_loc_q,inv_loc_v,S_el,res_el,p_el)
        return pot_el, kin_el, str_el
    else
        return inv_loc, S_el = condensed_element_matrix(iel,inv_loc_q,inv_loc_v,S_el)
    end


    

end

"""
    superelement_force
"""
function superelement_force(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector)

    mc =  Main.model_container
    cf =  Main.Frames.current_frames
    ec =  Main.element_container
    sec = Main.SuperElements.superelement_container
    sem = Main.SE_matrices


    iel = findfirst(x -> x == nbr, sec.numbers)[1]

    rotation = mc.uniform_rotation
    rotation == true && (Omega_d = mc.rotation_speed)
    grav = mc.gravity

    #
    # Getting dimensions of the superelement from superelement matrix file "SE_matrices"
    # Allocation of iteration matrix, residual and load vectors
    #

    imx = sec.matrix_sets[iel]

    Nrig = sem[imx].Nrig
    N_B = sem[imx].N_B
    N_I = sem[imx].N_I
    K = sem[imx].K
    M = sem[imx].M

    iq_B = [i for i in 1:N_B]
    iq = [i for i in 1:N_B+N_I]

    N_I > 0 && (init_qI = N_B;  iq_I = [i for i in 1:N_I]  .+ init_qI)

    rotation == true && (
        S_rot = sem[imx].S_rot;
        A_rot = sem[imx].A_rot)

    ndim = (N_B + N_I + Nrig)*2

    res_el  = Vector{Float64}(zeros(ndim))
    p_el  = Vector{Float64}(zeros(ndim))

    DqDw = Array{Float64,2}(zeros(N_B,N_B+6))


    #
    # Getting detailed characteristics  of the superelement from "superelement_container"
    #

    node_0 = sec.ref_node_orders[iel]
    bnodes = sec.boundary_node_orders[iel]
    Nbnodes = sec.boundary_node_numbers[iel]
    node_components = sec.boundary_node_components[iel]
    X_B = sec.local_node_coordinates[iel]
    psi_rel_B = sec.local_node_orientations[iel]
    Mrig = sem[imx].Mrig
    mass = Mat3(Mrig[1:3,1:3])
    J = Mat3(Mrig[4:6,4:6])
   
    nl_corr = sec.nl_correction[iel]

    #
    # Definition of arrays for selective treatment of nodes
    #

    ix_B_tr = Array{Int,2}(zeros(Nbnodes,3))
    ix_B_rot = Array{Int,2}(zeros(Nbnodes,3))
    irotB = Vector{Bool}(undef,Nbnodes)
    R0txB = Array{Float64,2}(zeros(3,Nbnodes))

    #
    # getting localization of the superelement from "element_container"
    #

    iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
    inv_loc_x = ec.inv_loc_x[iel]
    inv_loc_int = ec.inv_loc_int[iel]
    inv_loc_q = [inv_loc_x;inv_loc_int]
    inv_loc_v = ec.inv_loc_v[iel]
    inv_loc =[inv_loc_q;inv_loc_v.+ mc.Ndofs_q]


    #
    # pull current solution at element level
    #

    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)


    #
    # localization of degrees of freedom
    #

    ix_0 = [i for i in 1:3]
    itheta_0 = ix_0 .+ 3
    init_v  = N_B + N_I + 6
    iv_0 = ix_0 .+ init_v
    iomega_0 = itheta_0 .+ init_v

    
    ix_B = iq_B .+6
    ix = iq .+6
    iv = ix .+ init_v
    N_I > 0 &&  (ix_I = iq_I .+ 6)



    #
    # Rigib body contribution of the reference node
    #

    v_0 = Vec3(s[iv_0]+Ds[iv_0])
    vdot_0 = Vec3(sdot[iv_0])
    Omega_0= Vec3(s[iomega_0]+Ds[iomega_0])
    Omegadot_0= Vec3(sdot[iomega_0])
    Dx_0 = Vec3(Ds[ix_0])
    Dpsi_0 = RV3(Ds[itheta_0])
    x_0 = cf[nodes[3]].x + Dx_0
    psi_0 = RV3(cf[nodes[3]].psi,Dpsi_0)
    xdot_0 = Vec3(sdot[ix_0])
    R_0 = rot(psi_0)
    R_0T =  transpose(R_0)
    T_0 = tang(Dpsi_0)
    psidot_0 = Vec3(sdot[itheta_0])
    W_0 = T_0*psidot_0


    p_el[ix_0] = (mass*grav).v

if rotation == true 

    ROmega_d = rot(-psi_0,Omega_d)
    res_el[ix_0] =  (mass*(grav - vdot_0 - crossp(Omega_d,v_0))).v
    res_el[iv_0]= (v_0 - xdot_0 - crossp(Omega_d,x_0)).v 
    res_el[itheta_0] = -(J*Omegadot_0 + crossp(W_0+ROmega_d, J*Omega_0)).v
    res_el[iomega_0]= (Omega_0 - W_0 - ROmega_d).v
    
else 
    res_el[ix_0] =  (mass*(grav -vdot_0)).v
    res_el[itheta_0] = -(J*Omegadot_0 + crossp(W_0,J*Omega_0)).v
    res_el[iv_0]= (v_0 - xdot_0).v
    res_el[iomega_0]= (Omega_0 -W_0).v

end
    CM_kin_el = 1.0/2.0*(dotp(mass*v_0,v_0)+dotp(Omega_0,J*Omega_0))


    id_0 = [ix_0; itheta_0]





    #
    # Preparing contribution of internal modes
    #

    qdot_I = sdot[ix_I]
    q_I = s[ix_I] + Ds[ix_I]

    v = s[iv]+Ds[iv]
    vdot = sdot[iv]

    #
    # Preparing contribution of boundary  modes
    #

    init_xP = 6



    q_B    = Vector{Float64}[]
    qdot_B    = Vector{Float64}[]

    for iP = 1:Nbnodes

        X_P = X_B[iP]




        ncomp = node_components[iP]



        if ncomp == 6  
            psi_rel_P = psi_rel_B[iP]
            R_rel_P = rot(psi_rel_P)
            irotB[iP] = true 
        else
             irotB[iP] = false
        end

        node_P = bnodes[iP]

        ix_P  = [i for i in 1:ncomp] .+ init_xP
        iq_P  = ix_P .- 6
        ix_0P  = [ix_0; itheta_0; ix_P]
        ix_P_tr  = ix_P[1:3]
        ix_B_tr[iP,:] = ix_P_tr

        xdot_P = Vec3(sdot[ix_P_tr])

        Bp = R_0T

        x_P = cf[node_P].x + Vec3(Ds[ix_P_tr])

        Xpu = rot(-psi_0, (x_P - x_0))
        tXpu = tilde(Xpu)
        R0txB[:,iP] = Xpu.v

        u_P = Xpu - X_P
        udot_P= rot(-psi_0,(xdot_P - xdot_0)) + crossp(Xpu,W_0)

        if irotB[iP] == true

            itheta_P = ix_P[4:6]
            Dpsi_P = RV3(Ds[itheta_P])
            psidot_P = Vec3(sdot[itheta_P])
            psi_P = RV3(cf[node_P].psi,Dpsi_P) 
            T_P = tang(Dpsi_P)
            W_P = T_P*psidot_P


            R_P = rot(psi_P)

            phi_P = RV3(-psi_0, RV3(psi_P,-psi_rel_P))
            phidot_P = R_rel_P *W_P - W_0
            Tm_phi = invtang(phi_P)
            Tm_phiT = transpose(Tm_phi)
            ix_B_rot[iP,:] = copy(itheta_P)

            if nl_corr == true

                q_P = [(Tm_phiT*u_P).v; phi_P.v]
                qdot_P = [udot_P.v; phidot_P.v]

                DinvT = DinvtangT(u_P,phi_P)

                Ap = Tm_phiT*tXpu- DinvT
                Bp  = (Tm_phiT*Bp).mat

                 DqDw[iq_P,ix_0P] = [-Bp  Ap.mat Bp (DinvT*R_rel_P).mat; zeros(3,3) -Tm_phiT.mat zeros(3,3) (Tm_phi*R_rel_P).mat]


            else
                q_P = [u_P.v; phi_P.v]
                qdot_P = [udot_P.v; phidot_P.v]
                Ap = tXpu
                Bp  = Bp.mat

                DqDw[iq_P,ix_0P] =    [-Bp  Ap.mat Bp zeros(3,3); zeros(3,3) -Tm_phiT.mat zeros(3,3) (Tm_phi*R_rel_P).mat]

            end
###
        else
            tXpuT = tXpu*T_0
            q_P = copy(u_P.v)
            qdot_P =  copy(udot_P.v)
            DqDw[iq_P,ix_0P] = [-Bp.mat  tXpuT.mat  Bp.mat  ]
        end

        q_B = [q_B; q_P]
        qdot_B = [qdot_B; qdot_P]

        init_xP += ncomp

    end

    qdot = [qdot_B; qdot_I]
    q = [q_B; q_I]

# Contribution of boundary nodes in local coordinates

    res_el[iv] = v - qdot
    Fel = K*q
    mv = M*v
    elast_kin_el = 1.0/2.0*transpose(v)*mv
    res_el[ix] = - (Fel + M*vdot)

    Om = W_0 + Omega_d
    A_omega = Om.v'*A_rot
    G = - (Om[1]*S_rot[1]+Om[2]*S_rot[2]+Om[3]*S_rot[3])
    res_el[iv] -= A_omega*q
    res_el[ix] -= G*v

    str_el = 1.0/2.0*transpose(q)*Fel

#
#   Contribution of elatic modes to RB motion
#

if rotation == true
    s_qv = Vec3([q'*S_rot[i]*v  for i in 1:3]) 
    s_qvdot = [q'*S_rot[i]*vdot  for i in 1:3] 
    s_vqdot = [v'*S_rot[i]*qdot  for i in 1:3]
    res_el[itheta_0] -= (crossp(Om,s_qv)).v - s_vqdot + s_qvdot
end

#
# transformation of boundary node contribution to global coordinates
#

    res_el[id_0] += transpose(DqDw[:,id_0])*res_el[ix_B]
    res_el[ix_B] = transpose(DqDw[:,ix_B])*res_el[ix_B]


    res_el[ix_0] += p_el[ix_0]

    pot_el = - transpose(p_el[ix_0])*x_0.v

    kin_el = CM_kin_el + elast_kin_el

    for i = 1:ndim
        iloc = inv_loc[i] 
        iloc > 0  &&  (res[iloc] += res_el[i]; p[iloc] += p_el[i])
    end

    return pot_el, kin_el, str_el

end

