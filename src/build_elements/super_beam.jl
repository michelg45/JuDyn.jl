"""
    super_beam

        function constructing the iteration matrix residual vector  of a super-beam element. The element is constructed as a superelement resulting from the assembly of 2 beam elements with cubic interpolation of the deflection. Its stiffness and mass kernel has been construted by Herting's method by the 'super_beam_matrix_kernel.jl' function.

            Calling sequence: super_beam_force(nbr,y,Dy,ydot, res,p,alpha_stiff,gamma_p,matrix,itime,h)
"""
function super_beam(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector,
    alpha_stiff::Float64,gamma_p::Float64,matrix::Bool,itime::Int,h::Float64)

 nc = Main.node_container
 mc =  Main.model_container
 cf =  Main.Frames.current_frames
 ec =  Main.element_container
 sbc = SetElements.super_beam_container

 iel = findfirst(x -> x == nbr, sbc.numbers)[1]


 ndim = 36

 S_el = Array{Float64,2}(zeros(ndim,ndim))
 res_el  = Vector{Float64}(zeros(ndim))
 p_el  = Vector{Float64}(zeros(ndim))

 DqDx = Array{Float64,2}(zeros(12,18))
 DqDw = Array{Float64,2}(zeros(12,18))
 DqdotDx = Array{Float64,2}(zeros(12,18))

 #
 # Getting detailed characteristics  of the superelement from "superelement_container"
 #

 node_0 = sbc.ref_node_orders[iel]
 bnodes = sbc.boundary_node_orders[iel]
 nodes = [bnodes; node_0]

 length = sbc.length[iel]
 psi_rel = sbc.local_node_orientations[iel]
 R_rel = [rot(psi_rel[i]) for i in 1:2]

 mass = sbc.mass[iel]
 Jrot = Mat3(diagm(sbc.Jrot[iel]))

 
 rotation = mc.uniform_rotation
 rotation == true && (Omega_d = mc.rotation_speed; tOmega_d = tilde(Omega_d).mat; 
 S_elast = sbc.S[iel])

 if rotation == true
    time = itime*h
    phi_d = RV3(time*Omega_d)
    grav = rot(-phi_d, mc.gravity)
else
    grav = mc.gravity
end

 K_elast = sbc.K_elast[iel]
 M_elast = sbc.M_elast[iel]
 

 nl_correction = sbc.nl_correction[iel]

 #
 # getting localization of the superelement from "element_container"
 #

 iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
 inv_loc_x = ec.inv_loc_x[iel]
 inv_loc_v = ec.inv_loc_v[iel]

 inv_loc =[inv_loc_x; inv_loc_v.+ mc.Ndofs_q]

 #
 # pull current solution at element level
 #
 (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

 ix_P = Vector{Vector{Int}}(undef,2)
 x_P = Vector{Vec3}(undef,2)
 X_P = Vector{Vec3}(undef,2)
 Dx_P = Vector{Vec3}(undef,2)
 xdot_P = Vector{Vec3}(undef,2)
 itheta_P = Vector{Vector{Int}}(undef,2)
 psi_P = Vector{RV3}(undef,2)
 Dpsi_P = Vector{RV3}(undef,2)
 psidot_P = Vector{Vec3}(undef,2)

 R_P = Vector{Mat3}(undef,2)
 T_P = Vector{Mat3}(undef,2)
 W_P = Vector{Vec3}(undef,2)
 DW_P = Vector{Mat3}(undef,2)
 u_P = Vector{Vec3}(undef,2)
 u = Vector{Vec3}(undef,2)
 udot_P = Vector{Vec3}(undef,2)
 phi_P = Vector{RV3}(undef,2)
 Tm1P = Vector{Mat3}(undef,2)
 tXpu = Vector{Mat3}(undef,2)

 for ip = 1:2

     ix_P[ip] = [(ip-1)*6+i for i in 1:3]
     itheta_P[ip] = ix_P[ip] .+ 3
     xdot_P[ip] = Vec3(sdot[ix_P[ip]])
     psidot_P[ip] = Vec3(sdot[itheta_P[ip]])
     Dx_P[ip] = Vec3(Ds[ix_P[ip]])
     Dpsi_P[ip] = RV3(Ds[itheta_P[ip]])
     psi_P[ip]= RV3(cf[nodes[ip]].psi,Dpsi_P[ip])
     x_P[ip]= cf[nodes[ip]].x + Dx_P[ip]
     R_P[ip] = rot(psi_P[ip])
     T_P[ip] = tang(Dpsi_P[ip])
     W_P[ip] = T_P[ip]*psidot_P[ip]
     DW_P[ip] = Dtang(psidot_P[ip],Dpsi_P[ip])

 end

ix_0 = ix_P[1] .+ 12
iv_0 = ix_0 .+ 18
itheta_0 = itheta_P[1] .+ 12
iomega_0 = itheta_0 .+18

Dx_0 = Vec3(Ds[ix_0])
Dpsi_0 = RV3(Ds[itheta_0])
x_0 = cf[nodes[3]].x + Dx_0
psi_0 = RV3(cf[nodes[3]].psi,Dpsi_0)
xdot_0 = Vec3(sdot[ix_0])
psidot_0 = Vec3(sdot[itheta_0])
R_0 = rot(psi_0)
T_0 = tang(Dpsi_0)
W_0 = T_0*psidot_0
DW_0 = Dtang(psidot_0,Dpsi_0)

W0R0T = tilde(W_0)*transpose(R_0)



v_0 = Vec3(s[iv_0]+Ds[iv_0])
vdot_0 = Vec3(sdot[iv_0])
Omega_0 = Vec3(s[iomega_0]+Ds[iomega_0])
Omegadot_0 = Vec3(sdot[iomega_0])

iq_P = [ix_P[1]; itheta_P[1]; ix_P[2]; itheta_P[2]]
iv_P = iq_P .+ 18
ix_0P = [iq_P; ix_0; itheta_0]
v_elast = s[iv_P] + Ds[iv_P]
vdot_elast = sdot[iv_P]

p_el[ix_0] = (mass*grav).v

if rotation == true 
    ROmega_d = rot(-psi_0,Omega_d)
    res_el[ix_0] =  (mass*(grav - vdot_0 - crossp(Omega_d,v_0))).v
    res_el[iv_0]= (v_0 - xdot_0 - crossp(Omega_d,x_0)).v 
    res_el[itheta_0] = -(Jrot*Omegadot_0 + crossp(W_0+ROmega_d, Jrot*Omega_0)).v
    res_el[iomega_0]= (Omega_0 - W_0 - ROmega_d).v
    S_el[ix_0,iv_0]= mass*(gamma_p*eye(3) + alpha_stiff*tOmega_d)
    S_el[iv_0,ix_0]= gamma_p*eye(3) + alpha_stiff*tOmega_d
    A_0= (gamma_p*I3 + alpha_stiff*tilde(ROmega_d))*T_0 + alpha_stiff*DW_0
    S_el[itheta_0,iomega_0]= (gamma_p*Jrot + alpha_stiff*tilde(W_0+ROmega_d)*Jrot).mat
    
else 

    res_el[ix_0] =  mass*(grav - vdot_0).v
    res_el[itheta_0] = -(Jrot*Omegadot_0 + crossp(W_0,Jrot*Omega_0)).v
    res_el[iv_0]= (v_0 - xdot_0).v
    res_el[iomega_0]= (Omega_0 -W_0).v
    S_el[ix_0,iv_0]=gamma_p*mass*eye(3)
    S_el[iv_0,ix_0]= gamma_p*eye(3)
    A_0 = gamma_p*T_0 + alpha_stiff*DW_0
    S_el[itheta_0,iomega_0]= (gamma_p*Jrot+alpha_stiff*tilde(W_0)*Jrot).mat


end


S_el[iv_0,iv_0]=  - alpha_stiff*eye(3)
S_el[itheta_0,itheta_0]= -(tilde(Jrot*Omega_0)*A_0).mat
S_el[iomega_0,iomega_0]= - alpha_stiff*eye(3)
S_el[iomega_0,itheta_0]  =  A_0.mat

CM_kin_el = 1.0/2.0*(dotp(mass*v_0,v_0)+dotp(Omega_0,Jrot*Omega_0))
pot_el = - transpose(p_el[ix_0])*x_0.v

for ip = 1:2
     X_P[ip] = Vec3([0.5*length*(-1)^ip;  0.0; 0.0 ])
     Xpu = rot(-psi_0, (x_P[ip] - x_0))
     tXpu[ip] = tilde(Xpu)
     phi_P[ip] = RV3(-psi_0, RV3(psi_P[ip],-psi_rel[ip]))
     Tm1P[ip] = invtang(phi_P[ip])
     u[ip] = Xpu - X_P[ip]


     udot_P[ip]= rot(-psi_0,(xdot_P[ip] - xdot_0)) + crossp(Xpu,W_0)
     coef_P = (tilde(udot_P[ip]) + tilde(W_0)*tXpu[ip])*T_0  + tXpu[ip]*DW_0

     if nl_correction == true
         u_P[ip] = transpose(Tm1P[ip])*u[ip]
         DinvT_P = DinvtangT(u_P[ip],phi_P[ip])
#         u_P[ip] = transpose(Tm1P[ip])*Xpu - X_P[ip]
#         DinvT_P = DinvtangT(Xpu,phi_P[ip])
         A_0 = transpose(Tm1P[ip])*tXpu[ip]- DinvT_P*transpose(Tm1P[ip])
         A_P =  DinvT_P*Tm1P[ip]*R_rel[ip]
         B_P = transpose(R_0*Tm1P[ip])
     else
        u_P[ip] = u[ip]
#         u[ip] = Xpu - X_P[ip]
         A_0 = tXpu[ip]
         B_P = transpose(R_0)
     end

     DqDx[ix_P[ip],ix_0] = -B_P.mat
     DqDw[ix_P[ip],ix_0] = DqDx[ix_P[ip],ix_0]
     DqDx[ix_P[ip],itheta_0] = (A_0*T_0).mat
     DqDw[ix_P[ip],itheta_0] = A_0.mat
     DqDx[ix_P[ip],ix_P[ip]] = B_P.mat
     DqDw[ix_P[ip],ix_P[ip]] = DqDx[ix_P[ip],ix_P[ip]]
     if nl_correction == true
         DqDx[ix_P[ip],itheta_P[ip]] = A_P.mat
         DqDw[ix_P[ip],itheta_P[ip]] = (A_P*T_P[ip]).mat
     end
     DqDx[itheta_P[ip],itheta_0] = -(transpose(Tm1P[ip])*T_0).mat
     DqDw[itheta_P[ip],itheta_0] = -(transpose(Tm1P[ip])).mat
     DqDx[itheta_P[ip],itheta_P[ip]] = (Tm1P[ip]*R_rel[ip]*T_P[ip]).mat
     DqDw[itheta_P[ip],itheta_P[ip]] = (Tm1P[ip]*R_rel[ip]).mat

     DqdotDx[ix_P[ip],ix_0] = W0R0T.mat
     DqdotDx[ix_P[ip],itheta_0] = coef_P.mat
     DqdotDx[ix_P[ip],ix_P[ip]] = - W0R0T.mat
     DqdotDx[itheta_P[ip],itheta_0] = - DW_0.mat
     DqdotDx[itheta_P[ip],itheta_P[ip]] =  (R_rel[ip]*DW_P[ip]).mat
end

q_elast = [u_P[1].v; phi_P[1].v; u_P[2].v; phi_P[2].v]
qdot_elast  = [udot_P[1].v; (R_rel[1]*W_P[1] - W_0).v; udot_P[2].v; (R_rel[2]*W_P[2] - W_0).v]

f_elast = K_elast*q_elast

if rotation == true 
    Om = (W_0+Omega_d).v
    G_elast = - (Om[1]*S_elast[1]+Om[2]*S_elast[2]+Om[3]*S_elast[3])
    f_iner = M_elast*vdot_elast + G_elast*v_elast
    A_omega = zeros(12,12)
    A_omega[1:3,1:3] = tOmega_d
    A_omega[7:9,7:9] = tOmega_d
    res_el[iv_P] = v_elast - qdot_elast - A_omega*q_elast
    S_el[ix_0P,iv_P] += transpose(DqDw)*(gamma_p*M_elast + alpha_stiff*G_elast)
    S_el[iv_P,ix_0P] +=  alpha_stiff*DqdotDx + (gamma_p*eye(12)+ alpha_stiff*A_omega)*DqDx
    svq_dot = zeros(3)
    sqv_dot = zeros(3)
    sqv = zeros(3)
    for i = 1:3
        svq_dot[i] = -qdot_elast'*S_elast[i]*v_elast
        sqv_dot[i] = -vdot_elast'*S_elast[i]*q_elast
        sqv[i] = -v_elast'*S_elast[i]*q_elast
    end
    res_el[itheta_0] -=  cross(Om, sqv) + sqv_dot - svq_dot
else 
    f_iner = M_elast*vdot_elast
    res_el[iv_P] = v_elast - qdot_elast
    S_el[ix_0P,iv_P] += gamma_p*transpose(DqDw)*M_elast
    S_el[iv_P,ix_0P] +=  alpha_stiff*DqdotDx + gamma_p*DqDx
end

res_el[ix_0P] -=  transpose(DqDw)*(f_elast+f_iner)
S_el[ix_0P,ix_0P] += alpha_stiff*transpose(DqDw)*K_elast*DqDx
S_el[iv_P,iv_P] = - alpha_stiff*eye(12)


str_el = 1.0/2.0*transpose(q_elast)*f_elast

elast_kin_el = 1.0/2.0*(transpose(v_elast)*M_elast*v_elast)

kin_el = CM_kin_el + elast_kin_el

if matrix == false
    push_element_sparse(res,iel,inv_loc_x,inv_loc_v,S_el,res_el)
    return pot_el, kin_el, str_el
else
    return inv_loc, S_el = condensed_element_matrix(iel,inv_loc_x,inv_loc_v,S_el)

end

end

"""
    super_beam_force

        function constructing the residual vector  of a super-beam element. The element is constructed as a superelement resulting from the assembly of 2 beam elements with cubic interpolation of the deflection. Its stiffness and mass kernel has been construted by Herting's method by the 'super_beam_matrix_kernel.jl' function.

        Calling sequence: super_beam_force(nbr,y,Dy,ydot, res,p)


"""
function super_beam_force(nbr::Int,y::Vector,Dy::Vector,ydot::Vector,res::Vector,p::Vector)

 nc = Main.node_container
 mc =  Main.model_container
 cf =  Main.Frames.current_frames
 ec =  Main.element_container
 sbc = SetElements.super_beam_container

 iel = findfirst(x -> x == nbr, sbc.numbers)[1]


 ndim = 36

 res_el  = Vector{Float64}(zeros(ndim))
 p_el  = Vector{Float64}(zeros(ndim))

 # DqDx = Array{Float64,2}(zeros(12,18))
 DqDw = Array{Float64,2}(zeros(12,18))
 # DqdotDx = Array{Float64,2}(zeros(12,18))

 #
 # Getting detailed characteristics  of the superelement from "superelement_container"
 #

 node_0 = sbc.ref_node_orders[iel]
 bnodes = sbc.boundary_node_orders[iel]
 nodes = [bnodes; node_0]

 length = sbc.length[iel]
 psi_rel = sbc.local_node_orientations[iel]
 R_rel = [rot(psi_rel[i]) for i in 1:2]

 mass = sbc.mass[iel]
 Jrot = Mat3(diagm(sbc.Jrot[iel]))

 
 rotation = mc.uniform_rotation
 rotation == true && (Omega_d = mc.rotation_speed; tOmega_d = tilde(Omega_d).mat; 
 S_elast = sbc.S[iel])

if rotation == true
    time = itime*h
    phi_d = RV3(time*Omega_d)
    grav = rot(-phi_d, mc.gravity)
else
    grav = mc.gravity
end

 K_elast = sbc.K_elast[iel]
 M_elast = sbc.M_elast[iel]
 

 nl_correction = sbc.nl_correction[iel]

 #
 # getting localization of the superelement from "element_container"
 #

 iel = findfirst(x -> x == nbr, ec.element_numbers)[1]
 inv_loc_x = ec.inv_loc_x[iel]
 inv_loc_v = ec.inv_loc_v[iel]

 inv_loc =[inv_loc_x; inv_loc_v.+ mc.Ndofs_q]

 #
 # pull current solution at element level
 #
 (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

 ix_P = Vector{Vector{Int}}(undef,2)
 x_P = Vector{Vec3}(undef,2)
 X_P = Vector{Vec3}(undef,2)
 Dx_P = Vector{Vec3}(undef,2)
 xdot_P = Vector{Vec3}(undef,2)
 itheta_P = Vector{Vector{Int}}(undef,2)
 psi_P = Vector{RV3}(undef,2)
 Dpsi_P = Vector{RV3}(undef,2)
 psidot_P = Vector{Vec3}(undef,2)

 R_P = Vector{Mat3}(undef,2)
 T_P = Vector{Mat3}(undef,2)
 W_P = Vector{Vec3}(undef,2)
 DW_P = Vector{Mat3}(undef,2)
 u_P = Vector{Vec3}(undef,2)
 u = Vector{Vec3}(undef,2)
 udot_P = Vector{Vec3}(undef,2)
 phi_P = Vector{RV3}(undef,2)
 Tm1P = Vector{Mat3}(undef,2)
 tXpu = Vector{Mat3}(undef,2)

 for ip = 1:2

     ix_P[ip] = [(ip-1)*6+i for i in 1:3]
     itheta_P[ip] = ix_P[ip] .+ 3
     xdot_P[ip] = Vec3(sdot[ix_P[ip]])
     psidot_P[ip] = Vec3(sdot[itheta_P[ip]])
     Dx_P[ip] = Vec3(Ds[ix_P[ip]])
     Dpsi_P[ip] = RV3(Ds[itheta_P[ip]])
     psi_P[ip]= RV3(cf[nodes[ip]].psi,Dpsi_P[ip])
     x_P[ip]= cf[nodes[ip]].x + Dx_P[ip]
     R_P[ip] = rot(psi_P[ip])
     T_P[ip] = tang(Dpsi_P[ip])
     W_P[ip] = T_P[ip]*psidot_P[ip]
     DW_P[ip] = Dtang(psidot_P[ip],Dpsi_P[ip])

 end

ix_0 = ix_P[1] .+ 12
iv_0 = ix_0 .+ 18
itheta_0 = itheta_P[1] .+ 12
iomega_0 = itheta_0 .+18

Dx_0 = Vec3(Ds[ix_0])
Dpsi_0 = RV3(Ds[itheta_0])
x_0 = cf[nodes[3]].x + Dx_0
psi_0 = RV3(cf[nodes[3]].psi,Dpsi_0)
xdot_0 = Vec3(sdot[ix_0])
psidot_0 = Vec3(sdot[itheta_0])
R_0 = rot(psi_0)
T_0 = tang(Dpsi_0)
W_0 = T_0*psidot_0
DW_0 = Dtang(psidot_0,Dpsi_0)

# W0R0T = tilde(W_0)*transpose(R_0)



v_0 = Vec3(s[iv_0]+Ds[iv_0])
vdot_0 = Vec3(sdot[iv_0])
Omega_0 = Vec3(s[iomega_0]+Ds[iomega_0])
Omegadot_0 = Vec3(sdot[iomega_0])

iq_P = [ix_P[1]; itheta_P[1]; ix_P[2]; itheta_P[2]]
iv_P = iq_P .+ 18
ix_0P = [iq_P; ix_0; itheta_0]
v_elast = s[iv_P] + Ds[iv_P]
vdot_elast = sdot[iv_P]

p_el[ix_0] = (mass*grav).v

if rotation == true 
    ROmega_d = rot(-psi_0,Omega_d)
    res_el[ix_0] =  (mass*(grav - vdot_0 - crossp(Omega_d,v_0))).v
    res_el[iv_0]= (v_0 - xdot_0 - crossp(Omega_d,x_0)).v 
    res_el[itheta_0] = -(Jrot*Omegadot_0 + crossp(W_0+ROmega_d, Jrot*Omega_0)).v
    res_el[iomega_0]= (Omega_0 - W_0 - ROmega_d).v
    A_0= (gamma_p*I3 + alpha_stiff*tilde(ROmega_d))*T_0 + alpha_stiff*DW_0
else 
    res_el[ix_0] =  mass*(grav - vdot_0).v
    res_el[itheta_0] = -(Jrot*Omegadot_0 + crossp(W_0,Jrot*Omega_0)).v
    res_el[iv_0]= (v_0 - xdot_0).v
    res_el[iomega_0]= (Omega_0 -W_0).v
end

CM_kin_el = 1.0/2.0*(dotp(mass*v_0,v_0)+dotp(Omega_0,Jrot*Omega_0))
pot_el = - transpose(p_el[ix_0])*x_0.v

for ip = 1:2
     X_P[ip] = Vec3([0.5*length*(-1)^ip;  0.0; 0.0 ])
     Xpu = rot(-psi_0, (x_P[ip] - x_0))
     tXpu[ip] = tilde(Xpu)
     phi_P[ip] = RV3(-psi_0, RV3(psi_P[ip],-psi_rel[ip]))
     Tm1P[ip] = invtang(phi_P[ip])
     u[ip] = Xpu - X_P[ip]


     udot_P[ip]= rot(-psi_0,(xdot_P[ip] - xdot_0)) + crossp(Xpu,W_0)
#     coef_P = (tilde(udot_P[ip]) + tilde(W_0)*tXpu[ip])*T_0  + tXpu[ip]*DW_0

     if nl_correction == true
         u_P[ip] = transpose(Tm1P[ip])*u[ip]
         DinvT_P = DinvtangT(u_P[ip],phi_P[ip])
         A_0 = transpose(Tm1P[ip])*tXpu[ip]- DinvT_P*transpose(Tm1P[ip])
         A_P =  DinvT_P*Tm1P[ip]*R_rel[ip]
         B_P = transpose(R_0*Tm1P[ip])
     else
        u_P[ip] = u[ip]
         A_0 = tXpu[ip]
         B_P = transpose(R_0)
     end

     DqDw[ix_P[ip],ix_0] = -B_P.mat
     DqDw[ix_P[ip],itheta_0] = A_0.mat
     DqDw[ix_P[ip],ix_P[ip]] = B_P.mat
     if nl_correction == true
         DqDw[ix_P[ip],itheta_P[ip]] = (A_P*T_P[ip]).mat
     end
     DqDw[itheta_P[ip],itheta_0] = -(transpose(Tm1P[ip])).mat
     DqDw[itheta_P[ip],itheta_P[ip]] = (Tm1P[ip]*R_rel[ip]).mat
    end
q_elast = [u_P[1].v; phi_P[1].v; u_P[2].v; phi_P[2].v]
qdot_elast  = [udot_P[1].v; (R_rel[1]*W_P[1] - W_0).v; udot_P[2].v; (R_rel[2]*W_P[2] - W_0).v]

f_elast = K_elast*q_elast

if rotation == true 
    Om = (W_0+Omega_d).v
    G_elast = - (Om[1]*S_elast[1]+Om[2]*S_elast[2]+Om[3]*S_elast[3])
    f_iner = M_elast*vdot_elast + G_elast*v_elast
    A_omega = zeros(12,12)
    A_omega[1:3,1:3] = tOmega_d
    A_omega[7:9,7:9] = tOmega_d
    res_el[iv_P] = v_elast - qdot_elast - A_omega*q_elast
    svq_dot = zeros(3)
    sqv_dot = zeros(3)
    sqv = zeros(3)
    for i = 1:3
        svq_dot[i] = -qdot_elast'*S_elast[i]*v_elast
        sqv_dot[i] = -vdot_elast'*S_elast[i]*q_elast
        sqv[i] = -v_elast'*S_elast[i]*q_elast
    end
    res_el[itheta_0] -=  cross(Om, sqv) + sqv_dot - svq_dot
else 
    f_iner = M_elast*vdot_elast
    res_el[iv_P] = v_elast - qdot_elast
end

res_el[ix_0P] -=  transpose(DqDw)*(f_elast+f_iner)

str_el = 1.0/2.0*transpose(q_elast)*f_elast

elast_kin_el = 1.0/2.0*(transpose(v_elast)*M_elast*v_elast)

kin_el = CM_kin_el + elast_kin_el

for i = 1:ndim
    iloc = inv_loc[i] 
    iloc > 0  &&  (res[iloc] += res_el[i]) 
end

pot_el = 0.0

return kin_el, str_el, pot_el

end


