"""
    shell

Function constructing the residual vector and the iteration matrix of a shell element. The material can be elastic or viscoeslatic. It is then described by a viscoeslatic model of Mawell or Kelvin-Voigt types.

The material type is governed by the `visco_type` parameter provided by the `shell_container` array.
        
* `visco_type` = "none"       elastic material
* `visco_type` = "maxwell"    Maxwell viscoelastic material
* `visco_type` = "damped"   Kelvin-Voigt viscoelastic material

The _Maxwell_ material model is described as 
    
> ``\\sigma`` =  ``K_{\\infty} . \\epsilon + K_{b} . \\tau . \\dot\\alpha``
>
> ``\\alpha  + \\tau \\dot \\alpha = \\epsilon``
>
> with    
>        
> ``K_{\\infty} = r_{\\infty} K ``   and     ``K_{b} = (1.0 - r_{\\infty})K``
>
> and with the vector of time_constants defined as
> 
> `` \\tau  = [\\tau_{E}, \\, \\tau_S, \\, \\tau_S ] ``
>
> with
>
>``\\tau_{E} = \\frac{1}{3}/((1-2\\nu)\\tau_B + 2(1+\\nu)\\tau_S)``
>
> and
> 
> ``\\tau_B``  and ``\\tau_{S}`` being the time constants of the 3D-material in bulk deformation and  shear.  
        
The _Kelvin-Voigt_ material model is described as 
    
>
> ``\\sigma`` =  ``K_{\\infty} . \\epsilon + K_{b}  . \\dot\\epsilon``
>
           
Calling sequence:
>
> `kin_el`, `str_el`, `pot_el` = `shell_force`(nbr,`y_n`,Dy,`ydot_np1`,res,p,`alpha_stiff`,`theta_p`,matrix,itime,niter,h)
>
"""
function shell(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64}, alpha_stiff::Float64,theta_p::Float64,matrix::Bool,itime::Int,niter::Int,h::Float64)

    # nbc = Main.node_container
    cf =  Main.Frames.current_frames
    
    mc =  Main.model_container

    sc = Main.SetElements.shell_container
    ec = Main.element_container    

    iel = findfirst(x -> x == nbr, sc.numbers)[1]
    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]


    
    Nnodes = sc.npoints[iel]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]
    N_depl = Nnodes*6
    Nv = sum(inv_loc_v)
    N_vit = Nnodes*5
    Ndofs_q = mc.Ndofs_q
    Nv == 0 ?  (Ndim = N_depl; inv_loc= inv_loc_x) : (Ndim = N_depl+N_vit; inv_loc= [inv_loc_x; inv_loc_v.+Ndofs_q])
   


    global tangent_stiffness = zeros(N_depl,N_depl)
#   global M_NQ = zeros(N_vit,N_depl)
    global S_el = zeros(Ndim,Ndim)
    global res_el  = zeros(Ndim)
    global p_el  = zeros(Ndim)
    global f_int  = zeros(N_depl)
    global f_iner  = zeros(N_depl)
    global w  = zeros(N_depl)
    global str_el = 0.0
    global stresses = zeros(12)

    visco_type = sc.visco_type[iel]



    ngauss_points = sc.ngauss_points[iel]
    R_rel = sc.relative_rotations[iel]
    f_0g = sc.init_deformations[iel]
    K = sc.stiffness_matrix[iel]
    M_NN = sc.mass_matrix[iel]
    F = sc.shapes[iel]
    DF = sc.shape_derivatives[iel]
    J = sc.gauss_jacobians[iel]
    wg = sc.gauss_weights[iel]

    if visco_type == "maxwell"
        time_constants = sc.time_constants[iel]
        ratio_infty = sc.ratio_infty[iel]
        if itime  > 1
            Gamma_2 = broadcast(x -> h/x/2.0*exp(-h/2.0/x), time_constants)
            Gamma_1 = broadcast(x -> exp(-h/x), time_constants)
            # Gamma_1 = broadcast(x -> (2.0*x/h-1.0)/(1.0+2.0*x/h), time_constants)
            # Gamma_2 =  broadcast(x -> 1.0/(1.0+2.0*x/h), time_constants)
        else
            Gamma_1 = zeros(12)
            Gamma_2 = zeros(12)
        end
        if niter == 1 
            for ig = 1:ngauss_points
                sc.strains_g[iel][ig][:,1] = sc.strains_g[iel][ig][:,2]
                sc.visco_strains_g[iel][ig][:,1] = sc.visco_strains_g[iel][ig][:,2]
            end
        end
    end



    global penalty  = 0.0
    global psi_12 = zeros(12)

    if Nnodes > 4
        penalty  = sc.stiffness_properties[iel][1]*sc.thickness[iel]^3/12
        psi_12[6] = 1.0
        psi_12[12] = 1.0
    end



    
    global H = Vector{NodeFrame}(undef,Nnodes)
    global H_n = Vector{NodeFrame}(undef,Nnodes)
    global f_g = Vector{Array{Float64,2}}(undef,ngauss_points)
    global eps_g = Vector{Array{Float64,2}}(undef,ngauss_points)


    inode = sc.node_orders[iel]


   (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)
   # (s,sdot,Ds) = pull_vectors(inv_loc_v.+Ndofs_q,y,ydot,Dy)
   # (v,vdot,Dv) = pull_vectors(inv_loc_v.+Ndofs_q,y,ydot,Dy)

    global inx = []
    for k = 1:Nnodes
        inx = [inx; [(k-1)*6+i for i in 1:5]]
    end

    str_el = 0.0

    R_n = Vector{Matrix}(undef, Nnodes)
    invTt = Vector{Matrix}(undef, Nnodes)


    S_r = zeros(N_depl,N_depl)
    S_l = zeros(N_depl,N_depl)
    visco_type == "damped" && (S_rd = zeros(N_depl,N_depl))
    x = Vector{Vec3}(undef,Nnodes)
    psi = Vector{RV3}(undef,Nnodes)

    idepl = [i for i in 1:N_depl]

    visco_type == "damped" && (
        global tangent_damping = zeros(N_depl,N_depl);
        global time_constants = diagm(sc.time_constants[iel]);
        global qdot = sdot[idepl])

    for ip = 1:Nnodes
        ix = [(ip-1)*6+i for i in 1:3] 
        irot = ix .+ 3
        x_n = cf[inode[ip]].x
        Dx = Vec3(Ds[ix])
        x[ip] = x_n + Dx
        psi_n = cf[inode[ip]].psi
        Dpsi = RV3(Ds[irot])
        psi[ip] = RV3(psi_n, Dpsi)
       
        psi_p = RV3(psi[ip],-R_rel[ip])
        
        H[ip] = NodeFrame(x[ip], psi_p)

        S_r[ix,ix] =  (transpose(invtang(Dpsi))*rot(RV3(R_rel[ip],-psi_n))).mat
        S_r[irot,irot] = (rot(R_rel[ip])*tang(Dpsi)).mat 
        S_l[ix,ix] = rot(-psi_p).mat
        S_l[irot,irot] = rot(R_rel[ip]).mat


        visco_type == "damped" && (S_rd[ix,ix] = rot(-psi_p).mat; 
                                    S_rd[irot,irot] = S_r[irot,irot])

    end

    strains = zeros(12)

    if visco_type == "maxwell"

        for ing = 1:ngauss_points

            f_g, QN = frame_solve(H, F[ing], DF[ing], J[ing], Nnodes)
            eps_g = reshape((f_g - f_0g[ing]),(12,1))
            sc.visco_strains_g[iel][ing][:,2] = Gamma_2 .* (eps_g + sc.strains_g[iel][ing][:,1]) + Gamma_1 .* sc.visco_strains_g[iel][ing][:,1] 
            strains += eps_g*wg[ing]
            sc.strains_g[iel][ing][:,2] = eps_g*wg[ing]
            d_stress = K[ing]*eps_g - (1.0 - ratio_infty) *K[ing]*sc.visco_strains_g[iel][ing][:,2]
            stresses += d_stress
            f_int  += transpose(QN)*(d_stress + penalty*psi_12*psi_12'*eps_g)
            tangent_stiffness  += transpose(QN)*(K[ing] -(1.0 -ratio_infty) * Gamma_2 .* K[ing]  + penalty*psi_12*psi_12')*QN
        end

    else

        for ing = 1:ngauss_points

            f_g, QN = frame_solve(H, F[ing], DF[ing], J[ing], Nnodes)
            eps_g = reshape((f_g - f_0g[ing]),(12,1))
            str_el += 0.5*(transpose(eps_g)*K[ing]*eps_g)[1]
            strains += eps_g*wg[ing]
            stresses += K[ing]*eps_g
            f_int  += transpose(QN)*(K[ing]+penalty*psi_12*psi_12')*eps_g
            tangent_stiffness  += transpose(QN)*(K[ing]+penalty*psi_12*psi_12')*QN
            visco_type == "damped" && (tangent_damping  += transpose(QN)*time_constants*K[ing]*QN)
        end

    end
     
    sc.stresses[iel][:] = stresses/sc.area[iel]
    sc.strains[iel][:,2] = strains


    tangent_stiffness = S_l'*tangent_stiffness*S_r

    f_int = S_l'*f_int

    # gravity = (mc.gravity).v 

    res_el[idepl] = - f_int
    S_el[idepl,idepl] = alpha_stiff*tangent_stiffness

    visco_type == "damped" && (
        tangent_damping = S_l'*tangent_damping*S_rd;
        res_el[idepl] -=  tangent_damping*qdot;
        S_el[idepl,idepl] += theta_p*tangent_damping)


    kin_el = 0.0    

if  Nv > 0

    iv_x = reshape([(i-1)*6+k for k in 1:3, i in 1:Nnodes], (Nnodes*3,1))
    ivit = [Nnodes*6+i for i in 1:Nnodes*5]

    
    vdot = sdot[ivit]
    v = s[ivit] + Ds[ivit]


    Iom = eye(3)[1:2,1:3]

    mv = M_NN[:,inx]*v

    kin_el = 0.5*v'*M_NN[inx,inx]*v

    mvdot = M_NN[:,inx]*vdot

    global R_b = eye(6)
    
    for np = 1:Nnodes

        ix = [(np-1)*6+i for i in 1:3]
        iv_x = [N_depl+(np-1)*5+i for i in 1:3]
        iq = [(np-1)*6+i for i in 1:6]
        itheta = ix .+ 3
        iomega = [N_depl+(np-1)*5+i for i in 4:5]
        iv = [iv_x;iomega]
        inp = inode[np]

        psin = RV3(cf[inp].psi,-R_rel[np])
        Dpsi = RV3(Ds[itheta])
        psi = RV3(psin,Dpsi)
        R_p = rot(psi)

        
        psidot = Vec3(sdot[itheta])
        xdot = Vec3(sdot[ix])

        W = tang(Dpsi,psidot)
        T = tang(Dpsi)
        v_x = rot(-psi,xdot)
        w[ix] = v_x.v
        w[itheta] = W.v
        DW = Dtang(psidot,Dpsi)

        # W_p = [(tilde(W)).mat  zeros(3,3); zeros(3,3) zeros(3,3)]
        Z_p = [(R_p*tilde(W)).mat  zeros(3,3); tilde(v_x).mat tilde(W).mat]

        R_b[1:3,1:3] = R_p.mat
        R_b[4:6,4:6] = rot(-R_rel[np]).mat

        S_el[iv_x, ix] = theta_p*R_p.mat' 
        S_el[iomega, itheta] = Iom*(theta_p*T + alpha_stiff*DW).mat
        S_el[iv_x, itheta] =   alpha_stiff*(tilde(v_x)*T).mat

       
        # f_iner[iq]  = R_b*(mvdot[iq] + W_p*mv[iq])

        f_iner[iq]  = R_b*mvdot[iq] + Z_p*mv[iq]

        
        S_el[iq,ivit] = R_b*(theta_p*M_NN[iq, inx] + alpha_stiff*Z_p*M_NN[iq, inx])

        # u_p = mvdot[ix]
        # S_el[ix,itheta] -= (R_p*tilde(Vec3(u_p))*T).mat
    end
        S_el[ivit,ivit] = -alpha_stiff*eye(N_vit)
        res_el[ivit] = v - w[inx]
        res_el[idepl] -= f_iner
end


if matrix == false

    pot_el = 0.0


    push_element_sparse(res,p, iel2,inv_loc_x,inv_loc_v,S_el,res_el, p_el)

    return kin_el, str_el, pot_el
else
    return inv_loc, S_el = condensed_element_matrix(iel2,inv_loc_x,inv_loc_v,S_el)
end

end


"""
    shell_force

Function constructing the residual vector  of an elastic shell element. The material can be elastic or viscoeslatic. It is then described by a viscoeslatic model of Mawell or Kelvin-Voigt types as explained for the function `shell`.
        
Calling sequence:
>
> `kin_el`, `str_el`, `pot_el` = `shell_force`(nbr,`y_n`,Dy,`ydot_np1`,res,p,itime,niter,h)
>


        
"""
function shell_force(nbr::Int,y::Vector{Float64},Dy::Vector{Float64},ydot::Vector{Float64},
    res::Vector{Float64},p::Vector{Float64},itime::Int,niter::Int,h::Float64)

    cf =  Main.Frames.current_frames
    
    mc =  Main.model_container

    sc = Main.SetElements.shell_container
    ec = Main.element_container    

    iel = findfirst(x -> x == nbr, sc.numbers)[1]
    iel2 = findfirst(x -> x == nbr, ec.element_numbers)[1]


    
    Nnodes = sc.npoints[iel]

    inv_loc_x= ec.inv_loc_x[iel2]
    inv_loc_v = ec.inv_loc_v[iel2]
    N_depl = Nnodes*6
    Nv = sum(inv_loc_v)
    N_vit = Nnodes*5
    Ndofs_q = mc.Ndofs_q
    Nv == 0 ?  (Ndim = N_depl; inv_loc= inv_loc_x) : (Ndim = N_depl+N_vit; inv_loc= [inv_loc_x; inv_loc_v.+Ndofs_q])


    global res_el  = zeros(Ndim)
    global p_el  = zeros(Ndim)
    global f_int  = zeros(N_depl)
    global f_iner  = zeros(N_depl)
    global w  = zeros(N_depl)
    global str_el = 0.0
    global stresses = zeros(12)

    ngauss_points = sc.ngauss_points[iel]
    R_rel = sc.relative_rotations[iel]
    f_0g = sc.init_deformations[iel]
    K = sc.stiffness_matrix[iel]
    M_NN = sc.mass_matrix[iel]
    F = sc.shapes[iel]
    DF = sc.shape_derivatives[iel]
    J = sc.gauss_jacobians[iel]
    wg = sc.gauss_weights[iel]
    visco_type = sc.visco_type[iel]

    if visco_type == "maxwell"
        time_constants = sc.time_constants[iel]
        ratio_infty = sc.ratio_infty[iel]
        if itime  > 1
            Gamma_2 = broadcast(x -> h/x/2.0*exp(-h/2.0/x), time_constants)
            Gamma_1 = broadcast(x -> exp(-h/x), time_constants)
            # Gamma_1 = broadcast(x -> (2.0*x/h-1.0)/(1.0+2.0*x/h), time_constants)
            # Gamma_2 =  broadcast(x -> 1.0/(1.0+2.0*x/h), time_constants)
        else
            Gamma_1 = zeros(12)
            Gamma_2 = zeros(12)
        end
    end

    global penalty  = 0.0
    global psi_12 = zeros(12)

    if Nnodes > 4
        penalty  = sc.stiffness_properties[iel][1]*sc.thickness[iel]^3/12
        psi_12[6] = 1.0
        psi_12[12] = 1.0
    end

    
    global H = Vector{NodeFrame}(undef,Nnodes)
    global H_n = Vector{NodeFrame}(undef,Nnodes)
    global f_g = Vector{Array{Float64,2}}(undef,ngauss_points)
    global eps_g = Vector{Array{Float64,2}}(undef,ngauss_points)


    inode = sc.node_orders[iel]


    (s,sdot,Ds) = pull_vectors(inv_loc,y,ydot,Dy)

    global inx = []
    for k = 1:Nnodes
        inx = [inx; [(k-1)*6+i for i in 1:5]]
    end

    str_el = 0.0

    S_l = eye(N_depl)

    x = Vector{Vec3}(undef,Nnodes)
    psi = Vector{RV3}(undef,Nnodes)


    for ip = 1:Nnodes
        ix = [(ip-1)*6+i for i in 1:3] 
        irot = ix .+ 3
        x_n = cf[inode[ip]].x
        Dx = Vec3(Ds[ix])
        x[ip] = x_n + Dx
        psi_n = cf[inode[ip]].psi
        Dpsi = RV3(Ds[irot])
        psi[ip] = RV3(psi_n, Dpsi)
       
        psi_p = RV3(psi[ip],-R_rel[ip]) 
        
        H[ip] = NodeFrame(x[ip], psi_p)

   
        S_l[ix,ix] = rot(RV3(R_rel[ip],-psi[ip])).mat
        S_l[irot,irot] = rot(R_rel[ip]).mat

    end

    strains = zeros(12)

    if visco_type == "maxwell"

        for ing = 1:ngauss_points

            f_g, QN = frame_solve(H, F[ing], DF[ing], J[ing], Nnodes)
            eps_g = reshape((f_g - f_0g[ing]),(12,1))
            sc.visco_strains_g[iel][ing][:,2] = Gamma_2 .* (eps_g + sc.strains_g[iel][ing][:,1]) + Gamma_1 .* sc.visco_strains_g[iel][ing][:,1] 
            strains += eps_g*wg[ing]
            d_stress = K[ing]*eps_g - (1.0 - ratio_infty) *K[ing]*sc.visco_strains_g[iel][ing][:,2]
            stresses += d_stress
            str_el += 0.5*(transpose(eps_g)*d_stress)[1]
            f_int  += transpose(QN)*(d_stress + penalty*psi_12*psi_12'*eps_g)

        end
    else

        for ing = 1:ngauss_points

            f_g, QN = frame_solve(H, F[ing], DF[ing], J[ing], Nnodes)
            eps_g = reshape((f_g - f_0g[ing]),(12,1))
            str_el += 0.5*(transpose(eps_g)*K[ing]*eps_g)[1]
            stresses += K[ing]*eps_g
            strains += wg[ing]*eps_g
            f_int  += transpose(QN)*(K[ing]+penalty*psi_12*psi_12')*eps_g

        end
    end
   
     
    sc.stresses[iel][:] = stresses/sc.area[iel]
    sc.strains[iel][:,2] = strains

    f_int = S_l'*f_int



    # gravity = (mc.gravity).v 



    
    idepl = [i for i in 1:N_depl]


    res_el[idepl] = - f_int
    
    kin_el = 0.0


    if  Nv > 0

        iv_x = reshape([(i-1)*6+k for k in 1:3, i in 1:Nnodes], (Nnodes*3,1))
        ivit = [Nnodes*6+i for i in 1:Nnodes*5]
    
        
        vdot = sdot[ivit]
        v = s[ivit] + Ds[ivit]
    
    
        Iom = eye(3)[1:2,1:3]
    
        mv = M_NN[:,inx]*v

        kin_el = 0.5*v'*mv[ivit .- N_depl]
    
        mvdot = M_NN[:,inx]*vdot
    
        global R_b = eye(6)
        
        for np = 1:Nnodes
    
            ix = [(np-1)*6+i for i in 1:3]
            iv_x = [N_depl+(np-1)*5+i for i in 1:3]
            iq = [(np-1)*6+i for i in 1:6]
            itheta = ix .+ 3
            inp = inode[np]
    
            psin = RV3(cf[inp].psi,-R_rel[np])
            Dpsi = RV3(Ds[itheta])
            R_p = rot(psin)
    
            
            psidot = Vec3(sdot[itheta])
            xdot = Vec3(sdot[ix])
    
            W = tang(Dpsi,psidot)
            v_x = rot(-psin,xdot)
            w[ix] = v_x.v
            w[itheta] = W.v

    
    
            # W_p = [(tilde(W)).mat  zeros(3,3); zeros(3,3) zeros(3,3)]

            Z_p = [(R_p*tilde(W)).mat  zeros(3,3); tilde(v_x).mat tilde(W).mat]
    
            R_b[1:3,1:3] = R_p.mat
            R_b[4:6,4:6] = rot(-R_rel[np]).mat
    
          
            # f_iner[iq]  = R_b*(mvdot[iq] + W_p*mv[iq])

            f_iner[iq]  = R_b*mvdot[iq] + Z_p*mv[iq]
               
        end

            res_el[ivit] = v - w[inx]
            res_el[idepl] -= f_iner

end  


for i = 1:Ndim
    iloc = inv_loc[i] 
    iloc > 0  &&  (res[iloc] += res_el[i]) 
end

pot_el = 0.0

return kin_el, str_el, pot_el

end

