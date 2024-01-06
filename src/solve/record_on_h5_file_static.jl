
"""
    record_on_h5_file_static
"""
function record_on_h5_file_static(dsets::Vector{Any},itime::Int,itime_vals_saved::Int,times::Float64,niter::Int,y_n::Vector{Float64},
    ydot::Vector{Float64},p::Vector{Float64},pot_energy::Float64,str_energy::Float64,ext_work::Float64,
    eigvals::Bool,eig_freq::Int,vals::Vector,lambda_n::Vector{Float64},bounds::Vector{Float64})

    mc = Main.model_container

    
    Nnodes = mc.Nodes
    Nbounds = mc.Inequalities
    N_shells = mc.Shells
    N_beams = mc.Beams
  
    N_shells > 0 && (sc = Main.SetElements.shell_container)
    N_beams > 0 && (bc = Main.SetElements.beam_container)
    
    dsets[1][itime,1] = Float64(times)
    dsets[4][itime,1] = niter
    dsets[2][itime,:] = y_n
    dsets[3][itime,:] = ydot
    dsets[5][itime,:] = p
    dsets[7][itime,1] = pot_energy
    dsets[8][itime,1] = str_energy
    dsets[9][itime,1] = ext_work

    for i = 1:Nnodes
        dsets[13][itime,1:3,i] = (Main.Frames.current_frames[i].x).v
    end

    Nbounds > 0 && (dsets[14][itime,:] = lambda_n; dsets[15][itime,:] = bounds)
    N_shells > 0 && (for iel in 1:N_shells dsets[16][itime,:,iel] = sc.stresses[iel] end)
    N_beams > 0 && (for iel in 1:N_beams dsets[18][itime,:,iel] = bc.stresses[iel] end)
    N_shells > 0 && (for iel in 1:N_shells dsets[24][itime,:,iel] = sc.strains[iel][:,2] end)
    N_beams > 0 && (for iel in 1:N_beams dsets[25][itime,:,iel] = bc.strains[iel][:,2] end)
    
    if eigvals == true &&  (mod(itime,eig_freq) == 0) && (itime > 1)

        dsets[10][itime_vals_saved,1] = times
        dsets[11][itime_vals_saved,:] = real(vals)
        dsets[12][itime_vals_saved,:] = imag(vals)


    end

    return

end

