"""
    record_on_h5_file_dynamic
"""
function record_on_h5_file_dynamic(dsets::Vector{Any},itime::Int,itime_saved::Int,itime_vals_saved::Int,times::Float64,
    niter::Int,y_n::Vector{Float64},ydot::Vector{Float64},
    p::Vector{Float64},pot_energy::Float64,kin_energy::Float64,
    str_energy::Float64,ext_work::Float64,Npas::Int,save::Bool,save_freq::Int,
    eigvals::Bool,eig_freq::Int,vals::Vector)

    mc = Main.model_container

    
    Nnodes = mc.Nodes
    N_shells = mc.Shells
    N_beams = mc.Beams
    
    N_shells > 0 && (sc = Main.SetElements.shell_container)
    N_beams > 0 && (bc = Main.SetElements.beam_container)
    
    if save == true &&  (itime == 1 || mod(itime,save_freq) == 0 || itime == Npas)

        itime_saved += 1
        
        dsets[1][itime_saved,1] = Float64(times)
        dsets[4][itime_saved,1] = niter
        dsets[2][itime_saved,:] = y_n
        dsets[3][itime_saved,:] = ydot
        dsets[5][itime_saved,:] = p
        dsets[6][itime_saved,1] = kin_energy
        dsets[7][itime_saved,1] = pot_energy
        dsets[8][itime_saved,1] = str_energy
        dsets[9][itime_saved,1] = ext_work

        for i = 1:Nnodes
            dsets[13][itime_saved,1:3,i] = (Main.Frames.current_frames[i].x).v
        end

        N_shells > 0 && (for iel in 1:N_shells dsets[16][itime_saved,:,iel] = sc.stresses[iel] end)
        N_beams > 0 && (for iel in 1:N_beams dsets[18][itime_saved,:,iel] = bc.stresses[iel] end)
        N_shells > 0 && (for iel in 1:N_shells dsets[24][itime_saved,:,iel] = sc.strains[iel][:,2] end)
        N_beams > 0 && (for iel in 1:N_beams dsets[25][itime_saved,:,iel] = bc.strains[iel][:,2] end)
    
    end

    if eigvals == true &&  (mod(itime,eig_freq) == 0) && (itime > 1)

        itime_vals_saved += 1

        dsets[10][itime_vals_saved,1] = times
        dsets[11][itime_vals_saved,:] = real(vals)
        dsets[12][itime_vals_saved,:] = imag(vals)


    end

    return itime_saved, itime_vals_saved

end

