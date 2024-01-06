"""
    create_h5_file_dynamic

        Function creating the .h5 file for archiving the results of a dynamic analysis

        Input:
            * h5_file::String       file name
            * Npas::Int             number of time steps
            * save_freq::Int        archiving frequency
            * eig_vals::Bool        controls eigenvalue analysis
            * eig_freq::Int         frequency of eigenvalue analysis
            * max_vals::Int         number of computed eigenvaluses

        Calling sequence:
        
            create_h5_file_dynamic(h5_file,Npas,save_freq,eigvals,eig_freq,max_vals)

"""
function create_h5_file_dynamic(h5_file::String, Npas::Int,save_freq::Int,eigvals::Bool,eig_freq::Int,max_vals::Int)

mc = Main.model_container
mc.Shells > 0 && (sc = Main.SetElements.shell_container)
mc.Beams > 0 && (bc = Main.SetElements.beam_container)

Ndofs = mc.Ndofs
Nel = mc.Elements
Ndofs_q = mc.Ndofs_q
Ndofs_x = mc.Ndofs_x
Ndofs_int = mc.Ndofs_int
Ndofs_mult = mc.Ndofs_mult
Ndofs_v = mc.Ndofs_v
Nnodes = mc.Nodes
Nel = mc.Elements
Nbounds = mc.Inequalities
N_shells = mc.Shells
N_beams = mc.Beams

N_saved_steps = Int(floor(Npas/save_freq))+2
mod(Npas, save_freq) == 0 &&  (N_saved_steps -=1)

dsets = Vector{Any}(undef,30)
model_data = Vector{Int}(undef,20)

file=h5open(h5_file,"w")


g1 = create_group(file,"times")
g2 = create_group(file,"y")
g3 = create_group(file,"ydot")
g4 = create_group(file,"nitmax")
g5 = create_group(file,"p")
g6 = create_group(file,"kin_energy")
g7 = create_group(file,"pot_energy")
g8 = create_group(file,"str_energy")
g9 = create_group(file,"ext_work")
g13 = create_group(file,"node_positions")
g30 = create_group(file,"model_data")

dsets[1] = create_dataset(g1,"times",datatype(Float64), dataspace(N_saved_steps,1))
dsets[2] = create_dataset(g2,"y",datatype(Float64), dataspace(N_saved_steps,Ndofs))
dsets[3] = create_dataset(g3,"ydot",datatype(Float64), dataspace(N_saved_steps,Ndofs))
dsets[4] = create_dataset(g4,"nitmax",datatype(Int), dataspace(N_saved_steps,1))
dsets[5] = create_dataset(g5,"p",datatype(Float64), dataspace(N_saved_steps,Ndofs))
dsets[6] = create_dataset(g6,"kin_energy",datatype(Float64), dataspace(N_saved_steps,1))
dsets[7] = create_dataset(g7,"pot_energy",datatype(Float64), dataspace(N_saved_steps,1))
dsets[8] = create_dataset(g8,"str_energy",datatype(Float64), dataspace(N_saved_steps,1))
dsets[9] = create_dataset(g9,"ext_work",datatype(Float64), dataspace(N_saved_steps,1))
dsets[13] = create_dataset(g13,"node_positions",datatype(Float64), dataspace(N_saved_steps,3,Nnodes))
dsets[30] = create_dataset(g30,"model_data",datatype(Int), dataspace(20,1))



if eigvals == true
    mc.eigvals = eigvals
    N_eigvals_steps = Int(floor(Npas/eig_freq))+2
    mod(Npas, eig_freq) == 0 &&  (N_eigvals_steps -=1)
    max_vals > Ndofs_q && (max_vals = Ndofs_q)
    max_vals > 20 && (max_vals = 20)
    println("number of computed eigenvalues: ", max_vals)
    println("number of eigenvalue records: ", N_eigvals_steps)
    g10 = create_group(file,"time_vals")
    g11 = create_group(file,"re_vals")
    g12 = create_group(file,"im_vals")
    dsets[10] = create_dataset(g10,"time_vals",datatype(Float64), dataspace(N_eigvals_steps,1))
    dsets[11] = create_dataset(g11,"re_vals",datatype(Float64), dataspace(N_eigvals_steps,max_vals))
    dsets[12] = create_dataset(g12,"im_vals",datatype(Float64), dataspace(N_eigvals_steps,max_vals))

end
N_shells > 0 && (gshells = create_group(file,"shells"); 
dsets[16] =create_dataset(gshells,"stresses",datatype(Float64), dataspace(N_saved_steps,12,N_shells));
dsets[17] =create_dataset(gshells,"numbers", datatype(Int), dataspace(N_shells,1));
dsets[24] =create_dataset(gshells,"strains",datatype(Float64), dataspace(N_saved_steps,12,N_shells)))
N_beams > 0 && (gbeams = create_group(file,"beams"); 
dsets[18] =create_dataset(gbeams,"stresses",datatype(Float64), dataspace(N_saved_steps,6,N_beams));
dsets[19] =create_dataset(gbeams,"numbers", datatype(Int), dataspace(N_beams,1));
dsets[25] =create_dataset(gbeams,"strains",datatype(Float64), dataspace(N_saved_steps,6,N_beams)))

dsets[30][:,1] = model_data
N_shells > 0 && (dsets[17][:,1] = sc.numbers[:])
N_beams > 0 && (dsets[19][:,1] = bc.numbers[:])

model_data[1:10] = [Nnodes,Nel,Ndofs,Ndofs_x,Ndofs_q,Ndofs_int, Ndofs_mult,Ndofs_v,Nbounds,max_vals]
dsets[30][:,1] = model_data

return file, dsets

end