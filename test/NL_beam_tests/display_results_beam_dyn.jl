function display_results_beam_dyn()

solx = y[:,idx]
L = solx[1,2]
solx[:,1] = solx[:,1] .- 0.5*L
solx[:,2] = solx[:,2] .- L
soly = y[:,idy]
solz = y[:,idz]
println("max. x displacement at mid-length",  findmin(solx[:,1]))
println("max. y displacement at mid-length",  findmax(soly[:,1]))
println("max. z displacement at mid-length",  findmax(solz[:,1]))
println("max. x displacement at tip",  findmin(solx[:,2]))
println("max. y displacement at tip",  findmax(soly[:,2]))
println("max. z displacement at tip",  findmax(solz[:,2]))


titre = " lateral displacements "
lab = ["mid-length"  "tip"]
ylab = L"u_y"*" (m)"
xlab = "t (s)"
display(plot(times,soly,title = titre, label = lab, xlabel = xlab, ylabel = ylab))

titre = " vertical displacements "
lab = ["mid-length"  "tip"]
xlab = "t (s)"
ylab = L"u_z"*" (m)"
display(plot(times,solz,title = titre, label = lab, xlabel = xlab, ylabel = ylab))

titre = " axial displacements "
lab = ["mid-length"  "tip"]
ylab = L"u_x"*" (m)"
display(plot(times,solx,title = titre, label = lab,xlabel = xlab, ylabel = ylab))

titre = " mid-length displacements "
ylab = "(m)"
plot(times,solx[:,1],label = "x",xlabel = xlab, ylabel = ylab)
plot!(times,soly[:,1],label = "y")
display(plot!(times,solz[:,1],title = titre,  label = "z"))

titre = " tip displacements "
plot(times,solx[:,2],label = "x",xlabel = xlab, ylabel = ylab)
plot!(times,soly[:,2],label = "y")
display(plot!(times,solz[:,2],title = titre, label = "z"))

titre = "kinetic energy "
ylab = "KE (Nm)"
display(plot(times,kin_energy, title = titre, xlabel = xlab, ylabel = ylab,legend=:none))

titre = "strain  energy "
ylab = "SE (Nm)"
display(plot(times,str_energy, title = titre,xlabel = xlab, ylabel = ylab,legend=:none))

titre = "external_work "
ylab = "W (Nm)"
display(plot(times,ext_work, title = titre,xlabel = xlab, ylabel = ylab,legend=:none))

energy_norm = findmax(kin_energy+str_energy)[1]
relative_energy_balance =   (ext_work - kin_energy - str_energy)/energy_norm

titre =  "relative_energy balance "
ylab = " "
display(plot(times,relative_energy_balance, title = titre,xlabel = xlab, ylabel = ylab,legend=:none))

titre = "iterations per step "
ylab = "N"
display(plot(times,nitmax, title = titre, xlabel = xlab, ylabel = ylab,legend=:none))



end
