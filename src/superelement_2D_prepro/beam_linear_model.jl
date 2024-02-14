cd("src/super_element_prepro")
include("LinearStructureModel.jl")
include("../MyAlgebra/src/MyAlgebra.jl")
include("../Utils/src/Utils.jl")
include("../../src/LinearBeam.jl")
include("../../src/read_beam_properties.jl")
JSON_file = "../../src/beam_properties.json"
using LinearAlgebra
using VMLS
using JLD
using JSON
using ..MyAlgebra
using ..LinearBeam

beam_label = "Haug"
JLD_file = beam_label*"_linear_beam_model.jld"
length, stiffness_properties, mass_properties = read_beam_properties(JSON_file,beam_label)

Nel = 2
n_SE = 2
length = length/n_SE
Nodes = Nel+1
Ndofs = Nodes*6
model = LinearStructureModel(Nodes,Ndofs)

model.node_coordinates[1,:] = [(i-1)*length/Nel for i in 1:Nodes]
model.node_components .= 6
model.loc_vector = [(i-1)*6+j for i in 1:Nodes for j in 1:6]
model.ip_nodes = [(i-1)*6+1 for i in 1:Nodes]
L = length / Nel
K_el = beam_stiffness_3D_local(L,stiffness_properties)
M_el = beam_mass_3D_local(L,mass_properties)
S_el = beam_gyr_3D_local(l,m)
for iel = 1:Nel
    locel = [[(iel-1)*6+i for i in 1:6]; [iel*6+i for i in 1:6]]
    model.stiffness_matrix[locel,locel] += K_el
    model.mass_matrix[locel,locel] += M_el
    model.S_matrix[:][locel,locel] += S_el[:]  
end

file = jldopen(JLD_file,"w")
write(file,"initial_model",model)
close(file)

cd("../..")
