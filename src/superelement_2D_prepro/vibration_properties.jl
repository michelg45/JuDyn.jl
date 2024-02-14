function vibration_properties(JLD_file)


    Nodes, Ndofs, node_coordinates, node_components, ip_nodes,
    loc_vector, K, M = read_initial_model(JLD_file)

val, modes = eigen(K, M)

println("val ", val)

PREC = norm(K)/norm(M)*(sqrt(eps(Float64)))

U = zeros(Ndofs,6)
u_index = nonnuls((val .< PREC).*[i for i in 1:Ndofs])
phi_index = nonnuls((val .> PREC).*[i for i in 1:Ndofs])

println(u_index)
println(phi_index)

Nrig = size(u_index,1)
Nelast = size(phi_index,1)

println("Number of rigid body modes: ", Nrig)

if  Nrig != 6
    error("model has not the proper number of RB modes: Nrig = ", Nrig)
end

if (node_components[1] != 6)
    error("the current procedure requires 6 dofs on node 1")
end

A = [eye(3);zeros(3,3)]

U[:,1:Nrig] = modes[:,u_index]*inv(modes[1:Nrig,u_index])


Mrig = zeros(6,6)

Mrig[1:3,1:3] = transpose(U[:,1:3])*M*U[:,1:3]

U[:,4:6] -= U[:,1:3]*(Mrig[1:3,1:3]\(transpose(U[:,1:3])*M*U[:,4:6]))

Mrig[1:6,1:6] = transpose(U[:,1:6])*M*U[:,1:6]


X = zeros(Ndofs,3)

for i = 1:Nodes
    ip = ip_nodes[i]
    for j = 1:3
        X[ip+j-1,j] = node_coordinates[j,i]
    end
end

x_CM =  zeros(3)
for i = 1:3
    x_CM[i] = transpose(X[:,i])*M*U[:,i]/Mrig[i,i]
end


node_coordinates[1:3,:]  = node_coordinates[1:3,:] .- x_CM[1:3]

omega2 = val[phi_index]
Phi = modes[:,phi_index]


return Nrig, Nelast,node_coordinates, Mrig, U, omega2, Phi
end
