function linearized_rigid_body(k::Vec3, inertia_parameters::Vector{Float64})

   abs(norm2(k)-1.0) > eps(Float64)*10.0  && (error(" k must be a unit vector"))
   size(inertia_parameters,1) != 4 && (error("wrong number of inertia parameters"))

   m = inertia_parameters[1]

   J = inertia_parameters[2:4]

   M = diagm([m; m; m; J])

   unit_v = Vec3(ones(3))

   u = k - dotp(k,unit_v)*unit_v

   tk = tilde(k)

   G = zeros(6,6)

   G[1:3,1:3] = 2.0*m*tk.mat

   JJ = Mat3(M[4:6,4:6])

   G[4:6,4:6] = (JJ*tk + tk*JJ -tilde(JJ*k)).mat

   Kgyr =  zeros(6,6)

   Kgyr[1:3,1:3] = m*diagm(u.v)

   Kgyr[4:6,4:6] = (tk*JJ*tk  -tilde(JJ*k)*tk).mat

   return M, G, Kgyr
   
end
