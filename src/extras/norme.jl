norme(a::Vec3) = sqrt(sum(@.a.v^2))
norme(a::RV3) = sqrt(sum(@.a.v^2))
norme(a) = sqrt(sum(@.a^2))
