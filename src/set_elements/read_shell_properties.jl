function read_shell_properties(JSON_file::String,name::String)

  open(JSON_file,"r") do f
         global dicts
         dicts = read(f,String)
  end

  dict =  JSON.parse(dicts)

  thickness = dict[name]["General"]["thickness"]
 
    stiffness_properties = Vector{Float64}(undef,2)
    stiffness_properties[1] = dict[name]["Stiffness"]["E"]
    stiffness_properties[2] = dict[name]["Stiffness"]["nu"]

    mass_properties = Vector{Float64}(undef,1)
    mass_properties[1] = dict[name]["Mass"]["rho"]

  return thickness, stiffness_properties, mass_properties

end
