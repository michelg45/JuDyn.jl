"""
    read_beam_properties
"""
function read_beam_properties(JSON_file::String,name::String)

  open(JSON_file,"r") do f
         global dicts
         dicts = read(f,String)
  end

  dict =  JSON.parse(dicts)

  length = dict[name]["General"]["length"]
  shear = dict[name]["General"]["shear"]
  rotary_inertia = dict[name]["General"]["rotary_inertia"]
  visco = haskey(dict[name],"Visco")
 
  if shear == true
    stiffness_properties = Vector{Float64}(undef,6)
    stiffness_properties[1] = dict[name]["Stiffness"]["EA"]
    stiffness_properties[2] = dict[name]["Stiffness"]["kAG_y"]
    stiffness_properties[3] = dict[name]["Stiffness"]["kAG_z"]
    stiffness_properties[4] = dict[name]["Stiffness"]["GJ"]
    stiffness_properties[5] = dict[name]["Stiffness"]["EI_yy"]
    stiffness_properties[6] = dict[name]["Stiffness"]["EI_zz"]
  else
    stiffness_properties = Vector{Float64}(undef,4)
    stiffness_properties[1] = dict[name]["Stiffness"]["EA"]
    stiffness_properties[2] = dict[name]["Stiffness"]["GJ"]
    stiffness_properties[3] = dict[name]["Stiffness"]["EI_yy"]
    stiffness_properties[4] = dict[name]["Stiffness"]["EI_zz"]
  end

  if shear == true
    mass_properties = Vector{Float64}(undef,4)
    mass_properties[1] = dict[name]["Mass"]["m"]
    mass_properties[2] = dict[name]["Mass"]["J_xx"]
    mass_properties[3] = dict[name]["Mass"]["J_yy"]
    mass_properties[4] = dict[name]["Mass"]["J_zz"]
  else
    mass_properties = Vector{Float64}(undef,4)
    mass_properties[1] = dict[name]["Mass"]["m"]
    mass_properties[2] = dict[name]["Mass"]["J_zz"]
  end

  if visco == true
    tau_B = dict[name]["Visco"]["tau_B"]
    tau_S = dict[name]["Visco"]["tau_S"]
    nu    = dict[name]["Visco"]["nu"]
    ratio_infty = dict[name]["Visco"]["ratio_infty"]
    return length, stiffness_properties, mass_properties,tau_B,tau_S,nu, ratio_infty
  else
    return length, stiffness_properties, mass_properties
  end

end
