function read_strain_stresses(h5_file,el_type)

    file=h5open(h5_file,"r")
    if haskey(file,el_type) == true 
        stresses = read(file[el_type],"stresses")
        strains = read(file[el_type],"strains")
        numbers = read(file[el_type],"numbers")
        close(file)
        return numbers, strains, stresses
    else
        println("no strains / stresses for element type ", el_type,  " on ", h5_file)
        close(file)
        return
    end
end