function read_shell_stresses(h5_file)


    file=h5open(h5_file,"r")
    if haskey(file,"shells") == true 
        stresses = read(file["shells"],"stresses")
        numbers = read(file["shells"],"numbers")
        close(file)
        return numbers, stresses
    else
        println("no shell stresses on ", h5_file)
        close(file)
        return
    end
end
  