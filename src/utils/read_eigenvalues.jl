using HDF5


function read_eigenvalues(h5_file)

    """
        "read_eigenvalues" is a function reading on the HDF5 file
         the results of the eigenvalue analysis.

    """

file=h5open(h5_file,"r")
gtimes = file["time_vals"]
grvals=file["re_vals"]
gimvals=file["im_vals"]

times = read(gtimes, "time_vals")
re_vals = read(grvals, "re_vals")
im_vals = read(gimvals, "im_vals")

close(file)

return times, re_vals, im_vals

end