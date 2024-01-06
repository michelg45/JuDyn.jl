"""
   save_results

    function to save the temporary h5_file as a permanent one named "name.h5"

      calling sequence: save_results(h5_file,name)
"""
function save_results(h5_file::String, label::String)

  h5_new_name = label*".h5"

  if h5_new_name != h5_file
    mv(h5_file, h5_new_name, force = true)
  end

  println("results saved on ", h5_new_name)

  return h5_new_name

end
