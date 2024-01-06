function read_string(out_msg::String)

    strg = nothing

    while strg === nothing
        write(stdout,out_msg);
        try
            strg = readline()
        catch
        end
      end

    return strg

end
