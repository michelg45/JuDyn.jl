function read_number(out_msg::String)

  num = nothing

  while num === nothing
    write(stdout,out_msg);
    try
      num = tryparse(Int64,readline())
    catch
    end
  end

  return num
end
