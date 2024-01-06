function read_positive_number(out_msg::String)

  global num = nothing

  while num === nothing || num <= 0
    write(stdout,out_msg);
    try
     num = tryparse(Int64,readline())
    catch 
    end
  end
  return num
end
