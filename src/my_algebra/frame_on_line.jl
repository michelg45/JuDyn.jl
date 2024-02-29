"""
    frame_on_line

Function construction a rotation frame with its ``\\mathbf{n}_1`` axis aligned 
on the line ``(\\mathbf{x}_1, \\mathbf{x}_2)``.

| Input: | |
|:---------------|:------------------------|
| x1, x2::Vec3 | node positions. |
| | |
| **Output:**  | | 
| R::RV3 | frame rotation vector | 

Calling sequence:

    R = frame_on_line(x1,x2)

"""
function frame_on_line(x_1::Vec3,x_2::Vec3)

    dx = x_2 - x_1

    n_1 = dx/norm2(dx)
     
    cs = 1.0

    global n_2 = Vec3() 

    for i = 1:5
        rd = Vec3(rand(3))
        n_2 = copy(rd/norm2(rd))
        cs = dotp(n_1,n_2)

        abs(cs) <= 0.5 && break

    end

    n_3 = crossp(n_1,n_2)
    n_3 = n_3/norm2(n_3)
    n_2 = crossp(n_3,n_1)
    n_2 = n_2/norm2(n_2)

    return R = invrot(Mat3(n_1,n_2,n_3))

end