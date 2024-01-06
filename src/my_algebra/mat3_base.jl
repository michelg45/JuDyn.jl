function Base.getindex(M::Mat3,i::Int,j::Int)
    return M.mat[i,j]
end

function Base.setindex!(M::Mat3,val,i::Int,j::Int)
    M.mat[i,j] = val
end

Base.:/(a::Mat3, b::Float64) = Mat3(a.mat/b)
Base.:*(a::Mat3, b::Float64) = Mat3(a.mat*b)
Base.:*(a::Float64, b::Mat3) = Mat3(a*b.mat)
Base.:+(a::Mat3, b::Mat3) = Mat3(a.mat+b.mat)
Base.:-(a::Mat3, b::Mat3) = Mat3(a.mat-b.mat)
Base.:*(a::Mat3, b::Mat3) = Mat3(a.mat*b.mat)
Base.:*(a::Mat3, b::Vec3) = Vec3(a.mat*b.v)
Base.transpose(a::Mat3) = Mat3(transpose(a.mat))


# const I33=[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
const I3=Mat3([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
# const Z33=zeros(3,3)

Mat3(a::Vec3,b::Vec3)   = Mat3([a[1]*b[1] a[1]*b[2] a[1]*b[3];a[2]*b[1] a[2]*b[2] a[2]*b[3];a[3]*b[1] a[3]*b[2] a[3]*b[3]])
outerp(a::Vec3,b::Vec3) = Mat3([a[1]*b[1] a[1]*b[2] a[1]*b[3];a[2]*b[1] a[2]*b[2] a[2]*b[3];a[3]*b[1] a[3]*b[2] a[3]*b[3]])
tilde(psi::Vec3) = Mat3([0.0 -psi[3] psi[2]; psi[3] 0.0 -psi[1]; -psi[2] psi[1] 0.0])
    
