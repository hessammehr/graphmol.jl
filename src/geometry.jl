#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    rotation,
    interiorangle,
    utov,
    translate,
    trim_u,
    trim_v,
    trim_uv,
    isclockwise,
    radiantophase,
    transformmatrix


import LinearAlgebra: cross
using LinearAlgebra
using StaticArrays


rotation(a::Real) = [cos(a) -sin(a); sin(a) cos(a)]

cross(u::SVector{2}, v::SVector{2}) = u[1] * v[2] - u[2] * v[1]
interiorangle(u::SVector{2}, v::SVector{2}) = acos(
    dot(u, v) / (hypot(u...) * hypot(v...)))

utov(uv::SMatrix{2,2}) = uv[2, :] - uv[1, :]


function translate(uv::SMatrix{2,2}, angle::Real, dist::Real)
    rot = rotation(angle)
    move = rot * utov(uv) / hypot(utov(uv)...) * dist
    uv .+ transpose(move)
end

trim_u(uv::SMatrix{2,2}, k::Real) = uv + [transpose(utov(uv) * k); 0 0]
trim_v(uv::SMatrix{2,2}, k::Real) = uv + [0 0; transpose(utov(uv) * -k)]
trim_uv(uv::SMatrix{2,2}, k::Real) = uv + [
    transpose(utov(uv) * k / 2); transpose(utov(uv) * -k / 2)]


function isclockwise(vertices::AbstractArray{SVector{2}})
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycver = cat(vertices, vertices, dims=1)
    for i = 1:vlen
        (p, q, r) = cycver[i:i+2]
        cp = cross(q - p, r - p)
        intangle = interiorangle(p - q, r - q)
        if cp < 0
            clockwise += intangle
            counter += 2pi - intangle
        else
            clockwise += 2pi - intangle
            counter += intangle
        end
    end
    if round(clockwise / pi) == vlen - 2
        true
    elseif round(counter / pi) == vlen - 2
        false
    else
        NaN
    end
end


function radiantophase(angle::Real)
    mod((angle + 2pi) / 2pi)
end


function transformmatrix(scale::SVector{2}, rotvec::SVector{2},
                         translate::SVector{2})
    scale = SVector{3,3}[scale[1] 0 0; 0 scale[2] 0; 0 0 1]
    rot = [rotvec[1] -rotvec[2] 0; rotvec[2] rotvec[1] 0; 0 0 1]
    tl = [1 0 translate[1]; 0 1 translate[2]; 0 0 1]
    tl * rot * scale
end


function rotation(axis::SVector{3}, angle::Real)
    (x, y, z) = axis
    c = cos(angle)
    s = sin(angle)
    a12 = x * y * (1 - c)
    a13 = x * z * (1 - c)
    a23 = y * z * (1 - c)
    SVector{3,3}[
        (c + x^2 * (1 - c)) (a12 - z * s) (a13 + y * s);
        (a12 + z * s) (c + y^2 * (1 - c)) (a23 - x * s);
        (a13 - y * s) (a23 + x * s) (c + z^2 * (1 - c))
    ]
end
