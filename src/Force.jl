module Force

export AbstractForce, set_force!

abstract type AbstractForce end

# ========================================================================

struct NoForce <: AbstractForce end

function set_force!(
    f::NoForce,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])

    # Set the force to zero
    fx = zeros(type, nx, ny)
    fy = zeros(type, nx, ny)

    return fx, fy
end

# ========================================================================

Base.@kwdef struct TaylorGreenForce <: AbstractForce
    C = 0.01
    κ = 1
end

function set_force!(
    f::TaylorGreenForce,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])
    C = type(f.C)
    κ = type(f.κ)

    # Set the force to zero
    fx = zeros(type, nx, ny)
    fy = zeros(type, nx, ny)
    for j in 1:ny
        for i in 1:nx
            fx[i, j] = C * sin(κ * x[i]) * cos(κ * y[j])
            fy[i, j] = -C * cos(κ * x[i]) * sin(κ * y[j])
        end
    end

    return fx, fy
end

end
