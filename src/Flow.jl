module Flow

export AbstractFlow, set_flow!

using FFTW

abstract type AbstractFlow end

# ========================================================================

struct NoFlow <: AbstractFlow end

function Base.show(io::IO, f::NoFlow)
    print(io, "Zero flow: ux, uy = 0")
end

function set_flow!(
    f::NoFlow,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])

    ux = zeros(type, nx, ny)
    uy = zeros(type, nx, ny)

    return ux, uy
end

# ========================================================================

Base.@kwdef struct ShearFlow <: AbstractFlow
    U_max = 0.03
    α = 80.0
    δ = 0.05
    H = 2π
end

function Base.show(io::IO, f::ShearFlow)
    print(io, "Shear flow: ux = $(f.U_max) tanh($(f.α)*($(f.H)/4 - |y - 0.5*$(f.H)|)), uy = $(f.U_max)*$(f.δ)*sin(2π(x + 0.25*$(f.H))/$(f.H))")
end

function set_flow!(
    f::ShearFlow,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])
    ux = zeros(type, nx, ny)
    uy = zeros(type, nx, ny)

    U_max = type(f.U_max)
    α = type(f.α)
    δ = type(f.δ)
    H = type(f.H)

    for j in 1:ny
        for i in 1:nx
            ux[i, j] = U_max * tanh(α * (H / 4 - abs(y[j] - 0.5 * H)))
            uy[i, j] = U_max * δ * sin(2π * (x[i] + 0.25 * H) / H)
        end
    end

    return ux, uy
end


# ========================================================================

Base.@kwdef struct TaylorGreenFlow <: AbstractFlow
    C = 0.01
    κ = 1
end

function Base.show(io::IO, f::TaylorGreenFlow)
    print(io, "Taylor-Green flow ux = $(f.C)*sin(2πx * $(f.κ))*cos(2πy * $(f.κ)), uy = -$(f.C)*cos(2πx * $(f.κ))*sin(2πy * $(f.κ))")
end

function set_flow!(
    f::TaylorGreenFlow,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])
    ux = zeros(type, nx, ny)
    uy = zeros(type, nx, ny)

    C = type(f.C)
    κ = type(f.κ)

    for j in 1:ny
        for i in 1:nx
            ux[i, j] = C * sin(x[i] * κ) * cos(y[j] * κ)
            uy[i, j] = -C * cos(x[i] * κ) * sin(y[j] * κ)
        end
    end

    return ux, uy
end

# ========================================================================

Base.@kwdef struct VortexMerginFlow <: AbstractFlow
    fftfunc::Function = FFTW.fft
    ifftfunc::Function = FFTW.ifft
end

function Base.show(io::IO, f::VortexMerginFlow)
    print(io, "Vortex Merging flow : ω = exp(-π((x-3π/4)²+(y-π)²)) + exp(-π((x-5π/4)² + (y-π)^2))")
end

function set_flow!(
    f::VortexMerginFlow,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])
    type_c = Complex{type}
    ω = zeros(type, nx, ny)
    for j in 1:ny
        for i in 1:nx
            ω[i, j] = exp(-π * ((x[i] - 3π / 4)^2 + (y[j] - π)^2)) + exp(-π * ((x[i] - 5π / 4)^2 + (y[j] - π)^2))
        end
    end

    # Fourier transform
    ω̂ = f.fftfunc(ω)
    ϕ̂ = zeros(type_c, nx, ny)
    for j in 1:ny
        for i in 1:nx
            ϕ̂[i, j] = ω̂[i, j] / (kx[i]^2 + ky[j]^2)
            if abs(kx[i]) < 1.0e-6 && abs(ky[j]) < 1.0e-6
                ϕ̂[i, j] = 0.0f0 + im * 0.0f0
            end
        end
    end

    # calculate the velocity
    ux_f = zeros(type_c, nx, ny)
    uy_f = zeros(type_c, nx, ny)
    for j in 1:ny
        for i in 1:nx
            ux_f[i, j] = im * ky[j] * ϕ̂[i, j]
            uy_f[i, j] = -im * kx[i] * ϕ̂[i, j]
        end
    end

    # Inverse Fourier transform
    ux = f.ifftfunc(ux_f)
    uy = f.ifftfunc(uy_f)

    return ux, uy
end

# ========================================================================

Base.@kwdef struct PrescribedSpectrumFlow <: AbstractFlow
    s = 3
    kₚ = 12
    fftfunc::Function = FFTW.fft
    ifftfunc::Function = FFTW.ifft
end

function Base.show(io::IO, f::PrescribedSpectrumFlow)
    aₛ = (2 * f.s + 1)^(f.s + 1) / 2^f.s / factorial(f.s)
    print(io, "Prescribed Spectrum flow : E(k) = $(aₛ)/2 * 1/$(f.kₚ) * (k/$(f.kₚ))^(2*$(f.s)+1) exp[-($(f.s)+1/2)(k/$(f.kₚ))²]")
end

function set_flow!(
    f::PrescribedSpectrumFlow,
    nx::Int,
    ny::Int,
    n_total::Int,
    x::Vector,
    y::Vector,
    kx::Vector,
    ky::Vector,
)

    type = eltype(x[1])
    type_c = Complex{type}

    s = type(f.s)
    kₚ = type(f.kₚ)
    aₛ = (2 * s + 1)^(s + 1) / 2^s / factorial(s)

    E = zeros(type, nx, ny)
    for j in 1:ny
        for i in 1:nx
            k = sqrt(kx[i]^2 + ky[j]^2)
            if k > 0
                E[i, j] = aₛ / (2 * kₚ) * (k / kₚ)^(2 * s + 1) * exp(-(s + 1 / 2) * (k / kₚ)^2)
            else
                E[i, j] = 0.0f0
            end
        end
    end

    ω̂ = zeros(type_c, nx, ny)
    for j in 1:ny
        for i in 1:nx
            k = sqrt(kx[i]^2 + ky[j]^2)
            ω̂[i, j] = sqrt(E[i, j] * k / π)
        end
    end

    # Set up two random vectors ξ(kx, ky) and η(kx, ky) ∈ [0, 2π]
    ξ = rand(type, nx, ny) * 2π
    η = rand(type, nx, ny) * 2π

    # Enforce conjugate relations
    # ξ(-kx, ky) = -ξ(kx, ky)
    # ξ(-kx, -ky) = -ξ(kx, ky)
    # ξ(kx, -ky) = ξ(kx, ky)
    ξ[begin:end, ny-ny÷2+2:ny] = ξ[begin:end, ny÷2:-1:2]
    # η(-kx, ky) = η(kx, ky)
    # η(-kx, -ky) = -η(kx, ky)
    # η(kx, -ky) = -η(kx, ky)
    η[begin:end, ny-ny÷2+2:ny] = -η[begin:end, ny÷2:-1:2]

    # Set up phase function
    # exp(iζ(kx, ky)) where ζ(kx, ky) = ξ(kx, ky) + η(kx, ky)
    ζ = ξ .+ η
    ω̂ = @. ω̂ * exp(im * ζ)

    # Set up the velocity field
    ψ = zeros(ComplexF32, nx, ny)
    for j in 1:ny
        for i in 1:nx
            k² = kx[i]^2 + ky[j]^2
            if k² > 0
                ψ[i, j] = ω̂[i, j] / k²
            else
                ψ[i, j] = 0.0f0
            end
        end
    end

    # u = ∂ψ/∂y, v = -∂ψ/∂x
    û = zeros(ComplexF32, nx, ny)
    v̂ = zeros(ComplexF32, nx, ny)
    for j in 1:ny
        for i in 1:nx
            û[i, j] = im * ky[j] * ψ[i, j]
            v̂[i, j] = -im * kx[i] * ψ[i, j]
        end
    end

    ux = f.ifftfunc(û) * n_total
    uy = f.ifftfunc(v̂) * n_total

    return ux, uy

end

end
