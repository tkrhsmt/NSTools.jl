module Spectral

using FFTW

export make_spectral
export energy_spectral, prod_spectral

function make_spectral(
    E_f::Array{Complex{AbstractFloat},2},
)

    # wavenumber grid
    grid = size(E_f)
    kx = fftfreq(grid[1], grid[1])
    ky = fftfreq(grid[2], grid[2])
    kmax = Int(ceil(maximum(vcat(kx, ky))))
    k = [range(1.0, kmax);]

    # energy spectrum
    E = zeros(kmax)
    for j in 1:grid[2]
        for i in 1:grid[1]
            kx_local = kx[i]
            ky_local = ky[j]
            kk = Int(ceil(sqrt(kx_local^2 + ky_local^2)))

            if kk == 0
                kk = 1
            end
            if kmax >= kk
                E[kk] += real(E_f[i, j])
            end
        end
    end

    return k, E
end

function energy_spectral(
    ux::Array,
    uy::Array
)

    # fourier transform
    ux_f = FFTW.fft(ux)
    uy_f = FFTW.fft(uy)

    # energy in spectral space
    E_f = @. ux_f * conj(ux_f) / 2 + uy_f * conj(uy_f) / 2

    # make spectral object
    return make_spectral(E_f)
end

function prod_spectral(
    ux::Array,
    uy::Array
)

    # wavenumber grid
    grid = size(ux)
    kx = fftfreq(grid[1], grid[1])
    ky = fftfreq(grid[2], grid[2])

    # fourier transform
    ux_f = FFTW.fft(ux)
    uy_f = FFTW.fft(uy)

    # compute nonlinear term in spectral space
    uxdx_f = im .* ux_f
    uxdy_f = im .* ux_f
    uydx_f = im .* uy_f
    uydy_f = im .* uy_f

    for j in 1:grid[2]
        for i in 1:grid[1]
            uxdx_f[i, j] *= kx[i]
            uxdy_f[i, j] *= ky[j]
            uydx_f[i, j] *= kx[i]
            uydy_f[i, j] *= ky[j]
        end
    end

    T1 = ux .* FFTW.ifft(uxdx_f) + uy .* FFTW.ifft(uxdy_f)
    T2 = ux .* FFTW.ifft(uydx_f) + uy .* FFTW.ifft(uydy_f)

    # Fourier transform of the nonlinear term
    T1_f = FFTW.fft(T1)
    T2_f = FFTW.fft(T2)

    # nonlinear term in spectral space
    T_f = @. T1_f * conj(ux_f) / 2 + T2_f * conj(uy_f) / 2

    # make spectral object
    return make_spectral(T_f)
end

end
