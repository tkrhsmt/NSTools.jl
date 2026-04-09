module Spectral

using FFTW

export make_spectral
export energy_spectral, prod_spectral
export helmholz_decomp

function make_spectral(
    E_f::Array,
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

    # wavenumber grid
    grid = size(ux)
    n_total = grid[1] * grid[2]

    # fourier transform
    ux_f = FFTW.fft(ux) / n_total
    uy_f = FFTW.fft(uy) / n_total

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
    n_total = grid[1] * grid[2]

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
    T1_f = FFTW.fft(T1) / n_total
    T2_f = FFTW.fft(T2) / n_total
    # division for normalization
    ux_f = ux_f / n_total
    uy_f = uy_f / n_total

    # nonlinear term in spectral space
    T_f = @. T1_f * conj(ux_f) / 2 + T2_f * conj(uy_f) / 2

    # make spectral object
    return make_spectral(T_f)
end

function helmholz_decomp(
    ux::Array,
    uy::Array
)

    # wavenumber grid
    grid = size(ux)
    kx = fftfreq(grid[1], grid[1])
    ky = fftfreq(grid[2], grid[2])

    ux_f = FFTW.fft(ux)
    uy_f = FFTW.fft(uy)

    # Helmholtz decomposition in spectral space
    for j in 1:grid[2]
        for i in 1:grid[1]
            if kx[i] == 0 && ky[j] == 0
                ux_f[i, j] = 0.0 + 0.0im
                uy_f[i, j] = 0.0 + 0.0im
            else
                tmp1 = kx[i] * (kx[i] * ux_f[i, j] + ky[j] * uy_f[i, j]) / (kx[i]^2 + ky[j]^2)
                tmp2 = ky[j] * (kx[i] * ux_f[i, j] + ky[j] * uy_f[i, j]) / (kx[i]^2 + ky[j]^2)
                ux_f[i, j] = tmp1
                uy_f[i, j] = tmp2
            end
        end
    end

    ux_d = real(FFTW.ifft(ux_f))
    uy_d = real(FFTW.ifft(uy_f))

    ux_s = ux .- ux_d
    uy_s = uy .- uy_d

    return (ux_s, uy_s), (ux_d, uy_d)

end

end
