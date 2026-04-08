"""
    generate_ticks(x)

Generate ticks covering x spanned by π/2.
"""
function generate_ticks(x)
    step = π / 2
    ticks = x[1]:step:x[end]
    num_ticks = round(Int, x[end] / step)

    labels = []
    for i in 0:num_ticks
        # reduction
        gcd_value = gcd(i, 2)  # "2" comes from the denominator of the step π/2
        num_simplified = i ÷ gcd_value
        # println("i=", i, " num_simplified=", num_simplified, " gcd_value=", gcd_value)

        # Create label
        if num_simplified == 0
            push!(labels, L"0")  # 0π -> 0
        elseif num_simplified == 1 && gcd_value == 1
            push!(labels, L"\pi/2")  # 1π/2 -> π/2
        elseif num_simplified == 1 && gcd_value == 2
            push!(labels, L"\pi")  # 1π -> π
        elseif gcd_value == 2
            push!(labels, L"%$num_simplified\pi")  # nπ
        elseif gcd_value == 1
            push!(labels, L"%$num_simplified\pi/2")  # nπ/2
        end
    end

    return ticks, labels
end

function NSTools.plot_data(
    x::Vector,
    y::Vector,
    ω::Array,
    filename::String
)

    cmap = :balance
    vmin, vmax = minimum(ω), maximum(ω)
    climits = (-max(abs(vmin), abs(vmax)), max(abs(vmin), abs(vmax)))

    if vmin > 0.0f0 || vmax < 0.0f0 && (vmin != vmax)
        # If the minimum and maximum values are both positive or both negative, use a different colormap
        cmap = :plasma
        climits = (vmin, vmax)  # to enforce white (middle of the colourmap) with pp=0
    end

    xticks, xlabel = generate_ticks(x)
    yticks, ylabel = generate_ticks(y)

    with_theme(colormap=cmap) do
        fig = Figure(fontsize=22)
        ax = Axis(
            fig[1, 1];
            aspect=1,
            xlabel=L"$x$",
            ylabel=L"$y$",
            xticks=(xticks, xlabel),
            yticks=(yticks, ylabel)
        )
        pltobj = heatmap!(ax, x, y, ω; colormap=(cmap, 0.7), colorrange=climits)
        Colorbar(fig[1, 2], pltobj)
        colsize!(fig.layout, 1, Aspect(1, 1.0))
        resize_to_layout!(fig)
        save(filename, fig)

        println("    Saved: ", filename)

    end

end
