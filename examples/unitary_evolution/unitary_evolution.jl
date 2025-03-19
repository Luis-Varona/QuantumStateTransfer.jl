using QuantumStateTransfer: UnitaryEvolution, unitary_evolution
using Graphs: cycle_graph

using FileIO
using ImageIO, ImageMagick, ImageShow, Images, OpenCV
using Plots, Plots.PlotMeasures


"""
    plot_c4_amplitudes(c4_evolution::UnitaryEvolution, source::Int)

[To do]

# Arguments
- `c4_evolution::UnitaryEvolution`: [To do]
- `source::Int`: [To do]

# Returns
- `Matrix{RGB{N0f8}}`: [To do]

# Examples
[To do]
"""
function plot_c4_amps(c4_evolution::UnitaryEvolution, source::Int)
    time_steps = c4_evolution.time_steps
    transfer_amplitudes = c4_evolution.transfer_amplitudes[:, :, source]
    labels = permutedims(map(i -> replace("to node $i", "node $source" => "itself"), 1:4))
    lc = permutedims(circshift([:pink, :crimson, :green, :dodgerblue], source))
    lw = permutedims(circshift([8, 2, 2, 2], source))
    la = permutedims(circshift([0.5, 1, 2, 1], source))
    
    amplitude_plot = plot(
        time_steps,
        real.(transfer_amplitudes),
        imag.(transfer_amplitudes),
        xlabel="Time",
        ylabel="Real",
        zlabel="Imaginary",
        labels=labels,
        title="Transfer Amplitude from Node $source on the Cycle Graph C₄",
        lc=lc,
        lw=lw,
        la=la,
        legend=(0.8, 0.8),
        titlelocation=(0.48, 0.96),
        top_margin=-35px,
        bottom_margin=-55px,
        left_margin=-35px,
        size=(640, 640),
        dpi=300,
    )
    
    io = IOBuffer()
    png(amplitude_plot, io)
    return load(io)[26:1825, 61:1860]
end

function plot_c4_fids(c4_evolution::UnitaryEvolution, source::Int)
    time_steps = c4_evolution.time_steps
    transfer_fidelities = c4_evolution.transfer_fidelities[:, :, source]
    labels = permutedims(map(v -> replace("to node $v", "node $source" => "itself"), 1:4))
    lc = permutedims(circshift([:pink, :crimson, :green, :dodgerblue], source))
    lw = permutedims(circshift([8, 2, 2, 2], source))
    la = permutedims(circshift([0.5, 1, 2, 1], source))
    
    fidelity_plot = plot(
        time_steps,
        transfer_fidelities,
        xlabel="Time",
        ylabel="Fidelity",
        labels=labels,
        title="Transfer Fidelity from Node $source on the Cycle Graph C₄",
        lc=lc,
        lw=lw,
        la=la,
        legend=(0.88, 1.1),
        titlelocation=(0.38, 1.1),
        top_margin=35px,
        left_margin=10px,
        size=(760, 570),
        dpi=300,
    )
    
    io = IOBuffer()
    png(fidelity_plot, io)
    return load(io)
end


c4 = cycle_graph(4)
c4_evolution = unitary_evolution(c4)
c4_transfer_amps_plots = map(u -> plot_c4_amps(c4_evolution, u), 1:4)
c4_transfer_fids_plots = map(u -> plot_c4_fids(c4_evolution, u), 1:4)

dest_amps = joinpath(@__DIR__, "cycle_graph_4_plots", "transfer_amps")
dest_fids = joinpath(@__DIR__, "cycle_graph_4_plots", "transfer_fids")
mkpath(dest_amps)
mkpath(dest_fids)

for u in 1:4
    filename_amps = "c4_transfer_amps_from_node_$u.png"
    filename_fids = "c4_transfer_fids_from_node_$u.png"
    save(joinpath(dest_amps, filename_amps), c4_transfer_amps_plots[u])
    save(joinpath(dest_fids, filename_fids), c4_transfer_fids_plots[u])
end