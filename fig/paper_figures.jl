# paper_figures.jl
# Create all figures for the paper and save

output_folder = "paper-figures"

include("plotting.jl")
include("noise.jl")
include("tke.jl")
include("w_slices.jl")
include("wVSP.jl")
include("hovmoller.jl")
include("psi_balance.jl")
include("snapshots.jl")

foldername = "../scratch/filament-instability-redux/Ri/Ri00"
preinit_foldername = "../scratch/filament-instability-redux/preinit/m2"

Ri_titles = [L"Ri_\text{min}=0.0", L"Ri_\text{min}=0.1", L"Ri_\text{min}=0.2"]
Ri_foldernames = map(a->joinpath("../scratch/filament-instability-redux/Ri", a), ["Ri00", "Ri01", "Ri02"]);

amplitude_titles = [L"10^{-4}", L"10^{-6}", L"10^{-8}"]
amplitude_foldernames = map(a->joinpath("../scratch/filament-instability-redux/amplitude", a), ["m4", "m6", "m8"])

# Main text
filenames = [
    "noise.png",
    "tke.png",
    "w1.png",
    "w2.png",
    "snapshots.png",
    "hovmoller.png",
    "psi.png"
]

figures = [
    noise_figure(preinit_foldername),
    tke_figure(Ri_foldernames, Ri_titles),
    wVSP_figure(foldername, 1.8, -0.05; σ=0, σh=1),
    wVSP_figure(foldername, 3.6, -0.05; σ=0, σh=1),
    snapshots_figure(foldername, 2π .* [0.7, 1.4, 2.1]),
    hovmoller_figure(foldername; marked_times=[0.7, 1.4, 2.1]),
    psi_balance_figure(Ri_foldernames, Ri_titles)
]

map(filenames, figures) do filename, fig
    println(filename)
    save(joinpath(output_folder, filename), fig; px_per_unit=2)
end

# Appendix
filenames = [
    "amplitude-w.png",
    "amplitude-tke.png",
    "amplitude-hovmoller.png",
    "amplitude-psi.png"
]

figures = [
    w_slices_figure(amplitude_foldernames[end], 4, 4.5, -0.05),
    tke_figure(amplitude_foldernames, amplitude_titles; ax_kw=(; limits=(0, 4, -0.2, 1.7)), tke₀=true),
    hovmoller_figure(amplitude_foldernames[end]),
    psi_balance_figure(amplitude_foldernames, amplitude_titles; σ=3)
]

map(filenames, figures) do filename, fig
    println(filename)
    save(joinpath(output_folder, filename), fig; px_per_unit=2)
end