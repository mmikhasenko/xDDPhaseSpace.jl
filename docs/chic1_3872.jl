using xDDPhaseSpace
using Plots
theme(:boxed)

# Define particle masses (GeV)
const mπ⁰ = 0.1349768
const mD⁰ = 1.86483
const mDˣ⁰ = 2.00685
const ΓDˣ⁰ = 55.2e-6
const thr = mD⁰ + mDˣ⁰
# 
const mγ = 0.0
# to match 2b and 3b phase space
const e_matching = WithThrE(30, thr) # 500 times the width away

# πDD: Pion decay to D meson pairs
decay_channel_pi =
    πDD((m1 = mπ⁰, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰))

# γDD: Photon decay to D meson pairs with electromagnetic couplings
μ₀ = -3.77  # Transition matrix element < D0 | mu | D*0 >
decay_channel_gamma = γDD(
    (m1 = mγ, m2 = mD⁰, m3 = mD⁰),  # m1 = 0.0 for photon
    BW(m = mDˣ⁰, Γ = ΓDˣ⁰),
    BW(m = mDˣ⁰, Γ = ΓDˣ⁰),
    μ₀,
    -μ₀,
)

channels = (decay_channel_pi, decay_channel_gamma)
labs = ("π", "γ")
asymptotics = [ρ_thr(ch, e_matching) for ch in channels]

# 
let
    plot(xlab = "Δm(Dˣ⁰D⁰) [GeV]", ylab = "ρ(m)")
    for (ch, asy, lab) in zip(channels, asymptotics, labs)
        plot!(e -> ρ_thr(ch, WithThrE(e, thr)) / asy, -2, 2; lab)
    end
    plot!(yscale = :log10)
end
