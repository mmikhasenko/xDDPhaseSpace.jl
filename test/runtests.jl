using xDDPhaseSpace
using Test
using Parameters
using QuadGK

@testset "scalar production" begin
    s, s12, s13, msq = 1, 2, 3, (2, 3, 5)
    v = (; s, s12, s13, msq)
    @test xDDPhaseSpace.p12_p3(v) + v.s12 == xDDPhaseSpace.p_p12(v)
    @test xDDPhaseSpace.p13_p2(v) + v.s13 == xDDPhaseSpace.p_p13(v)
    @test xDDPhaseSpace.p12_p2(v) + xDDPhaseSpace.p2_p3(v) == xDDPhaseSpace.p_p2(v)
    @test xDDPhaseSpace.p13_p2(v) + msq[2] == xDDPhaseSpace.p_p2(v)
end

@testset "ABCDEG" begin
    s, s12, s13, msq = 1, 2, 3, (2, 3, 3)
    v = (; s, s12, s13, msq)
    @test xDDPhaseSpace.A(v) != 0
    @test xDDPhaseSpace.B(v) != 0
    @test xDDPhaseSpace.C(v) != 0
    @test xDDPhaseSpace.C((; s, s12, s13, msq)) ≈
          xDDPhaseSpace.C((; s12 = s13, s13 = s12, s, msq))
    @test xDDPhaseSpace.D(v) != 0
    @test xDDPhaseSpace.E(v) != 0
    @test xDDPhaseSpace.G(v) != 0
end


const mπ⁰ = 0.1349768
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mD⁺ = 1.86965
const mDˣ⁺ = 2.0102558
const ΓDˣ⁺ = 83.4e-6

# mapping and integrating Dalitz

struct TestCh{T} <: AbstractxDD
    ms::T
end
import xDDPhaseSpace: decay_matrix_element_squared
decay_matrix_element_squared(d::TestCh, s, σ3, σ2) = 1.0

@testset "implementation of the phase-space integral" begin

    e_test = WithThrE(0.0, mD⁰ + mDˣ⁺)
    ms_test = (m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰)

    t = TestCh(ms_test)
    ρ_2dim = ρ_thr(t, e2m(e_test))

    s = e2m(e_test)^2
    ms = ms_test
    ρ_1dim =
        quadgk(
            σ1 ->
                sqrt(
                    xDDPhaseSpace.λ(s, ms.m1^2, σ1) * xDDPhaseSpace.λ(σ1, ms.m2^2, ms.m3^2),
                ) / σ1,
            (ms.m2 + ms.m3)^2,
            (sqrt(s) - ms.m1)^2,
        )[1] / (8π)^2 / (2π * s)

    @test abs(ρ_2dim - ρ_1dim) / ρ_1dim < 1e-4
end



import xDDPhaseSpace: σ3of1_pm, σ3of1, σ2of3_pm

# complex integal
@testset "σ3of1_pm and σ3of1" begin
    m = e2m(WithThrE(1.1, mD⁰ + mDˣ⁺))
    # 
    ms = (m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰)
    σ1 = (ms.m2 + ms.m3)^2 + rand() * ((m - ms.m1)^2 - (ms.m2 + ms.m3)^2)
    #
    zm1, z1 = σ3of1_pm(σ1, ms^2, m^2)
    @test zm1 ≈ σ3of1(σ1, -1, ms^2, m^2)
    @test z1 ≈ σ3of1(σ1, +1, ms^2, m^2)
    #
    σ1_x(x) = (ms.m2 + ms.m3)^2 + x * ((m - ms.m1)^2 - (ms.m2 + ms.m3)^2)
    #
    x = rand()
    @show x
    @test prod((ms.m1 + ms.m2)^2 .≤ σ3of1_pm(σ1_x(x), ms^2, m^2) .≤ (m - ms.m3)^2)
    #
    σ3_x(x) = (ms.m1 + ms.m2)^2 + x * ((m - ms.m3)^2 - (ms.m1 + ms.m2)^2)
    #
    x = rand()
    @show x
    @test prod((ms.m3 + ms.m1)^2 .≤ σ2of3_pm(σ3_x(x), ms^2, m^2) .≤ (m - ms.m2)^2)
end



decay_channel =
    πDD((m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁺, Γ = ΓDˣ⁺), BW(m = mDˣ⁺, Γ = ΓDˣ⁺))

ρ_thr(decay_channel, WithThrE(1.0, mD⁰ + mDˣ⁺) |> e2m) # 1MeV from threshold

@testset "πDD decay channel with Breit-Wigner resonances" begin
    # Create decay channel with πDD and two D*+ resonances
    decay_channel =
        πDD((m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁺, Γ = ΓDˣ⁺), BW(m = mDˣ⁺, Γ = ΓDˣ⁺))

    # Test threshold calculation 1MeV above threshold
    ρ_value = ρ_thr(decay_channel, WithThrE(1.0, mD⁰ + mDˣ⁺) |> e2m)

    # Verify the result is a positive number
    @test ρ_value > 0
    @test typeof(ρ_value) <: Real

    # Test that the decay channel has the expected structure
    @test decay_channel.ms.m1 == mπ⁺
    @test decay_channel.ms.m2 == mD⁰
    @test decay_channel.ms.m3 == mD⁰
    @test decay_channel.R12.m == mDˣ⁺
    @test decay_channel.R12.Γ == ΓDˣ⁺
    @test decay_channel.R13.m == mDˣ⁺
    @test decay_channel.R13.Γ == ΓDˣ⁺

    # Test branch points calculation
    branch_pts = branch_points(decay_channel)
    @test length(branch_pts) == 2
    @test all(real.(branch_pts) .> 0)
    @test all(imag.(branch_pts) .< 0)  # Branch points should have negative imaginary parts due to decay width
end

@testset "γDD decay channel with Breit-Wigner resonances" begin
    # Define transition matrix elements for γDD
    μ12 = -3.77  # Transition matrix element < D | mu | D* >
    μ13 = +3.77 # Transition matrix element < D | mu | D* >

    # Create decay channel with γDD and two D*+ resonances
    decay_channel = γDD(
        (m1 = 0.0, m2 = mD⁰, m3 = mD⁰),
        BW(m = mDˣ⁺, Γ = ΓDˣ⁺),
        BW(m = mDˣ⁺, Γ = ΓDˣ⁺),
        μ12,
        μ13,
    )

    # Test threshold calculation 1MeV above threshold
    ρ_value = ρ_thr(decay_channel, WithThrE(1.0, mD⁰ + mDˣ⁺) |> e2m)

    # Verify the result is a positive number
    @test ρ_value > 0
    @test typeof(ρ_value) <: Real

    # Test branch points calculation
    branch_pts = branch_points(decay_channel)
    @test length(branch_pts) == 2
    @test all(real.(branch_pts) .> 0)
    @test all(imag.(branch_pts) .< 0)  # Branch points should have negative imaginary parts due to decay width

    # Test that the decay channel can be converted to NamedTuple
    nt = obj2nt(decay_channel)
    @test startswith(nt.type, "γDD")
    @test nt.μ12 == μ12
    @test nt.μ13 == μ13
end
