include("../common_setup.jl")

# Extra packages required for Nd model only
using MAT
using OceanBasins
using GeoRegions
using NCDatasets
using UnitfulMoles
using UnitfulMoles: molNd


#==================#
# Model parameters #
#==================#
const εunit = u"pertenthousand"
const Ndunit = pM
@initial_value @units @flattenable @limits @description struct Params{Tp} <: AbstractParameters{Tp}
    α_a::Tp            |  1.0   | NoUnits      | true  |  (0,20)  | "Curvature of Nd release enhancement parabola"
    α_c::Tp            | -10.0  | εunit        | true  | (-20,0)  | "Center of Nd release enhancement parabola"
    α_GIC::Tp          |  2.0   | NoUnits      | true  |   (0,∞)  | "Geenland Nd release enhancement"
    σ_ε::Tp            |  3.0   | εunit        | true  |   (0,∞)  | "Per-pixel variance (std) of εNd"
    σ_river::Tp        |  6.0   | Mmol/yr      | true  |   (0,∞)  | "Magnitude of riverine source"
    σ_gw::Tp           |  2.0   | Mmol/yr      | true  |   (0,∞)  | "Magnitude of groundwater source"
    σ_hydro::Tp        |  1.0   | Mmol/yr      | true  |   (0,∞)  | "Magnitude of hydrothermal source"
    ϕ_0::Tp            |  20.0  | pmol/cm^2/yr | true  |   (0,∞)  | "Sedimentary flux at surface"
    ϕ_∞::Tp            |  10.0  | pmol/cm^2/yr | true  |   (0,∞)  | "Sedimentary flux at infinite depth"
    z_0::Tp            | 200.0  | m            | true  |   (0,∞)  | "Sedimentary-flux depth attenuation"
    ε_EAsia_dust::Tp   |  -8.0  | εunit        | true  | (-12,-2) | "εNd of East Asian mineral dust"
    ε_NEAf_dust::Tp    | -12.0  | εunit        | true  | (-15,-9) | "εNd of North East African mineral dust"
    ε_NWAf_dust::Tp    | -12.0  | εunit        | true  | (-15,-9) | "εNd of North West African mineral dust"
    ε_NAm_dust::Tp     |  -8.0  | εunit        | true  | (-12,-4) | "εNd of North American mineral dust"
    ε_SAf_dust::Tp     | -10.0  | εunit        | true  | (-25,-6) | "εNd of South African mineral dust"
    ε_SAm_dust::Tp     |  -3.0  | εunit        | true  | ( -7, 0) | "εNd of South American mineral dust"
    ε_MECA_dust::Tp    |  -2.0  | εunit        | true  | ( -5, 3) | "εNd of Arabian mineral dust"
    ε_Aus_dust::Tp     |  -4.0  | εunit        | true  | ( -7,-1) | "εNd of Australian mineral dust"
    ε_Sahel_dust::Tp   | -12.0  | εunit        | true  | (-15,-9) | "εNd of Sahel mineral dust"
    sol_volc::Tp       |   90.0 | u"percent"   | true  | (0, 100) | "Volcanic Nd solubility"
    sol_EAsia_dust::Tp |   10.0 | u"percent"   | true  | (0, 100) | "EAsia dust Nd solubility"
    sol_NEAf_dust::Tp  |   10.0 | u"percent"   | true  | (0, 100) | "NEAf dust Nd solubility"
    sol_NWAf_dust::Tp  |   10.0 | u"percent"   | true  | (0, 100) | "NWAf dust Nd solubility"
    sol_NAm_dust::Tp   |   10.0 | u"percent"   | true  | (0, 100) | "NAm dust Nd solubility"
    sol_SAf_dust::Tp   |   10.0 | u"percent"   | true  | (0, 100) | "SAf dust Nd solubility"
    sol_SAm_dust::Tp   |   10.0 | u"percent"   | true  | (0, 100) | "SAm dust Nd solubility"
    sol_MECA_dust::Tp  |   10.0 | u"percent"   | true  | (0, 100) | "MECA dust Nd solubility"
    sol_Aus_dust::Tp   |   10.0 | u"percent"   | true  | (0, 100) | "Aus dust Nd solubility"
    sol_Sahel_dust::Tp |   10.0 | u"percent"   | true  | (0, 100) | "Sahel dust Nd solubility"
    ε_volc::Tp         |  10.0  | εunit        | true  | (  0,15) | "εNd of volcanic ash"
    ε_hydro::Tp        |  10.0  | εunit        | true  | (-10,15) | "εNd of hydrothermal vents"
    K_prec::Tp         | 0.01   | NoUnits      | true  |   (0,∞)  | "Precipitation reaction constant"
    f_prec::Tp         | 0.4    | NoUnits      | true  |   (0,1)  | "Fraction of non-buried precipitated Nd"
    w₀_prec::Tp        | 0.7    | km/yr        | false |   (0,∞)  | "Settling velocity of precipitated Nd"
    K_POC::Tp          | 3e13   | NoUnits      | true  |   (0,∞)  | "POC-scavenging reaction constant"
    f_POC::Tp          | 0.78   | NoUnits      | true  |   (0,1)  | "Fraction of non-buried POC-scavenged Nd"
    w₀_POC::Tp         | 40.0   | m/d          | false |   (0,∞)  | "Settling velocity of POC-scavenged Nd"
    K_bSi::Tp          | 3e13   | NoUnits      | true  |   (0,∞)  | "bSi-scavenging reaction constant"
    f_bSi::Tp          |  0.5   | NoUnits      | true  |   (0,1)  | "Fraction of non-buried bSi-scavenged Nd"
    w₀_bSi::Tp         | 714.069| m/d          | false |   (0,∞)  | "Settling velocity of bSi-scavenged Nd"
    K_dust::Tp         | 2e15   | NoUnits      | true  |   (0,∞)  | "Dust-scavenging reaction constant"
    f_dust::Tp         | 0.073  | NoUnits      | true  |   (0,1)  | "Fraction of non-buried dust-scavenged Nd"
    w₀_dust::Tp        | 1.0    | km/yr        | false |   (0,∞)  | "Settling velocity of dust-scavenged Nd"
    τ_ns::Tp           | 100.0  | yr           | true  |   (0,∞)  | "Benthic sink"
end


#========================#
# Constants of the model #
#========================#
# Nd isotopic ratio from the chondritic uniform reservoir (CHUR)
const R_CHUR = 0.512638
# ε = (R / R_CHUR - 1) when ε is `upreferred`
ε(R) = R / R_CHUR - 1
R(ε) = R_CHUR * (ε + 1)


#===============================================#
# Transport operator (Circulation + Scavenging) #
#===============================================#
# Dissolved Cadmium (DNd) is transported by circulation and reversible scavenging
# Reversible scavenging is done by particles, for which we need concentration fields,
# which are dust, POC, bSi, and "prec" (for homogeneous particle concentration)
@enum ScavenginParticle _prec _dust _POC _bSi

# Auxiliary function to normalize fields (such that ∫x⋅dV = 1)
# E.g., this way, I can use σ_k as the total source in Gmol/yr
# TODO: Check that I still need vnormalize!
const volvec = ustrip.(vector_of_volumes(grd))
vnormalize(x) = x ./ (volvec'x)

# POC from Weber and John 2018 (from AWESOME OCIM)
# MAT files are read with the MAT.jl package
# And then regridded to the current circulation
# Note that for the GNOM paper, we used the OCIM2 was used,
# so that regridding is superfluous. But having regridding setup
# will allow easy swap of the circulation!
const AO_path = AO.download_and_unpack()
const POC = let
    ao = matread(joinpath(AO_path, "AWESOME-OCIM-master", "data", "ao.mat"))
    AO_lat = vec(ao["ao"]["lat"])
    AO_lon = vec(ao["ao"]["lon"])
    AO_depth = vec(ao["ao"]["depth"])
    AO_POC_3D = matread(joinpath(AO_path, "AWESOME-OCIM-master", "data", "WJ18", "POC_WJ18.mat"))["POC_WJ18"]
    POC_3D = regrid(AO_POC_3D, AO_lat, AO_lon, AO_depth, grd)
    POC = vectorize(POC_3D, grd)
    vnormalize(POC) # TODO: remove normalization?
end
# Dust from Chien et al available from AIBECS
const DustNd = 40.0u"g"/u"Mg"
const AEOL_Chienetal = let
    s_A_2D = AeolianSources.load("Chien")
    tmp = Any[]
    for k in AeolianSources.Chien_AEROSOLTYPE_NAMES
        v_2D = s_A_2D[k]
        # Take annual mean and permute dims
        v_2D_annual = permutedims(dropdims(mean(v_2D, dims=3), dims=3), (2,1))
        # Regrid to OCIM2 grid
        v_2D_annual_regridded = regrid(v_2D_annual, s_A_2D[:lat], s_A_2D[:lon], grd)
        # Paint the top layer only, with units
        OneHot3rdDim = reshape(grd.depth .== grd.depth[1], (1,1,size(grd)[3]))
        v_3D = (v_2D_annual_regridded * u"kg/m^2/yr" * DustNd / grd.δdepth[1]) .* OneHot3rdDim
        # Convert to moles of Nd and vectorize
        v = vectorize(ustrip.(u"molNd/m^3/s", v_3D), grd)
        push!(tmp, (k, v))
    end
    Dict(tmp)
end
# Dust from Kok et al available from AIBECS
const AEOL_Koketal = let
    s_A_2D = AeolianSources.load("Kok")
    tmp = Any[]
    for r in AeolianSources.Kok_REGIONS_NAMES
        v_2D = s_A_2D[r]
        # Permute dims
        v_2D_annual = permutedims(v_2D, (2,1))
        # Regrid to OCIM2 grid
        v_2D_annual_regridded = regrid(v_2D_annual, s_A_2D[:lat], s_A_2D[:lon], grd)
        # Paint the top layer only, with units
        OneHot3rdDim = reshape(grd.depth .== grd.depth[1], (1,1,size(grd)[3]))
        v_3D = (v_2D_annual_regridded * u"kg/m^2/yr" * DustNd / grd.δdepth[1]) .* OneHot3rdDim
        # Convert to moles of Nd and vectorize
        v = vectorize(ustrip.(u"molNd/m^3/s", v_3D), grd)
        push!(tmp, (r, v))
    end
    Dict(tmp)
end
const dustparticles = let # Repeat the dust field throughout the water column
    dust_tot_Kok = sum(v for (r,v) in pairs(AEOL_Koketal))
    dust3D = repeat(rearrange_into_3Darray(dust_tot_Kok, grd)[:,:,1], outer=(1,1,24))
    vnormalize(vectorize(dust3D, grd))
end
# Opal (bSi) from side Si-cycle model
# TODO: Add option to run it here, so that changes can be made to the grid.
const bSi = jldopen(joinpath(output_path, "optimized_Simodel_$circname.jld2")) do file
    vnormalize(file["PSi"])
end
# sparsitty pattern of scavenging transport operators
const iscav, jscav, _ = findnz(buildIabove(grd) + I)
# Instead of simply building the matrices with tranposrtoperator, to speed things up, these are constant
# vectors of the values that go into the matrices. This is much more complicated now, but it is much faster, too
# so you know...
const v_scav_noremin_dict, v_scav_justremin_dict = let
    ijscav = findall(!iszero, buildIabove(grd) + I)
    δz = ustrip.(grd.δz_3D[iswet(grd)])
    Iabove = AIBECS.buildIabove(grd.wet3D, findall(vec(iswet(grd))))
    w_base(z) = 1.0
    T_w_noremin = transportoperator(grd, w_base; z_top=z_top, z_bot=z_bot, frac_seafloor=f_topo, Iabove=Iabove, δz=δz, fsedremin=false)
    T_w_fullremin = transportoperator(grd, w_base; z_top=z_top, z_bot=z_bot, frac_seafloor=f_topo, Iabove=Iabove, δz=δz, fsedremin=true)
    T_w_justremin = T_w_fullremin - T_w_noremin
    diagonals = Dict(
        _prec => I,
        _dust => Diagonal(dustparticles),
        _POC => Diagonal(POC),
        _bSi => Diagonal(bSi)
    )
    v_scav_noremin = Dict((t, collect((T_w_noremin * diagonals[t])[ijscav])) for t in instances(ScavenginParticle))
    v_scav_justremin = Dict((t, collect((T_w_justremin * diagonals[t])[ijscav])) for t in instances(ScavenginParticle))
    # Check that the transport operators are going to be correct
    for t in instances(ScavenginParticle)
        (K, w₀, f) = rand(3)
        (K * w₀ * (v_scav_noremin[t] + f * v_scav_justremin[t]) ≈ collect((transportoperator(grd, z -> w₀ ; z_top=z_top, z_bot=z_bot, frac_seafloor=f_topo, fsedremin=f) * K * diagonals[t])[ijscav])) || error("Nah not good for $t")
    end
    v_scav_noremin, v_scav_justremin
end
Kwf(t, p) = AIBECS.UnPack.unpack(p, Val(Symbol(:K, t))),
             AIBECS.UnPack.unpack(p, Val(Symbol(:w₀, t))),
             AIBECS.UnPack.unpack(p, Val(Symbol(:f, t)))
function v_scav(t, p)
    K, w₀, f = Kwf(t, p)
    K .* w₀ .* (v_scav_noremin_dict[t] .+ f .* v_scav_justremin_dict[t])
end
# Using LinearOperators
const colptrT0, rowvalsT0 = let
    T0 = buildIabove(grd) + I
    SparseArrays.getcolptr(T0), rowvals(T0)
end
T_D(t, p) = SparseMatrixCSC(nb, nb, colptrT0, rowvalsT0, v_scav(t, p))
T_D(p) = LinearOperators((T, T_D(_prec, p), T_D(_dust, p), T_D(_POC, p), T_D(_bSi, p)))

#=========================#
# Local sources and sinks #
#=========================#
#
# - sources
#     - from dust
#     - from rivers
#     - a boundary source on the continental margin
# - sinks
#     - reversible scavenging
#     "The particles involved in the reversible scavenging of Nd are CaCO3, opal, POC, and dust."
#         Gu et al.

# The total source is the sum of all sources
s_tot(p) = sum(sₖ(p) for sₖ in (s_dust, s_aeol, s_sed, s_river, s_gw, s_hydro))
s_tot_iso(p) = sum(sₖ(p) for sₖ in (s_dust_iso, s_aeol_iso, s_sed_iso, s_river_iso, s_gw_iso, s_hydro_iso))

# Aeolian sources

# 1. Dust sources
# Load dust-deposition fields from Hamilton data
include(joinpath(root_path, "src", "Julia", "data_processing", "Hamilton_dust_dep.jl"))
const σ_raw_dust = sum(∫dV(yearly_dust_dep(r, grd), grd) for r in DUST_SOURCE_REGIONS[!,:region])
# Julia way of regions
const DUST_REGIONS = AeolianSources.Kok_REGIONS_NAMES
ε_dust(r, p) = AIBECS.UnPack.unpack(p, Val(Symbol(:ε_, r, :_dust)))
α_dust(r, p) = AIBECS.UnPack.unpack(p, Val(Symbol(:α_, r, :_dust)))
ε_dust(p) = [ε_dust(r, p) for r in DUST_REGIONS]
α_dust(p) = [α_dust(r, p) for r in DUST_REGIONS]


σ_dust(p) = AIBECS.UnPack.unpack(p, Val(:σ_dust)) # only one σ_dust
to_dust_region_string(r) = filter(:reg2 => isequal(Symbol(r)), DUST_SOURCE_REGIONS)[1,:region]
s_dust_normalized_tmp(r) = ustrip.(yearly_dust_dep(to_dust_region_string(r), grd) / σ_raw_dust)
const s_dust_normalized_arr = reduce(hcat, s_dust_normalized_tmp(r) for r in instances(DustRegion))
const s_dust_tot_arr = dropdims(sum(s_dust_normalized_arr, dims=2), dims=2)
#s_dust(p) = s_dust_tot_arr * σ_dust(p)
s_dust(p) = s_dust_normalized_arr * (α_dust(p) * σ_dust(p))
s_dust_iso(p) = s_dust_normalized_arr * (R.(ε_dust(p)) .* α_dust(p) * σ_dust(p))
# not fast but for diagnosis:
s_dust(r, p) = σ_dust(p) * α_dust(r, p) * view(s_dust_normalized_arr, :, Int(r))
s_dust_iso(r, p) = R(ε_dust(r, p)) * s_dust(r, p)


# 2. Volc source
# Each aeolian source
@enum AerosolType volc=1 # Just volcano used here
const σ_t_aeol_val = [Val(Symbol(:σ_, t)) for t in instances(AerosolType)]
σ_aeol(t, p) = AIBECS.UnPack.unpack(p, σ_t_aeol_val[Int(t)])
v_aeol(t) = AEOL_Chienetal[Symbol(t)]
s_aeol(t, p) = σ_aeol(t, p) * v_aeol(t)
s_aeol(p) = sum(s_aeol(t, p) for t in instances(AerosolType))
# isotopes
ε_aeol(t, p) = AIBECS.UnPack.unpack(p, Val(Symbol(:ε_, t)))
s_aeol_iso(t, p) = R(ε_aeol(t, p)) * s_aeol(t, p)
s_aeol_iso(p) = sum(s_aeol_iso(t, p) for t in instances(AerosolType))

# 3. Sedimentary source
# Using the Robinson et al data for sedimentary εNd
const ε_sed = let
    eNdmap_data_path, k = joinpath(data_path, "globalSedimentEpsilon.nc"), "globalsedimentepsilonnd"
    #eNdmap_data_path, k = joinpath(data_path, "seafloorEpsilonNd.nc"), "seafloor_end" # <- "raw" seafloor eNd
    lat, lon, eNd = Dataset(eNdmap_data_path, "r") do ds
        @show keys(ds)
        @info "Reordering lat because Robinson data goes from +90 down to -90 and that's bad for Plots.jl"
        lat = ds["lat"][:]
        ilat = sortperm(lat)
        lat = sort(lat)
        lon, eNd = ds["lon"][:], ds[k][:,ilat]
        @info "permuting dims of eNd"
        lat, lon, permutedims(eNd)
    end
    @info "painting eNd from Robinson et al to fill in missing values"
    eNd_painted = inpaint(eNd)
    @info "Regridding eNd from Robinson et alto $Circulation grid"
    eNd_paintregrid = regrid(eNd_painted, lat, lon, grd)
    @info "Copying eNd from Robinson et al at all depths"
    eNd_allz = repeat(eNd_paintregrid, outer=(1,1,size(grd)[3]))
    @warn "Imposing minimum -5 Pacific value"
    eNd_vec = vectorize(eNd_allz, grd)
    OCN = oceanpolygons()
    idx = ispacific2(latvec(grd), lonvec(grd), OCN) .& (eNd_vec .< -5)
    eNd_vec[idx] .= -5.0
    ustrip.(upreferred.(eNd_vec * εunit))
end

# bundle the 1/δz and f_topo together to convert from
# `s_sed_per_area` to `s_sed` (which is per volume)
const v_sed_multiplier = f_topo ./ ustrip.(vectorize(grd.δz_3D, grd))

# Reactivity α as a quadratic function of ε
const α_scaling = upreferred(10εunit)
function α_quad(ε, p)
    @unpack α_a, α_c = p
    @. α_a * ((ε - α_c) / α_scaling)^2 + 1
end
α_quad(p) = α_quad(ε_sed, p)
# Enhanced Greenland Nd release
const GIC_mask = let
    GREENLAND = GeoRegion("AR6_GIC")
    BitVector([isPointinGeoRegion(Point2(lon,lat), GREENLAND, throw=false) for (lon,lat) in zip(lonvec(grd),latvec(grd))])
end
function α_GIC(p)
    @unpack α_GIC = p
    @. α_GIC * GIC_mask + 1.0 * !GIC_mask
end
# Shifting of effective ε released as per notebook in extras/
shifted_ε(μ, σ, a, c, α_scaling) = (a * (μ - 2 * (c-μ)) * σ^2 + (a * (c-μ)^2 + α_scaling^2) * μ)/(a * σ^2 + a * (c - μ)^2 + α_scaling^2)
function shifted_ε_sed(p)
    @unpack α_a, α_c, σ_ε = p
    shifted_ε.(ε_sed, σ_ε, α_a, α_c, α_scaling)
end
R_sed(p) = R.(shifted_ε_sed(p))

function s_sed_per_area_fun(p)
    @unpack ϕ_0, ϕ_∞, z_0 = p
    z -> (ϕ_0 - ϕ_∞) * exp(-z / z_0) + ϕ_∞
end
s_sed_per_area(p) = s_sed_per_area_fun(p).(z_bot)

s_sed(p) = α_quad(p) .* α_GIC(p) .* v_sed_multiplier .* s_sed_per_area(p)
s_sed_iso(p) = R_sed(p) .* α_quad(p) .* α_GIC(p) .* v_sed_multiplier .* s_sed_per_area(p)

# 4. Hydrothermal source
# In case the circulation is not OCIM2,
# load the He flux from OCIM2
# and regrid it to the circulation.
const ϕ_He = let
    grd_OCIM2, _, ϕ_He, _ = OCIM2.load()
    (debug || Circulation ≠ OCIM2) ? ϕ_He = regrid(ϕ_He, grd_OCIM2, grd) : ϕ_He
end
const v_hydro = vnormalize(@. ϕ_He * (z_top > 0))
function s_hydro(p)
    @unpack σ_hydro = p
    return σ_hydro * v_hydro
end
function s_hydro_iso(p)
    @unpack ε_hydro = p
    return R(ε_hydro) * s_hydro(p)
end

# 5. Riverine source
# normalized river source
const s_river0 = let
    RIVERS = Rivers.load()
    rivers = regrid(RIVERS, grd)
    S = smooth_operator(grd, T)
    vnormalize(S * rivers) # S smoothes the singular river source points
end
# Riverine source scaled by global magnitude σ
function s_river(p)
    @unpack σ_river = p
    return σ_river .* α_quad(p) .* s_river0
end
# Isotope river source
s_river_iso(p) = R_sed(p) .* s_river(p)

# 5. Groundwater source
# Groundwater source scaled by global magnitude σ
const s_gw0 = let
    gw_raw = GroundWaters.load()
    gw = ustrip.(u"m^3/s", regrid(gw_raw, grd))
    vnormalize(gw)
end
function s_gw(p)
    @unpack σ_gw = p
    return σ_gw * s_gw0
end
# Isotope groundwater source
const R_gw = let
    Z2D = matread(joinpath(data_path, "Jeandel_DIVAnd_interpolated_topcore_1200km.mat"))["interpolated_Jeandel_eNd"]
    Z3D = zeros(size(grd))
    Z3D[:,:,1] .= Z2D
    εNd_gw_tmp = vectorize(Z3D, grd)
    R.(ustrip.(upreferred.(εNd_gw_tmp * εunit)))
end
s_gw_iso(p) = R_gw .* s_gw(p)


#================================================
Bottom sink
================================================#
function neph_sink(x, p)
    @unpack τ_ns = p
    @. f_topo * x / τ_ns
end



#====================================#
# Right-hand side of tracer equation #
#====================================#
G_Nd(DNd, ⁱDNd, p) = s_tot(p) - neph_sink(DNd, p)
G_Nd_iso(DNd, ⁱDNd, p) = s_tot_iso(p) - neph_sink(ⁱDNd, p)
#G_Nd(dx, DNd, ⁱDNd, p) = dx .= s_tot(p)
#G_Nd_iso(dx, DNd, ⁱDNd, p) = dx .= s_tot_iso(p)

Gs = (G_Nd, G_Nd_iso)

#===============#
# Problem setup #
#===============#
# Initialize parameters with initial values
p = Params()
# initial guess
x = ustrip.(u"mol/m^3", 10u"pM") * ones(nb)
x = [x; x]
# state function and its Jacobian
fun = AIBECSFunction(T_D, Gs, nb, Params)
F, ∇ₓF = F_and_∇ₓF(fun)
# problem
prob = SteadyStateProblem(fun, x, p)



