# Take WOA18 for Si obs
const Si_obs = let
    Si_obs_3D, σ²Si_obs_3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "silicate", 0, "1°", "an")
    Si_obs_3D_painted = inpaint(Si_obs_3D, 0)
    Si_obs_3D_painted[iswet(grd)] * upreferred(mol/g) * ρSW .|> upreferred .|> ustrip
end
