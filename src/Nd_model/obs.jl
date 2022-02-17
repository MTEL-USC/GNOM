
# TODO separate this and put it into data_processing!
using GEOTRACES
using Statistics
using OceanBasins
using Unitful
using XLSX
using DataFrames
using DataDeps
using Unitful: kg, L



# GEOTRACES data

const DNdobs, εNdobs = let
    # Flierdt et al. (2016) data
    println("Adding Flierdt et al. (2016) data...")
    DNdobs2, εNdobs2 = let
        register(DataDep(
            "Flierdt_etal_2016_Nd",
            """
            Dataset: Global Database from Neodymium in the oceans: a global database, a regional comparison and implications for palaeoceanographic research
            Website: https://figshare.com/articles/Global_Database_from_Neodymium_in_the_oceans_a_global_database_a_regional_comparison_and_implications_for_palaeoceanographic_research/3980064
            Author: Tina van de Flierdt et al.
            Date of Publication: 2016-10-04T07:53:38Z
            License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)

            The neodymium (Nd) isotopic composition of seawater has been used extensively to reconstruct ocean circulation on a variety of timescales. However, dissolved neodymium concentrations and isotopes do not always behave conservatively, and quantitative deconvolution of this non-conservative component can be used to detect trace metal inputs and isotopic exchange at ocean–sediment interfaces. In order to facilitate such comparisons for historical datasets, we here provide an extended global database for Nd isotopes and concentrations in the context of hydrography and nutrients. Since 2010, combined datasets for a large range of trace elements and isotopes are collected on international GEOTRACES section cruises, alongside classical nutrient and hydrography measurements. Here, we take a first step towards exploiting these datasets by comparing high-resolution Nd sections for the western and eastern North Atlantic in the context of hydrography, nutrients and aluminium (Al) concentrations. Evaluating those data in tracer–tracer space reveals that North Atlantic seawater Nd isotopes and concentrations generally follow the patterns of advection, as do Al concentrations. Deviations from water mass mixing are observed locally, associated with the addition or removal of trace metals in benthic nepheloid layers, exchange with ocean margins (i.e. boundary exchange) and/or exchange with particulate phases (i.e. reversible scavenging). We emphasize that the complexity of some of the new datasets cautions against a quantitative interpretation of individual palaeo Nd isotope records, and indicates the importance of spatial reconstructions for a more balanced approach to deciphering past ocean changes.

            Please cite this dataset: van de Flierdt, Tina; Griffiths, Alexander M.; Lambelet, Myriam; Little, Susan H.; Stichel, Torben; J. Wilson, David (2016): Global Database from Neodymium in the oceans: a global database, a regional comparison and implications for palaeoceanographic research. The Royal Society. Dataset. https://doi.org/10.6084/m9.figshare.3980064.v1
            """,
            "https://ndownloader.figshare.com/files/6294558",
            "ba250b1b64b8c43d1130a7aee39674e9ee5936ebb49050c4c224c740bff588b0",
            post_fetch_method = x -> mv(x, "rsta20150293_si_001.xlsx")
        ))
        # Then you can access the file via
        df = XLSX.openxlsx(joinpath(datadep"Flierdt_etal_2016_Nd", "rsta20150293_si_001.xlsx")) do f
            DataFrame(XLSX.gettable(f["Nd isotope database"])...)
        end
        # Note some Nd values are negative, some depths are not... numbers... or too deep...
        # Preprocess and remove them here
        # Also remove arctic and mediterranean obs
        # Actually I just keep what's in PAC, IND, or ATL.
        cleandepth(x::Number) = x
        function cleandepth(x::String)
            out = occursin("+/-", x) ? parse.(Float64, split(x, "+/-")[1]) : # take a in a ± b
            occursin("/", x) ? mean(parse.(Float64, split(x, "/"))) :  # take (a+b)/2 in a/b
            occursin("-", x) ? mean(parse.(Float64, split(x, "-"))) :  # take (a+b)/2 in a–b
            -999.0
            println("  Processing \"$x\" as $out")
            return out
        end
        # Rename columns, convert to floats, clean depth data, filter nonsense, and attach units
        DNdobs2 = let
            X = select(df,
               Symbol("Lat [°N]") => (x->Float64.(x)) => :lat,
               Symbol("Long [°E]") => (x->Float64.(x)) => :lon,
               Symbol("Depth [m]") => (x->cleandepth.(x)) => :depth,
               Symbol("Nd [pmol/kg]") => (x->Float64.(x)) => :Nd,
            )
            X = X[(-999 .< X.Nd .< 75) .& (-999 .< X.depth .< 6000) .& (-999 .< X.lon) .& (-999 .< X.lat), :]
            select(X, :lat, :lon=>(x->mod.(x,360))=>:lon, :depth=>(x->x*u"m")=>:depth, :Nd=>(x->x*u"pmol/kg")=>:Nd)
        end
        εNdobs2 = let
            X = select(df,
               Symbol("Lat [°N]") => (x->Float64.(x)) => :lat,
               Symbol("Long [°E]") => (x->Float64.(x)) => :lon,
               Symbol("Depth [m]") => (x->cleandepth.(x)) => :depth,
               Symbol("εNd") => (x->Float64.(x)) => :εNd
            )
            X = X[(-999 .< X.εNd) .& (-999 .< X.depth .< 6000) .& (-999 .< X.lon) .& (-999 .< X.lat), :]
            select(X, :lat, :lon=>(x->mod.(x,360))=>:lon, :depth=>(x->x*u"m")=>:depth, :εNd=>(x->x*per10000)=>:εNd)
        end
        # Add data source column
        DNdobs2.source = fill("van de Flierdt", size(DNdobs2,1))
        εNdobs2.source = fill("van de Flierdt", size(εNdobs2,1))
        # AIBECS mismatch expects "value" column
        DNdobs2.value = ustrip.(upreferred.(DNdobs2.Nd * ρSW))
        εNdobs2.value = ustrip.(upreferred.(εNdobs2.εNd))
        DNdobs2, εNdobs2
    end

    println("Adding GEOTRACES IDP17 data...")
    DNdobs1, εNdobs1 = let
        DNdobs1 = GEOTRACES.observations("Nd")
        εNdobs1 = GEOTRACES.observations("εNd")
        # Add data source column
        DNdobs1.source = fill("GEOTRACES IDP17", size(DNdobs1,1))
        εNdobs1.source = fill("GEOTRACES IDP17", size(εNdobs1,1))
        # AIBECS mismatch expects "value" column
        DNdobs1.value = ustrip.(upreferred.(DNdobs1.Nd * ρSW))
        εNdobs1.value = ustrip.(upreferred.(εNdobs1.εNd))
        DNdobs1, εNdobs1
    end

    println("Adding post GEOTRACES IDP17 data...")
    DNdobs3, εNdobs3 = let
        # register our compilation of post-IDP17 Nd and εNd data
        register(
            DataDep(
                "posd-IDP17_Nd_data",
                """
                Marine neodymium and epsilon data taken from the literature after GEOTRACES IDP-17.

                Compilation of data from the Indian Ocean (Amakawaet al., 2019), the Barents Sea (Laukert et al., 2018; Lauk-ert et al., 2019), the northern Iceland Basin (Morrison et al.,2019), the Northwestern Pacific (Che and Zhang, 2018), the Kerguelen Plateau (Grenier et al., 2018), the southeastern At-lantic Ocean (GA08, Rahlf et al., 2020; Rahlf et al., 2019;Rahlf et al., 2021; Rahlf et al., 2020), the Bay of Biscay (Dausmann et al., 2020; Dausmann et al., 2019), the Western North Atlantic (Stichel et al., 2020; Stichel et al., 2020), the arctic (Laukert et al., 2017; Laukert et al., 2017a, d), and theBermuda Atlantic Time-series Study (BATS; Laukert et al.,2017; Laukert et al., 2017b, c)

                See Pasquier, Hines, et al. (2021) for more details
                """,
                "https://ndownloader.figshare.com/files/28958076",
                "aba03e69706f5fc88727d4a5762ed970bd4531f9cb3dbc5974acdd6608d56bb3",
                fetch_method = fallback_download,
                post_fetch_method = unpack
            )
        )
        # path to post-IDP17 Nd and εNd data
        data3_dir = datadep"posd-IDP17_Nd_data"
        DNdobs3 = Vector{Float64}[]
        εNdobs3 = Vector{Float64}[]
        lat = Vector{Float64}[]
        lon = Vector{Float64}[]
        depth = Vector{Float64}[]
        εlat = Vector{Float64}[]
        εlon = Vector{Float64}[]
        εdepth = Vector{Float64}[]
        for f in readdir(data3_dir)
            print("  $f")
            Nd_df = XLSX.openxlsx(joinpath(data3_dir, f)) do xf
                DataFrame(XLSX.gettable(xf["Nd_LLD"])...)
            end
            Nd_label = names(Nd_df)[findfirst(occursin.("Nd", names(Nd_df)))]
            Ndlat, Ndlon, Nddepth = Float64.(Nd_df.Lat), Float64.(Nd_df.Lon), Float64.(Nd_df.Depth)
            obs_Nd_unit = uparse(split(Nd_label, '(')[2][1:end-1])
            Nd = Float64.(Nd_df[!, Nd_label])
            obs_Nd_unit ≠ u"pmol/kg" && error("No good unit $obs_Nd_unit")
            ikeep = Nd .< 75
            push!(DNdobs3, Nd[ikeep])
            push!(lat, Ndlat[ikeep])
            push!(lon, Ndlon[ikeep])
            push!(depth, Nddepth[ikeep])
            print(" +$(length(Nd[ikeep])) Nd obs,")

            εNd_df = XLSX.openxlsx(joinpath(data3_dir, f)) do xf
                DataFrame(XLSX.gettable(xf["eNd_LLD"])...)
            end
            εNd_label = names(εNd_df)[findfirst(occursin.("εNd", names(εNd_df)))]
            εNdlat, εNdlon, εNddepth = Float64.(εNd_df.Lat), Float64.(εNd_df.Lon), Float64.(εNd_df.Depth)
            εNd = Float64.(εNd_df[!, εNd_label])
            push!(εNdobs3, εNd)
            push!(εlat, εNdlat)
            push!(εlon, εNdlon)
            push!(εdepth, εNddepth)
            println("and +$(length(εNd)) εNd obs!")
        end
        DNdobs3 = vcat(DNdobs3...)
        εNdobs3 = vcat(εNdobs3...)
        lat = vcat(lat...)
        lon = vcat(lon...)
        depth = vcat(depth...)
        εlat = vcat(εlat...)
        εlon = vcat(εlon...)
        εdepth = vcat(εdepth...)

        # Make a DataFrame
        DNdobs3 = DataFrame(lat=lat, lon=mod.(lon, 360), depth=depth*u"m", Nd=DNdobs3*u"pmol/kg")
        εNdobs3 = DataFrame(lat=εlat, lon=mod.(εlon, 360), depth=εdepth*u"m", εNd=εNdobs3*per10000)
        # Add data source column
        DNdobs3.source = fill("post IDP17", size(DNdobs3,1))
        εNdobs3.source = fill("post IDP17", size(εNdobs3,1))
        # AIBECS mismatch expects "value" column
        DNdobs3.value = ustrip.(upreferred.(DNdobs3.Nd * ρSW))
        εNdobs3.value = ustrip.(upreferred.(εNdobs3.εNd))
        DNdobs3, εNdobs3
    end

    # Concatenate obs from all sources into single DataFrame
    DNdobs = vcat(DNdobs1, DNdobs2, DNdobs3, cols=[:lat, :lon, :depth, :Nd, :value, :source])
    εNdobs = vcat(εNdobs1, εNdobs2, εNdobs3, cols=[:lat, :lon, :depth, :εNd, :value, :source])

    # Show how many obs are in each dataset for both Nd and εNd
    println("Number of observations for each dataset:")
    println(vcat(DataFrame(:variable=>:Nd, (Symbol(x)=>count(DNdobs.source .== x) for x in ("van de Flierdt", "GEOTRACES IDP17", "post IDP17"))...),
                 DataFrame(:variable=>:εNd, (Symbol(x)=>count(εNdobs.source .== x) for x in ("van de Flierdt", "GEOTRACES IDP17", "post IDP17"))...)))

    DNdobs, εNdobs
end

# mask out observations not in the 3 major basins (Mediterranean and Arctic for Nd)
#const OCEANS = oceanpolygons()
#DNdobs = begin
#    ikeep = findall(isatlantic(metadata(DNdobs), OCEANS) .|
#                    isindian(metadata(DNdobs),   OCEANS) .|
#                    ispacific(metadata(DNdobs),  OCEANS))
#    onlykeep(DNdobs, ikeep)
#end
nothing
