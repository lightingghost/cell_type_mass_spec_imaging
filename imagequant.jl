using DelimitedFiles, Images, FileIO, ImageDraw

include("imagequantfxns.jl")

conf = Dict(
            "fillsize" => 20,
            "gaussianblur" => 1,
            "σscales" => (2.2.^[2]),
            "minimal_blob_amp" => 0.01, 
             "rev_channeldict" => Dict("dapi" => 1, "gfap" => 2, "ki67" => 3), 
            "channeldict" => Dict(1 => :blue, 2 => :green, 3=> :red)
          )

datadir = joinpath(ENV["HOME"], "Dropbox/Lab/Xiaoai/MSI/WithUpdated2.11.19/UpdatedPixelData")

folders = [x for x in readdir(datadir) if !occursin(r"^[.]", x) && occursin(r"^Y", x)]

finalarray = []
for dir in folders # dir = folders[1]
            fullpath = joinpath(datadir, dir)
            fullpathout = joinpath(datadir,"out", dir)
            mkpath(fullpathout)
            arrayout =process_tifs(dir, fullpath,fullpathout, conf)
			push!(finalarray, arrayout)
end

finalarray = vcat(finalarray...)
finalarray = vcat(finalarray...)
header = ["file" "percent_gfap_dapi" "percent_dapi_alone" "nuclei_count" "GFAP_positive_count" "GFAP_negative_count"]
out = vcat(header,finalarray )
percentfile = joinpath(datadir, "gfap_nongfap_percents_counts.csv")
writedlm(percentfile, out, ',')

