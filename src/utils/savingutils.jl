using DrWatson, Arrow, DataFrames

DrWatson._wsave(filename, df::DataFrame; kwargs...) = Arrow.write(filename, df; kwargs...)

function save_kspectrum(spect_data, params; folder="spectra")
    output_stream = datadir(folder, savename("kspectrum", params, "arrow"; sigdigits = 4))
    #println(output_stream)
    df_spect =  DataFrame(k=spect_data.k, ten=spect_data.ten, control=spect_data.control)
    safesave(output_stream, df_spect)
end


function load_kspectrum(params; folder="spectra")
    output_stream = datadir(folder, savename("kspectrum", params, "arrow"; sigdigits = 4))
    #println(output_stream)
    df_spect =  DataFrame(k=spect_data.k, ten=spect_data.ten, control=spect_data.control)
    safesave(output_stream, df_spect)
end
