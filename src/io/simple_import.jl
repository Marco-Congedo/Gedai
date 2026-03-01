
using MAT

"""
    load_leadfield(labels, leadfield_path)

Loads and subsets the precomputed leadfield matrix matching the given labels.
"""
function load_leadfield(labels::Vector{String})
    # Ensure the leadfield covariance file exists
    isfile(leadfield_cov_path) ||
        throw(AssertionError("Leadfile covariance file not found. Please reinstall the package"))

    try
        L = matread(leadfield_cov_path)

        # Check required keys
        if !(haskey(L, "electrodes") && haskey(L, "refCOV"))
            @error("Warning: Missing required keys in .mat file.")
            return nothing
        end

        # Build a lookup table for electrode names
        template = Dict(lab => i for (i, lab) in enumerate(lowercase.(L["electrodes"])))

        # Map input labels to indices, erroring if not found
        indices = map(lowercase.(labels)) do lab
            get(template, lab) do
                @error("Warning: Channel $lab not found in the available Leadfield sensor list.")
                return nothing
            end
        end

        indices === nothing && return nothing

        # Subset the Gram/Covariance matrix and return as a Symmetric type
        # return Symmetric(L["refCOV"][indices, indices])
        #return (L["refCOV"][indices, indices])
        return Symmetric(L["refCOV"][indices, indices])

    catch e
        @error("Warning: Failed to load leadfield from .mat file ($e).")
        return nothing
    end
end


function simple_import_set(file_name::String; verbose::Bool=true)

    dataset = matread(file_name)
    # println("DEBUG: Keys in .set file: ", keys(dataset))

    # Check if packaged under "EEG" struct or flat
    dataset = ("EEG" in keys(dataset) && length(keys(dataset)) == 1) ? dataset["EEG"] : dataset

    # Handle .fdt logic if needed, but import_set has it.
    # The user provided a complex script that depends on NeuroAnalyzer types.
    # We just want the raw data matrix and srate for our GEDAI test.

    # Extract Data
    if haskey(dataset, "data")
        raw_val = dataset["data"]
        if raw_val isa String
            verbose && println("DEBUG: 'data' is a String (FDT path?): ", raw_val)
            # FDT file logic...
            fdt_path = joinpath(dirname(file_name), raw_val)
            verbose && println("DEBUG: Constructing FDT path: ", fdt_path)

            n_chans = Int(dataset["nbchan"])
            n_pnts = Int(dataset["pnts"]) # or length(times)
            n_trials = Int(dataset["trials"])

            # Read FDT
            f = open(fdt_path, "r")
            # Float32 usually
            raw_bytes = read(f)
            close(f)
            data_vec = reinterpret(Float32, raw_bytes)
            # Return T x N (samples x channels) for Julia column-major efficiency
            data = collect(transpose(reshape(Float64.(data_vec), n_chans, n_pnts * n_trials)))
        else
            # raw_val is typically N x T from MATLAB, convert to T x N
            data = collect(transpose(Float64.(raw_val)))
        end
    else
        error("Dataset missing 'data' field")
    end

    # Srate
    srate = Float64(dataset["srate"])

    # Chanlocs
    labels = String[]
    if haskey(dataset, "chanlocs") && !isempty(dataset["chanlocs"])
        chanlocs = dataset["chanlocs"]
        # chanlocs is usually a struct array (Dictionary of list of values in MAT.jl?)
        # MAT.jl reads struct array as Dict("field" => CellArray or Matrix)
        if haskey(chanlocs, "labels")
            raw_labels = chanlocs["labels"]
            # Convert to String vector
            if raw_labels isa Matrix
                # If char array?
                # Usually Cell Array of strings in MATLAB -> Vector{Any} or Matrix{Any} in Julia
                labels = string.(vec(raw_labels))
            elseif raw_labels isa Vector
                labels = string.(raw_labels)
            else
                # Fallback - data is T x N, so n_chans = size(data, 2)
                labels = ["Ch$i" for i in 1:size(data, 2)]
            end
        end
    end

    # Fallback if empty - data is T x N, so n_chans = size(data, 2)
    labels = isempty(labels) ? ["Ch$i" for i in 1:size(data, 2)] : labels

    return data, srate, labels
end

function read_example_data(name::String; verbose = true)

    verbose && println(font1color, "LOADING DATA...", fontwhite)
    file = endswith(name, ".set") ? name : name * ".set" 
    filepath = joinpath(abspath(@__DIR__, "..", ".."), "Example Data", file)
    isfile(filepath) || throw(ArgumentError("The example data file name is mispelled. Possible names are 
        `artifact_jumps`, `empirical_NOISE_EOG_EMG`, `synthetic_bad_channels`, `blinks and bad channels` and `CAUEEG`"))
    data, srate, labels = simple_import_set(filepath; verbose=verbose);
    
    if verbose
        println(font2color, "Loaded Data: ", fontwhite, "$(size(data,1)) Samples x $(size(data,2)) Chans @ $srate Hz")
        print(font2color, "Labels ($(length(labels))): ", fontwhite)
        for (i, l) in enumerate(labels)
            i == length(labels) ? println(l) : print(l, ", ")
            # i % 16 == 0 ? println() : nothing
        end
    end
    
    return data, srate, labels
end