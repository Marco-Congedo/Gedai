using Serialization
using NPZ

# Define the path to the .jls file
jls_file = joinpath(@__DIR__, "io", "fsavLEADFIELD_4_GEDAI.jls")
npz_file = joinpath(@__DIR__, "io", "fsavLEADFIELD_4_GEDAI.npz")

# Check if the .jls file exists
if !isfile(jls_file)
    println("Error: .jls file not found at $jls_file")
    exit(1)
end

try
    # Deserialize the .jls file
    println("Reading .jls file: $jls_file")
    data = open(jls_file, "r") do io
        deserialize(io)
    end

    # Ensure the data contains the required keys
    if !haskey(data, "electrodes") || !haskey(data, "refCOV")
        println("Error: .jls file does not contain required keys 'electrodes' and 'refCOV'")
        exit(1)
    end

    # Extract the data
    electrodes = data["electrodes"]
    refCOV = data["refCOV"]

    # Write the data to a .npz file
    println("Writing to .npz file: $npz_file")
    NPZ.npzwrite(npz_file, Dict(
        "electrodes" => electrodes,
        "refCOV" => refCOV
    ))
    println("Conversion successful: $npz_file")
catch e
    println("Error during conversion: $e")
    exit(1)
end