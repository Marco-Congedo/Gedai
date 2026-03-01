# must have installed Revise, EEGPlot and GLMakie in main julia environment
using EEGPlot
using GLMakie, Revise

push!(LOAD_PATH, abspath(@__DIR__, ".."))
using Gedai
                    
example_data = "CAUEEG"
# Other examples:   "artifact_jumps", "empirical_NOISE_EOG_EMG", 
#                   "synthetic_bad_channels", "blinks and bad channels"

# 1. Read Data
data, srate, labels = read_example_data(example_data);

# 2. Denoise Data
clean, data_ref, score, elapsed = denoise(data, srate, labels);

# 3. Plot Data
eegplot(clean, srate, labels; 
            overlay=data_ref, 
            overlay_color = :burlywood, 
            Y = data_ref-clean,
            Y_color = :sienna2,
            X_title = "Input (brown) and Denoised (dark grey) EEG",
            Y_title = "Artifacts (Input - Denoised)",
        )

# for possible colors see: 
# https://juliagraphics.github.io/Colors.jl/stable/namedcolors/

