using Statistics
using LinearAlgebra

"""
    squeeze_dims(A)

Drops all singleton dimensions from array A.
"""
function squeeze_dims(A)
    return dropdims(A, dims=tuple(findall(size(A) .== 1)...))
end

"""
    modwt_custom(data, wavelet_type, level)

Custom implementation of MODWT (Maximal Overlap Discrete Wavelet Transform).
Currently supports 'haar'.

Input:
  data: (Samples x Channels)
  wavelet_type: "haar"
  level: Int

Output:
  wpt: (Level+1 x Samples x Channels) array
       Order: [W1, W2, ..., WJ, VJ] (Details ... Approx)
"""
function modwt_custom(data::Matrix{T}, wavelet_type::String, level::Int) where T<:AbstractFloat
    if lowercase(wavelet_type) != "haar"
        error("modwt_custom currently only supports 'haar' wavelet.")
    end

    n_samples, n_channels = size(data)
    inv_sqrt2 = 1.0 / sqrt(2.0)

    # Pre-allocate output: (Level+1, Samples, Channels)
    # Keeping MATLAB's dimension order (Level, Samples, Channels)
    wpt = zeros(T, level + 1, n_samples, n_channels)

    current_approx = copy(data)

    for j in 1:level
        step = 2^(j - 1)

        # Circular shift (along dim 1 = samples)
        shifted_approx = circshift(current_approx, (step, 0))

        # Compute Approx (Low Pass): (x[n] + x[n-step]) * c
        next_approx = (current_approx .+ shifted_approx) .* inv_sqrt2

        # Compute Detail (High Pass): (x[n-step] - x[n]) * c
        detail = (shifted_approx .- current_approx) .* inv_sqrt2

        # Store Detail Wj
        wpt[j, :, :] = detail

        # Update
        current_approx = next_approx
    end

    # Store Final Approx VJ
    wpt[level+1, :, :] = current_approx

    return wpt
end

"""
    modwtmra_custom(wpt, wavelet_type)

Inverse Stationary Wavelet Transform (MRA).

Input:
  wpt: (Bands x Samples x Channels)
       Order: [W1, ..., WJ, VJ]

Output:
  mra: (Bands x Samples x Channels)
       Reconstructed time-domain signal for each band.
"""
function modwtmra_custom(wpt::Array{T,3}, wavelet_type::String) where T<:AbstractFloat
    if lowercase(wavelet_type) != "haar"
        error("modwtmra_custom currently only supports 'haar' wavelet.")
    end

    n_bands, n_samples, n_channels = size(wpt)
    level = n_bands - 1
    inv_sqrt2 = 1.0 / sqrt(2.0)

    mra = zeros(T, n_bands, n_samples, n_channels)

    for band_idx in 1:n_bands
        # Initialize Approx at Level J
        if band_idx == n_bands
            current_approx_recon = wpt[n_bands, :, :]
        else
            current_approx_recon = zeros(T, n_samples, n_channels)
        end

        # Iterate backwards
        for j in level:-1:1
            step = 2^(j - 1)

            if band_idx == j
                D_j = wpt[j, :, :]
            else
                D_j = zeros(T, n_samples, n_channels)
            end

            # Reconstruction formula for Haar
            A_shifted = circshift(current_approx_recon, (-step, 0))
            D_shifted = circshift(D_j, (-step, 0))

            current_approx_recon = 0.5 * inv_sqrt2 .* ((current_approx_recon .+ A_shifted) .+ (D_shifted .- D_j))
        end

        mra[band_idx, :, :] = current_approx_recon
    end

    return mra
end


"""
    gedai_wavelet(data, srate, refCOV, threshold, threaded=false)

Performs the Gedai wavelet pipeline:
1. MODWT Decomposition (Haar)
2. Per-band denoising using `gedai_per_band` (excluding bands < 0.5 Hz)
3. Reconstruction using MODWTMRA
"""
function gedai_wavelet( data::Matrix{T}, 
                        srate::Real, 
                        refCOV::SymOrHerm, 
                        𝘞::Union{LowerTriangular{T}, SymOrHerm, Matrix{T}, Nothing}, 
                        𝘞⁻ᵗ::Union{UpperTriangular{T}, SymOrHerm,Matrix{T},  Nothing};
                            wavelet_levels::Int=9,
                            high_pass::Union{Real, Nothing} = 0.5,
                            top_PCs::Int = 3,  
                            cov_mean_type::Union{Int, Nothing} = nothing,
                            threshold::Real = 1.0,
                            brent_tol::Real = 0.01,
                            t_range::Tuple{Real, Real} = (-3.0, 12.0),
                            gevd_method::Symbol = :cholesky,
                            threaded::Bool=false,
                            verbose::Bool=true) where T<: Real

    try    

        num_of_wavelet_levels = wavelet_levels

        num_of_wavelet_bands = num_of_wavelet_levels
        wavelet_type = "haar" # Only Haar supported

        # 1. Decomposition
        #verbose && println("Performing MODWT (Level $num_of_wavelet_bands)...")
        # modwt_custom expects (Samples x Channels), data is already T x N
        wpt_EEG = modwt_custom(data, wavelet_type, num_of_wavelet_bands)

        # 2. Reconstruction (MRA)
        #verbose && println("Calculating MRA...")
        wpt_mra = modwtmra_custom(wpt_EEG, wavelet_type)

        n_bands, n_samples_trace, n_chans_trace = size(wpt_mra)

        # --- Frequency Calculation & Exclusion Logic (Matching MATLAB) ---
        center_freq = zeros(T, n_bands)
        lower_freq = zeros(T, n_bands)
        upper_freq = zeros(T, n_bands)

        # Fixed epoch size logic (Adaptive)
        num_cycles_in_window = 12.0
        epoch_lengths = zeros(Float64, n_bands)

        if verbose
            println("considering bands:", fontgrey)
            println("   ―――――――――――――――――――――――――――――――――――――――――――")
            println("   Band | Freq Range (Hz)   | Epoch Length (s)")
            println("   ―――――|―――――――――――――――――――|―――――――――――――――――")  
        end

        for f in 1:n_bands
            # The passband for MRA band 'f' is approx [Fs/(2^(f+1)), Fs/(2^f)]

            upper_bound = srate / (2^f)
            calc_lower_bound = srate / (2^(f + 1)) # Used for epoch calculation

            lower_bound = f == n_bands ?  0.0 : calc_lower_bound

            center_freq[f] = (lower_bound + upper_bound) / 2
            lower_freq[f] = lower_bound
            upper_freq[f] = upper_bound

            # Calculate Adaptive Epoch Size
            # Formula: Cycles / LowerBound (using the formula bound to avoid Inf)
            raw_epoch_length = num_cycles_in_window / calc_lower_bound

            # Enforce Even Number of Samples (Matlab parity)
            ideal_samples = raw_epoch_length * srate
            rounded_samples = round(ideal_samples)
            if rem(rounded_samples, 2) != 0
                if abs(ideal_samples - (rounded_samples - 1)) < abs(ideal_samples - (rounded_samples + 1))
                    final_samples = rounded_samples - 1.0
                else
                    final_samples = rounded_samples + 1.0
                end
            else
                final_samples = rounded_samples
            end

            epoch_lengths[f] = final_samples / srate

            # Print Table Row
            if verbose
                f_str = rpad(f, 4)
                freq_str = rpad("$(round(lower_bound, digits=2)) - $(round(upper_bound, digits=2))", 17)
                ep_str = "$(round(epoch_lengths[f], digits=2))"
                println("   $f_str | $freq_str | $ep_str")
            end  
        end
        verbose && println("   ―――――――――――――――――――――――――――――――――――――――――――", fontwhite)

        # Count how many bands have UPPER FREQ <= 0.5 Hz
        lowest_wavelet_bands_to_exclude = count(x -> x <= high_pass, upper_freq)

        num_bands_to_process = n_bands - lowest_wavelet_bands_to_exclude

        if lowest_wavelet_bands_to_exclude > 0
            verbose && println("Excluding $lowest_wavelet_bands_to_exclude bands (< $(high_pass)Hz). Processing 1 to $num_bands_to_process...")
        else
            verbose && println("Processing all $n_bands bands.")
        end

        # 3. Denoise Each Band
        #verbose && println("Denoising $num_bands_to_process bands...")

        cleaned_mra = zeros(Float64, n_bands, n_samples_trace, n_chans_trace)
        sensai_scores = zeros(Float64, n_bands)
        thresholds = zeros(Float64, n_bands)

        function process_band(f)
            # wpt_mra[f,:,:] is (n_samples, n_chans) = T x N
            # squeeze_dims ensures we get T x N
            band_data = squeeze_dims(wpt_mra[f, :, :])
            current_epoch_length = epoch_lengths[f]

            # Check if data length is sufficient for at least one epoch
            if current_epoch_length * srate > size(band_data, 1)
                # println("Warning: Band $f epoch size $(round(current_epoch_length, digits=2))s > Data duration. Skipping (setting to 0).")
                cleaned_mra[f, :, :] .= 0.0
                sensai_scores[f] = 0.0 # Or NaN?
                thresholds[f] = 0.0
                return
            end

            # Denoise (cl is T x N now)
            cl, sen, th = gedai_per_band(band_data, srate, current_epoch_length, refCOV, 𝘞, 𝘞⁻ᵗ; 
                                    top_PCs, 
                                    cov_mean_type,
                                    threshold,
                                    brent_tol,
                                    t_range,
                                    gevd_method,
                                    threaded = false) # in case, multi-thread across wavelet bands, not here

            # cleaned_mra[f, :, :] = cl (T x N matches)
            cleaned_mra[f, :, :] = cl
            sensai_scores[f] = sen
            thresholds[f] = th
        end

        # Loop ONLY up to num_bands_to_process
        if threaded
            Threads.@threads for f in 1:num_bands_to_process
                process_band(f)
            end
        else
            @simd for f in 1:num_bands_to_process
                @inbounds process_band(f)
            end
        end

        # 4. Final Reconstruction (Sum of bands)
        # Summing cleaned_mra includes the zeros for excluded bands, effectively removing them.
        # cleaned_mra is (n_bands, n_samples, n_chans), sum gives (1, n_samples, n_chans)

        # dropdims(sum..., dims=1) results in (n_samples, n_chans) -> T x N
        processed_data = dropdims(sum(cleaned_mra, dims=1), dims=1)

        return processed_data, sensai_scores, thresholds
    catch e
        println("ERROR IN WAVELET PIPE:")
        Base.showerror(stdout, e, catch_backtrace())
        println()
        rethrow(e)
    end

end
