using LinearAlgebra
using Statistics

"""
    run_broadband(eeg_data, srate, refCOV; threaded=false)

Performs the broadband pass of Gedai denoising using `gedai_per_band`.
Broadband parameters: epoch_length=1.0, threshold_type="auto"
"""
gedai_broadband(eeg_data::Matrix{T}, 
                srate::Real, 
                refCOV::SymOrHerm, 
                𝘞::Union{LowerTriangular{T}, SymOrHerm, Matrix{T}, Nothing},
                𝘞⁻ᵗ::Union{UpperTriangular{T}, SymOrHerm, Matrix{T}, Nothing};  
                    epoch_length::Real = 1.0,
                    top_PCs::Int = 3,
                    cov_mean_type::Union{Int, Nothing} = nothing,
                    threshold::Real = 1.0,
                    brent_tol::Real = 0.01,
                    t_range::Tuple{Real, Real} = (-3.0, 12.0),
                    gevd_method::Symbol = :cholesky,
                    threaded::Bool=false) where T<:Real = 
    gedai_per_band( eeg_data, srate, epoch_length, refCOV, 𝘞, 𝘞⁻ᵗ; 
                        top_PCs, 
                        cov_mean_type, 
                        threshold, 
                        brent_tol, 
                        t_range,
                        gevd_method, 
                        threaded)


"""
    gedai_per_band(data, srate, epoch_length, refCOV, threshold_type, threaded=false)

Performs Gedai denoising on a single band (or broadband) signal.
Data: Samples x Channels (T x N)
"""
function gedai_per_band(data::AbstractMatrix{T}, 
                        srate::Real, 
                        epoch_length::Real, 
                        refCOV::SymOrHerm, 
                        𝘞::Union{LowerTriangular{T}, SymOrHerm, Matrix{T}, Nothing},
                        𝘞⁻ᵗ::Union{UpperTriangular{T}, SymOrHerm, Matrix{T}, Nothing};
                            top_PCs::Int = 3,
                            cov_mean_type::Union{Int, Nothing} = nothing,
                            threshold::Real = 1.0,
                            brent_tol::Real = 0.01,
                            gevd_method::Symbol = :cholesky,
                            t_range::Tuple{Real, Real} = (-3.0, 12.0),
                            threaded::Bool=false) where T<:Real

    # Function to process streams
    # stream_data is assumed T x N. Iterate epochs (S x N).
    function process_stream(stream_data::AbstractMatrix{T}) where T<:Real
        n_samps_stream = size(stream_data, 1)
        n_ep = n_samps_stream ÷ epoch_samples

        evals = zeros(T, n_chans, n_ep)
        evecs = zeros(T, n_chans, n_chans, n_ep)

        function process_one_epoch(i)
            start_idx = (i - 1) * epoch_samples + 1
            end_idx = i * epoch_samples
            # Covariance of chunk (S x N). cov returns N x N.
            cov_mat = Symmetric(covmat(@view stream_data[start_idx:end_idx, :]; cov_mean_type))

            if gevd_method ∈ (:cholesky, :invsqrt)
                evals[:, i], evecs[:, :, i] = eigen(Symmetric(𝘞⁻ᵗ' * (cov_mat * 𝘞⁻ᵗ)))
            elseif gevd_method === :gevd
                evals[:, i], evecs[:, :, i] = eigen(cov_mat, refCOV)
            elseif gevd_method === :prewhite
                evals[:, i], evecs[:, :, i] = eigen(cov_mat)
            else
                throw(ArgumentError("broadband.jl: `gevd_method` not valid"))
            end
        end

        if threaded
            Threads.@threads for i in 1:n_ep
                process_one_epoch(i)
            end
        else
            @simd for i in 1:n_ep
                @inbounds process_one_epoch(i)
            end
        end

        return evals, evecs
    end

    n_samples, n_chans = size(data)
    epoch_samples = round(Int, srate * epoch_length)

    # Pad if needed (data is T x N, where T=# samples and N=# channels)
    # NOTE: pad! returns a NEW array when padding is needed; must capture return value.
    data = pad!(data, n_samples, epoch_samples)

    # Create second stream of data (overlapping by 50%)
    shifting = epoch_samples ÷ 2
    shifted_data = shif_data(data, epoch_samples, shifting) 

    # Process Streams
    Eval1, Evec1 = process_stream(data)
    Eval2, Evec2 = process_stream(shifted_data)

   
    # Determine Threshold and sensai_val (using Stream 1). Read value from `threshold` argument
    # In process_stream, we did not return the total covariance matrices,
    # so we compute them quickly here for stream 1 since clean_sensai needs them.
    n_ep1 = size(data, 1) ÷ epoch_samples
    cov_total1 = zeros(T, n_chans, n_chans, n_ep1)
    for i in 1:n_ep1
        start_idx = (i - 1) * epoch_samples + 1
        end_idx = i * epoch_samples
        cov_total1[:, :, i] = Symmetric(covmat(@view data[start_idx:end_idx, :]; cov_mean_type))
    end

    # Determine Threshold and sensai_val (using Stream 1). Read value from `threshold` argument
    function objective(th) # ???
        cov_signal, cov_noise = clean_sensai(th, gevd_method, refCOV, 𝘞⁻ᵗ, Eval1, Evec1; cov_total=cov_total1, threaded)
        score = sensai_metric_fast(cov_signal, cov_noise, refCOV, threshold; top_PCs, threaded)
        return -score
    end

    # Run optimization
    final_threshold, neg_score = local_fminbnd(objective, t_range[1], t_range[2], brent_tol)
    sensai_val = -neg_score
    
    # Clean Streams (Returns T x N)
    clean1, artifact1 = clean_eeg(data, srate, epoch_length, final_threshold, gevd_method, refCOV, 𝘞, 𝘞⁻ᵗ, Eval1, Evec1; 
                                threaded)
    clean2, artifact2 = clean_eeg(shifted_data, srate, epoch_length, final_threshold, gevd_method, refCOV, 𝘞, 𝘞⁻ᵗ, Eval2, Evec2; 
                                threaded)

    # Combine Streams; weights vector (S x 1).
    weights = create_cosine_weights(n_chans, epoch_samples)

    final_clean = merge_streams(clean1, artifact1, clean2, artifact2, weights, shifting, n_samples) 

    return final_clean, sensai_val, final_threshold
end
