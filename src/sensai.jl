using LinearAlgebra
using Statistics
using StatsBase

"""
    subspace_angles(A, B)

Calculates principal angles between subspaces defined by orthonormal bases A and B.
Returns angles in radians.
"""
function subspace_angles(A::Matrix{T}, B::Matrix{T}) where T<:Real
    # SVD of A' * B gives singular values = cos(theta)
    λ = svd(transpose(A) * B).S
    angles = acos.(clamp.(λ, -1.0, 1.0)) # Clamp to valid range for acos
    return sort(angles)
end

"""
    create_cosine_weights(T)

Creates a window of cosine weights (T x 1) for broadcasting.
"""
create_cosine_weights(T::Int) = reshape([0.5 - 0.5 * cos(2 * pi * u / T) for u in 1:T], T, 1)

# Backward compatibility
create_cosine_weights(N::Int, T::Int) = create_cosine_weights(T)


"""
    clean_eeg(stream_data, srate, epoch_length, threshold, Eval, Evec)

Reconstructs signal after removing artifactual components based on GEVD.
Input:
    stream_data: (Samples x Channels) matrix [T x N] (or view)
    Eval, Evec: GEVD results (Channels x Channels x Epochs)
    threaded: Bool (enable multi-threading)
Output:
    cleaned_stream: (T x N)
    artifacts_stream: (T x N)
"""
function clean_eeg( stream_data::AbstractMatrix{T},
                    srate::Real,
                    epoch_length::Real,
                    threshold::Real,
                    gevd_method::Symbol,
                    refCOV::SymOrHerm,
                    𝘞::Union{LowerTriangular{T}, SymOrHerm, Matrix{T}, Nothing}, 
                    𝘞⁻ᵗ::Union{UpperTriangular{T}, SymOrHerm, Matrix{T}, Nothing},
                    Eval::Array{T,2},
                    Evec::Array{T,3};
                        threaded::Bool=false) where T<:Real

    n_chans = size(stream_data, 2)
    epoch_samples = Int(round(srate * epoch_length))
    n_epochs = size(Eval, 2)

    # 1. Global Threshold Calculation (PIT)
    cutoff = get_cutoff(Eval, threshold)

    # 2. Clean each epoch
    cleaned_stream = similar(stream_data)
    artifacts_stream = similar(stream_data)

    weights = create_cosine_weights(n_chans, epoch_samples)
    half_epoch = epoch_samples ÷ 2

    function process_epoch(i)

        start_idx = (i - 1) * epoch_samples + 1
        end_idx = i * epoch_samples
        epoch_block = @view stream_data[start_idx:end_idx, :] # S x N

        mask = abs.(Eval[:, i]) .> cutoff # indices of the components to retain

        # Remove brain components (Keep only Artifacts) and back project to sensor space
        if gevd_method ∈ (:cholesky, :invsqrt)
            U = (𝘞⁻ᵗ * Evec[:, :, i][:, mask]) # generalized eigenvectors  
            U⁻¹ = U' * refCOV # inverse of U 
            signal_to_remove = (epoch_block * U)*U⁻¹
        elseif gevd_method === :gevd
            U = Evec[:, :, i][:, mask] # generalized eigenvectors  
            U⁻¹ = U' * refCOV # inverse of U    
            signal_to_remove = (epoch_block * U)*U⁻¹
        elseif gevd_method === :prewhite
            U = Evec[:, :, i][:, mask] # generalized eigenvectors  
            signal_to_remove = (epoch_block * U)*U'
        else
            throw(ArgumentError("sensai.jl: `gevd_method` not valid"))            
        end

        #signal_to_remove = epoch_block * (U * U⁻¹)

        cleaned_epoch = epoch_block - signal_to_remove

        # Apply weights uses broadcasting.
        # Epoch 1: stream2 starts here, so fade second half.
        # Last epoch: stream2 ends here, so fade first half only.
        # Second-to-last epoch (n_epochs-1): stream2's last epoch covers this region with its 
        #   fade-out weights (near zero at the end), so its contribution to the second half 
        #   of epoch n_epochs-1 is near zero — same situation as the last epoch.
        #   Only fade this epoch's first half so the second half isn't attenuated to zero.
        # Middle epochs: both halves covered by stream2, apply full cosine window.
        if i == 1
            cleaned_epoch[(half_epoch+1):end, :] .*= weights[(half_epoch+1):end, :]
        elseif i == n_epochs || i == n_epochs - 1
            cleaned_epoch[1:half_epoch, :] .*= weights[1:half_epoch, :]
        else
            cleaned_epoch .*= weights
        end

        cleaned_stream[start_idx:end_idx, :] = cleaned_epoch
        artifacts_stream[start_idx:end_idx, :] = signal_to_remove 
    end

    if threaded
        Threads.@threads for i in 1:n_epochs
            process_epoch(i)
        end
    else
        for i in 1:n_epochs
            @inbounds process_epoch(i)
        end
    end

    return cleaned_stream, artifacts_stream
end

"""
    sensai_metric(cleaned_data, artifact_data, srate, epoch_length, refCOV, threshold; threaded=false)
"""
function sensai_metric(cleaned_data::AbstractMatrix{T},
                        artifact_data::AbstractMatrix{T},
                        srate::Real, 
                        epoch_length::Real,
                        refCOV::SymOrHerm,
                        threshold::Real;
            top_PCs::Int =3,  
            cov_mean_type::Union{Int, Nothing} = nothing,                              
            threaded::Bool=false) where T<: Real

    n_samples_total, n_chans = size(cleaned_data)
    epoch_samples = Int(floor(srate * epoch_length))
    n_epochs = n_samples_total ÷ epoch_samples

    n_epochs < 1 && (return 0.0)

    # ref_vecs_local = Matrix{Float64}(I, size(refCOV)...)
    ref_vecs_local = princ_eigvecs(refCOV, top_PCs)

    # Pre-allocate for threading
    sig_sims = zeros(Float64, n_epochs)
    noise_sims = zeros(Float64, n_epochs)

    function process_metric(i)
        start_idx = (i - 1) * epoch_samples + 1
        end_idx = i * epoch_samples

        c_block = cleaned_data[start_idx:end_idx, :]
        #=
        
        sig_noise = Symmetric(covmat(c_block; cov_mean_type))
        F_sig = eigen(sig_noise)
        perm_sig = sortperm(F_sig.values, rev=true)
        sig_vecs = F_sig.vectors[:, perm_sig[1:top_PCs]]
        =#
        
        sig_vecs = princ_eigvecs(Symmetric(covmat(c_block; cov_mean_type)), top_PCs)
        angles = subspace_angles(sig_vecs, ref_vecs_local)
        sig_sims[i] = abs(prod(cos.(angles)))

        # Just use eigenvalues!
        #sig_sims[i] = prod( eigvals(Symmetric(covmat(c_block; cov_mean_type)))[end-top_PCs+1:end] )


        a_block = artifact_data[start_idx:end_idx, :]

        #=
        cov_noise = Symmetric(covmat(a_block; cov_mean_type))
        F_noise = eigen(cov_noise)
        perm_noise = sortperm(F_noise.values, rev=true)
        noise_vecs = F_noise.vectors[:, perm_noise[1:top_PCs]]
        =#
        noise_vecs = princ_eigvecs(Symmetric(covmat(a_block; cov_mean_type)), top_PCs)
        angles = subspace_angles(noise_vecs, ref_vecs_local)
        noise_sims[i] = abs(prod(cos.(angles)))
        

        # Just use eigevalues!
        # noise_sims[i] = prod( eigvals(Symmetric(covmat(a_block; cov_mean_type)))[end-top_PCs+1:end] )
    end

    if threaded
        Threads.@threads for i in 1:n_epochs
            process_metric(i)
        end
    else
        @simd for i in 1:n_epochs
            @inbounds process_metric(i)
        end
    end

    avg_sig = 100 * mean(sig_sims)
    avg_noise = 100 * mean(noise_sims)

    return avg_sig - (threshold * avg_noise)
end

"""
    clean_sensai(threshold, gevd_method, refCOV, 𝘞⁻ᵗ, Eval, Evec; threaded=false)

Calculates signal and noise covariance matrices directly from GEVD components.
Returns:
    cov_signal_epoched: Array of signal covariances (Channels x Channels x Epochs)
    cov_noise_epoched: Array of noise covariances (Channels x Channels x Epochs)
"""
function clean_sensai(threshold::Real,
                      gevd_method::Symbol,
                      refCOV::SymOrHerm,
                      𝘞⁻ᵗ::Union{UpperTriangular{T}, SymOrHerm, Matrix{T}, Nothing},
                      Eval::Array{T,2},
                      Evec::Array{T,3};
                      cov_total::Union{Array{T,3}, Nothing}=nothing,
                      threaded::Bool=false) where T<:Real
    
    isnothing(cov_total) && throw(ArgumentError("clean_sensai: `cov_total` must be provided and not nothing"))

    n_chans = size(Eval, 1)
    n_epochs = size(Eval, 2)

    cov_signal_epoched = zeros(T, n_chans, n_chans, n_epochs)
    cov_noise_epoched = zeros(T, n_chans, n_chans, n_epochs)

    cutoff = get_cutoff(Eval, threshold)

    function process_epoch(i)
        # mask = true for components to retain (signal), false for artifact
        mask = abs.(Eval[:, i]) .> cutoff
        artifact_mask = .!mask
        
        # We need U⁻¹ to project the decoupled variances back to sensor space
        if gevd_method ∈ (:cholesky, :invsqrt)
            U_bad = 𝘞⁻ᵗ * Evec[:, :, i][:, mask]
            U⁻¹ = U_bad' * refCOV
        elseif gevd_method === :gevd
            U_bad = Evec[:, :, i][:, mask] 
            U⁻¹ = U_bad' * refCOV
        elseif gevd_method === :prewhite
            U_bad = Evec[:, :, i][:, mask] 
            U⁻¹ = U_bad'
        else
            throw(ArgumentError("sensai.jl: `gevd_method` not valid in clean_sensai"))
        end
        
        # Noise covariance reconstruction:
        # signal_to_remove = data * U_bad * U⁻¹  so cov(signal_to_remove) = U⁻¹' * U_bad' * C * U_bad * U⁻¹
        # where C = cov_total. And cov(cleaned) = (I - P)' * C * (I - P) where P = U_bad * U⁻¹.
        if any(mask)
            projection_matrix = U_bad * U⁻¹  # N x N
            cov_noise = projection_matrix' * Symmetric(cov_total[:, :, i]) * projection_matrix
            cov_noise_epoched[:, :, i] = Symmetric(cov_noise)

            I_minus_P = I - projection_matrix
            cov_signal = I_minus_P' * Symmetric(cov_total[:, :, i]) * I_minus_P
            cov_signal_epoched[:, :, i] = Symmetric(cov_signal)
        else
            # No artifacts to remove: signal covariance = total covariance
            cov_signal_epoched[:, :, i] = Symmetric(cov_total[:, :, i])
        end
    end

    if threaded
        Threads.@threads for i in 1:n_epochs
            process_epoch(i)
        end
    else
        for i in 1:n_epochs
            @inbounds process_epoch(i)
        end
    end

    return cov_signal_epoched, cov_noise_epoched
end

"""
    sensai_metric_fast(cov_signal_epoched, cov_noise_epoched, refCOV, threshold)

Fast version of sensai_metric that directly takes epoch covariance matrices.
"""
function sensai_metric_fast(cov_signal_epoched::Array{T,3},
                            cov_noise_epoched::Array{T,3},
                            refCOV::SymOrHerm,
                            threshold::Real;
                            top_PCs::Int = 3,
                            threaded::Bool=false) where T<:Real
    n_epochs = size(cov_signal_epoched, 3)
    n_epochs < 1 && (return 0.0)

    ref_vecs_local = princ_eigvecs(refCOV, top_PCs)

    sig_sims = zeros(Float64, n_epochs)
    noise_sims = zeros(Float64, n_epochs)

    function process_metric(i)
        cov_sig = Symmetric(cov_signal_epoched[:, :, i])
        sig_vecs = princ_eigvecs(cov_sig, top_PCs)
        angles_sig = subspace_angles(sig_vecs, ref_vecs_local)
        sig_sims[i] = abs(prod(cos.(angles_sig)))

        cov_noise = Symmetric(cov_noise_epoched[:, :, i])
        noise_vecs = princ_eigvecs(cov_noise, top_PCs)
        angles_noise = subspace_angles(noise_vecs, ref_vecs_local)
        noise_sims[i] = abs(prod(cos.(angles_noise)))
    end

    if threaded
        Threads.@threads for i in 1:n_epochs
            process_metric(i)
        end
    else
        @simd for i in 1:n_epochs
            @inbounds process_metric(i)
        end
    end

    avg_sig = 100 * mean(sig_sims)
    avg_noise = 100 * mean(noise_sims)

    return avg_sig - (threshold * avg_noise)
end

