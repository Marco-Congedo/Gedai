"""
    denoise()

Entry point for Gedai denoising.
"""
function denoise(data::Matrix{T}, sr::Union{Float64, Int}, labels::Vector{String};
                    wavelet_levels  ::Int                           = 9,
                    high_pass       ::Union{Real, Nothing}          = 0.5,
                    lambda          ::Float64                       = 0.05,
                    epoch_length    ::Float64                       = 1.0,
                    top_PCs         ::Int                           = 3, 
                    cov_mean_type   ::Union{Int, Nothing}           = nothing, # or `0` 
                    threshold       ::Float64                       = 1.0,
                    brent_tol       ::Float64                       = 0.01,
                    t_range         ::Tuple{Float64, Float64}       = (0.0, 12.0),
                    gevd_method     ::Symbol                        = :cholesky,
                    refCOV          ::Union{SymOrHerm, Nothing}     = nothing,
                    precomp         ::Union{Chols, Whits, Nothing}  = nothing,
                    threaded        ::Bool                          = true,
                    BLAS_threaded   ::Bool                          = true,
                    verbose         ::Bool                          = true) where T<:Real

    t_start = time()

    verbose && println(font1color, "Gedai Denoising...", fontwhite)
    
    srate = Float64(sr)
    gevd_method ∈ (:cholesky, :invsqrt, :gevd, :prewhite) || throw(ArgumentError("Package Gedai.jl; function `denoise`: unknown `gevd_method` argument value"))

    multiband = wavelet_levels > 1 # if multiband is true, do Gedai Wavelets

    # Filter data and apply Pseudo Common Average Reference (pseudocar)
    X = isnothing(high_pass) ?  pseudo_car(data) :
                                pseudo_car(filter_data(data, srate, high_pass)) 

    if !(gevd_method === :prewhite) || (gevd_method === :prewhite && isnothing(precomp))
        refCOV = isnothing(refCOV) ? refcov(labels, lambda) : refCOV
    end

    # compute matrices that are used all the time if not provided by argument `precomp`
    𝘞, 𝘞⁻ᵗ = precompute(refCOV, gevd_method; warning = false)
    Xcopy = gevd_method === :prewhite ? deepcopy(X) : nothing

    # We can even 'pre-GEDAI' the data using all data instead of just pre-whitening ,
    # but the SENSAI score methodology must be revised in this case  
    #=
    if gevd_method === :prewhite
        U = eigvecs(Symmetric(𝘞⁻ᵗ * covmat(X; cov_mean_type) * 𝘞⁻ᵗ))
        𝘞⁻ᵗ = 𝘞⁻ᵗ * U'
        𝘞 = U * 𝘞
    end
    =#

    gevd_method === :prewhite && (X *= 𝘞⁻ᵗ)

    try
        # Manage BLAS threads to avoid over-subscription during multi-threaded execution
        BLAS_threads_before = manage_BLAS_threads(BLAS_threaded, threaded)

        # 2. Broadband GEDAI. Use a different threshold if it is the only pass or not
        verbose && println("Broadband Pass...")
        clean, _ = gedai_broadband(X, srate, refCOV, 𝘞, 𝘞⁻ᵗ; 
                                        top_PCs, 
                                        cov_mean_type, 
                                        epoch_length, 
                                        threshold = multiband ? threshold*2 : threshold, 
                                        brent_tol,
                                        t_range,
                                        gevd_method, 
                                        threaded)

        # 3. Wavelet Decomposition (if multiband is true)
        if multiband

            verbose && print("Wavelet Pass... ")
            clean, _ =  gedai_wavelet(clean, srate, refCOV, 𝘞, 𝘞⁻ᵗ;
                                        wavelet_levels, 
                                        high_pass,
                                        top_PCs,
                                        cov_mean_type, 
                                        threshold,
                                        brent_tol,
                                        t_range,
                                        gevd_method,
                                        threaded,
                                        verbose)
        end

        # Calculate global SENSAI on final result; use 1.0s epoch length as standard for scoring
        sensai_out = sensai_metric(clean, X - clean, srate, 1.0, refCOV, threshold;
                                    top_PCs, 
                                    cov_mean_type, 
                                    threaded)


        # Calculate Timing
        elapsed = time() - t_start

        # Logical CPU and threads Info
        verbose && println( font3color, "\nYour (System, Julia, BLAS) Threads: ", fontwhite, "(", 
                            fontwhite, Sys.CPU_THREADS, ",", Threads.nthreads(), ",", LinearAlgebra.BLAS.get_num_threads(), ")")

        # Output Information
        verbose && println(font1color, "\nFinal SENSAI score: ", font2color, round(sensai_out, digits=2),
        fontwhite, "; ", font1color, "Time: ", font2color, round(elapsed, digits=3), "s")


        # Return 4-tuple (clean data, referenced data, sensai score, elapsed time)
        # NB: returning clean, X, ... when gevd_method === :prewhite, X shows the components that have  been removed!
        if gevd_method === :prewhite
            return clean*𝘞, Xcopy, sensai_out, elapsed
        else
            return clean, X, sensai_out, elapsed
        end
    finally
        # Restore the number of BLAS threads as before entering the function
        BLAS_threaded || LinearAlgebra.BLAS.set_num_threads(BLAS_threads_before)
    end
end
