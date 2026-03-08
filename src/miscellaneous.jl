using Statistics
using LinearAlgebra
using StatsBase
using Leadfields


"""
 Return the leadfield covariance matrix for sensor in `labels`, with regularization `lambda`
"""
refcov(labels::Vector{String}, lambda::Float64= 0.05) =
    regularize(load_leadfield(labels), lambda)

"""
 Return the leadfield covariance matrix for sensor in `labels`, with regularization `lambda`
"""
function refcov(labels::Vector{String}, reference::String, lambda::Float64= 0.05;
                cov_mean_type::Union{Int, Nothing} = nothing) # or `0` )

    K, _ = leadfield(labels; reference = reference)

    # compute and regularize (lambda) the model covariance matrix of the leadfield as it is done in Gedai.jl
    return regularize(Symmetric(cov(SimpleCovariance(), K'; mean=cov_mean_type)), lambda)
end


"""
    chol_factors(labels)
    chol_factors(refCOV)

Given the labels of the reference covariance matrix, return the 2-tuple with elements:
1. The Cholesky factor L 
2. The Inverse transpose of L
"""
function chol_factors(refCOV::SymOrHerm) 
    if size(refCOV, 1) < 200 # faster by Gaussian elimination (choInv) when BLAS is multi-threaded
        return PosDefManifold.choInv!(LowerTriangular(Matrix(refCOV)))
    else
        LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS)
        L = cholesky(refCOV).L # standard BLAS cholesky factorization
        return L, L \ I
    end
end


"""
    precompute(refCOV)

Given the reference covariance matrix, return the 2-tuple with elements:
2. The Cholesky factor L 
3. The Inverse transpose of L
"""
function precompute(refCOV::SymOrHerm, gevd_method::Symbol; warning::Bool = true)
    gevd_method ∈ (:cholesky, :invsqrt, :gevd, :prewhite) || throw(ArgumentError("Package Gedai.jl; function `precompute`: unknown `gevd_method` argument value"))
    warning && gevd_method === :gevd && @warn("Package Gedai.jl; function `precompute`: using `gevd` as `gevd_method` there is no point of doing precomputations; the 2-tuple (nothing, nothing) is returned ")
    gevd_method === :gevd && return(nothing, nothing)
    
    if gevd_method === :cholesky
        return chol_factors(refCOV)
    elseif gevd_method ∈ (:invsqrt, :prewhite)
        return Tuple(PosDefManifold.pow(Hermitian(refCOV), 0.5, -0.5))
    end 
end



filter_data(X::Matrix{T}, srate::Float64, hp::Float64) where T <: Real = 
        DSP.filtfilt(DSP.digitalfilter(DSP.Highpass(hp), DSP.Butterworth(4); fs=srate), X)

"""
    non_rank_deficient_average_ref(data)

Performs a "Non-Rank Deficient" re-referencing to the average potential.
This method treats the implicit original reference channel (zeros) as part of the average.
Formula: new_data = data - sum(data) / (n_channels + 1)

Input:
  data: (Samples x Channels) matrix (T x N)

Output:
  ref_data: (Samples x Channels) re-referenced data
"""
function non_rank_deficient_average_ref(data::Matrix{Float64})
    n_chans = size(data, 2)
    
    # Calculate the average including the implicit reference (0)
    # Sum across channels (dim 2)
    # Divide by (n_chans + 1)
    avg_pot = sum(data, dims=2) ./ (n_chans + 1)
    
    # Subtract this average from all channels
    return data .- avg_pot
end

pseudo_car = non_rank_deficient_average_ref

# get the first n principal eigenvectors of symmetric or hermitian matrix C
princ_eigvecs(C::SymOrHerm, n) = eigvecs(C)[:, end-n+1:end]


regularize(C::Union{Nothing, SymOrHerm}, lambda::Float64 = 0.05) = 
    isnothing(C) ? nothing : Symmetric((1.0-lambda) * C + Diagonal(repeat([lambda * tr(C)/size(C, 2)], size(C, 2))))

# covariance matrix computations.
# mean = nothing subtracts the mean; mean = 0 doesn't
covmat(X; cov_mean_type = nothing) = StatsBase.cov(SimpleCovariance(), X; mean=cov_mean_type)

# Pad if needed (data is T x N, where T=# samples and N=# channels)
function pad!(data::Matrix{T}, n_samples::Int, epoch_samples::Int) where T<: Real
    remainder = n_samples % epoch_samples
    if remainder != 0
        samples_to_pad = epoch_samples - remainder
        @views padding = data[end-samples_to_pad+1:end, :]
        newdata = similar(data, size(data,1) + samples_to_pad, size(data,2))
        @views newdata[1:end-samples_to_pad, :] .= data
        @views newdata[end-samples_to_pad+1:end, :] .= reverse(padding, dims=1)
        data = newdata
    end
    return data
end

# Create a second stream of data with 50% epoch-length overlapping 
# Exclude first and last half-epochs (shifting). Truncate to exact epochs.
# Returns a view on data. Does not materialize anything
function shif_data(data::Matrix{T}, epoch_samples::Int, shifting::Int) where T<: Real
    nrows = size(data, 1)
    usable = nrows - 2shifting
    n_epochs_2 = usable ÷ epoch_samples
    last_row = shifting + n_epochs_2 * epoch_samples
    first_row = shifting + 1
    @view data[first_row:last_row, :]
end


function merge_streams( clean1::AbstractMatrix{T},
                        artifact1::AbstractMatrix{T},
                        clean2::AbstractMatrix{T},
                        artifact2::AbstractMatrix{T},
                        weights::AbstractMatrix{T},
                        shifting::Int,
                        n_samples::Int) where {T<:Real}

    n_samples_2 = size(clean2, 1)
    # Allocate output sized to n_samples (original, pre-padding)
    final_clean = similar(clean1, n_samples, size(clean1, 2))

    # Fade in/out only if we have enough samples
    if n_samples_2 >= shifting
        @views begin
            # fade in (front)
            w_front = weights[1:shifting, :]
            clean2[1:shifting, :]    .*= w_front

            # fade out (tail)
            tail_start   = n_samples_2 - shifting + 1
            w_tail_start = size(weights, 1) - shifting + 1
            w_tail       = weights[w_tail_start:end, :]

            clean2[tail_start:end, :]    .*= w_tail
        end
    end

    # Initialize output from stream 1 (already cropped to n_samples)
    @views final_clean[:, :]    .= clean1[1:n_samples, :]

    start_idx = shifting + 1
    end_idx   = min(start_idx + n_samples_2 - 1, n_samples)  # clamp to output size
    n_to_add  = end_idx - start_idx + 1
    @views final_clean[start_idx:end_idx, :]    .+= clean2[1:n_to_add, :]

    return final_clean
end


"""
    local_fminbnd(fun, ax, cx, tol=0.001)

Finds the local minimum of function `fun` in interval `[ax, cx]` using Brent's method 
(Golden Section Search + Parabolic Interpolation).
Ported from standard numerical recipes / verified MATLAB implementation.
"""
function local_fminbnd(fun::Function, ax::Real, cx::Real, tol::Real=0.01)
    # Phi is the golden ratio conjugate
    phi = (3 - sqrt(5)) / 2
    
    # Initialize
    a = min(ax, cx)
    b = max(ax, cx)
    v = a + phi * (b - a)
    w = v
    x = v
    
    fv = fun(v)
    fw = fv
    fx = fv
    
    # Golden section step size
    d = 0.0
    e = 0.0
    
    iter = 0
    max_iter = 100
    
    while iter < max_iter
        iter += 1
        
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-10
        tol2 = 2.0 * tol1
        
        # Check for convergence
        if abs(x - xm) <= (tol2 - 0.5 * (b - a))
            break
        end
        
        if abs(e) > tol1
            # Attempt parabolic interpolation
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0
                p = -p
            end
            q = abs(q)
            etemp = e
            e = d
            
            # Check if parabolic step is acceptable
            if abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)
                # Reject parabolic, use Golden Section
                if x >= xm
                    e = a - x
                else
                    e = b - x
                end
                d = phi * e
            else
                # Accept parabolic step
                d = p / q
                u = x + d
                if (u - a) < tol2 || (b - u) < tol2
                   d = sign(xm - x) * tol1
                end
            end
        else
            # Golden Section step
            if x >= xm
                e = a - x
            else
                e = b - x
            end
            d = phi * e
        end
        
        # Perform step
        if abs(d) >= tol1
            u = x + d
        else
            u = x + sign(d) * tol1
        end
        
        fu = fun(u)
        
        # Update book-keeping variables
        if fu <= fx
            if u >= x
                a = x
            else
                b = x
            end
            v = w; fv = fw
            w = x; fw = fx
            x = u; fx = fu
        else
            if u < x
                a = u
            else
                b = u
            end
            if fu <= fw || w == x
                v = w; fv = fw
                w = u; fw = fu
            elseif fu <= fv || v == x || v == w
                v = u; fv = fu
            end
        end
    end
    
    return x, fx
end


function manage_BLAS_threads(BLAS_threaded, threaded)
    BLAS_threads_before = LinearAlgebra.BLAS.get_num_threads()
    pinthreads(:compact) # pack threads as closely as possible onto neighboring cores
    LinearAlgebra.BLAS.set_num_threads(!threaded ? Sys.CPU_THREADS : 
                                        BLAS_threaded ? Sys.CPU_THREADS ÷ 2 : 1)
    return BLAS_threads_before
end


function get_cutoff(Eval::Array{T,2}, threshold::T) where T<: Real

    orig_data = Vector{T}(undef, length(Eval))
    k = 0

    @inbounds @simd for v in Eval
        a = abs(v)
        if a > 0
            k += 1
            orig_data[k] = log(a + 1e-20) + 100
        end
    end

    unique!(sort!(resize!(orig_data, k)))
    T1 = (105 - threshold) / 100.0

    if isempty(orig_data)
        Treshold1 = 100.0
    else
        ecdf_func = ecdf(orig_data)
        probs = ecdf_func(orig_data)
        upper_PIT = 0.95
        outliers = orig_data[probs .> upper_PIT]
        min_outlier = isempty(outliers) ? maximum(orig_data) : minimum(outliers)
        Treshold1 = T1 * min_outlier
    end

    return exp(Treshold1 - 100.0)
end