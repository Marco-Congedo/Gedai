# must have installed Revise, EEGPlot and GLMakie in main julia environment
using BenchmarkTools

push!(LOAD_PATH, abspath(@__DIR__, ".."))
using Gedai

# NB: all timining exclude HIGH-PASS FILTERING
args = (verbose = false, high_pass = nothing)

example_data = "CAUEEG"
data, srate, labels = read_example_data(example_data);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :prewhite, args...)


@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :prewhite, args...)


@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);
@benchmark denoise($data, srate, labels; threaded=true, refCOV = $refCOV, precomp = $precomp, args...)
# Time  (median): 4.952 s (old code) -> 1.690 (1thr, no precomp, cho) -> 1.854 (4thr, no precomp, cho) -> 1.447 (4thr, no precomp, gevd) -> 1.789 (4thr, precomp, cho)

example_data = "artifact_jumps"
data, srate, labels = read_example_data(example_data);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :prewhite, args...)

@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky,refCOV = $refCOV, precomp = $precomp, args...)

# Time  (median): 3.887 s (old code) -> 2.190 (1thr, no precomp, cho) -> 2.222 (4thr, no precomp, cho) -> 2.227 (4thr, precomp, gevd) -> 2.406 (4thr, precomp, cho)

example_data = "empirical_NOISE_EOG_EMG"
data, srate, labels = read_example_data(example_data);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :prewhite, args...)

@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
# Time (median): 9.8s (old code) -> 2.755 (1thr, no precomp, cho) -> 2.059 (4thr, no precomp, cho) -> 2.136 (4thr, precomp, gevd) -> 2.085 (4thr, precomp, cho)

example_data = "synthetic_bad_channels"
data, srate, labels = read_example_data(example_data);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :prewhite, args...)

@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
# Time (median): ?s (old code) -> 27.138 (1thr, no precomp, cho) -> 19.902 (4thr, no precomp, cho) -> 19.489 (4thr, precomp, gevd) -> 19.544 (4thr, precomp, cho)

example_data = "blinks and bad channels"
data, srate, labels = read_example_data(example_data);
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :gevd, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :invsqrt, args...)
@benchmark denoise($data, srate, labels; threaded=false, gevd_method = :prewhite, args...)

@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :cholesky, args...)
@benchmark denoise($data, srate, labels; threaded=true, gevd_method = :gevd, args...)
# Time (median): ?s (old code) -> 14.981 (1thr, no precomp, cho) -> 15.760 (4thr, no precomp, cho) -> 13.428 (4thr, precomp, gevd) -> 14.784 (4thr, precomp, cho)

