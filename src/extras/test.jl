push!(LOAD_PATH, abspath(@__DIR__, ".."))
using LinearAlgebra, PosDefManifold, Revise, Gedai

example_data = "CAUEEG";
data, srate, labels = read_example_data(example_data);

refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels);

# no multi-threading, no factors pre-computation : 76.43 , 75.66, 77.86, 77.2
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; threaded=false);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; threaded=false, gevd_method = :gevd);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; threaded=false, gevd_method = :invsqrt);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; threaded=false, gevd_method = :prewhite);

# multi-threading, no factors pre-computation : 78.05
clean2, data_ref2, score2, t2 = denoise(data, srate, labels); # by default threaded=true

# no multi-threading, factors pre-computation : 76.43
refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);
clean3, data_ref3, score3, t3 = denoise(data, srate, labels; threaded=false, refCOV, precomp);

# multi-threading, factors pre-computation (This is the fastest) : 78.05
clean4, data_ref4, score4, t4 = denoise(data, srate, labels; refCOV, precomp);

using EEGPlot, GLMakie
args=(overlay_color = :burlywood, Y_color = :sienna2);

eegplot(clean1, srate, labels; overlay=data_ref1, Y = data_ref1-clean1, args...)
eegplot(clean2, srate, labels; overlay=data_ref2, Y = data_ref2-clean2, args...)
eegplot(clean3, srate, labels; overlay=data_ref3, Y = data_ref3-clean3, args...)
eegplot(clean4, srate, labels; overlay=data_ref4, Y = data_ref4-clean4, args...)


# Other example files :

example_data = "artifact_jumps"
data, srate, labels = read_example_data(example_data);

refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels);


clean1, data_ref1, score1, t1 = denoise(data, srate, labels);
clean3, data_ref3, score3, t3 = denoise(data, srate, labels; threaded=false);

# SENSAI SCORE : 78.29

example_data = "empirical_NOISE_EOG_EMG"
data, srate, labels = read_example_data(example_data);

refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels);

clean1, data_ref1, score1, t1 =  denoise(data, srate, labels);
clean2, data_ref2, score2, t2 = denoise(data, srate, labels; refCOV, precomp);
# This will give a slight different result because we do not multithread BLAS
clean3, data_ref3, score3, t3 = denoise(data, srate, labels; threaded=false, BLAS_threaded=false);
clean3, data_ref3, score3, t3 = denoise(data, srate, labels; threaded=false, BLAS_threaded=true);
# SENSAI SCORE : 64.1

example_data = "synthetic_bad_channels"
data, srate, labels = read_example_data(example_data);

refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp);

precomp = precompute(refCOV, :prewhite);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp, gevd_method=:prewhite);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels);



example_data = "blinks and bad channels"
data, srate, labels = read_example_data(example_data);

refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :prewhite);
clean1, data_ref1, score1, t1 = denoise(data, srate, labels; refCOV, precomp, gevd_method=:prewhite);

clean1, data_ref1, score1, t1 = denoise(data, srate, labels);


# Suppose you want to clean a database of many files
# and that you have an array of data matrices `datas`
# all with the same electrode montage (labels) and sampling rate (srate):

clean = similar(data) # to store the cleaned data
refCOV = refcov(labels, 0.05);
precomp = precompute(refCOV, :cholesky);
@threads for (d, data) in enumerate(datas) # multi-threading across files 
            clean[d], rest = denoise(data, srate, labels; threaded=false, refCOV, precomp);
end
