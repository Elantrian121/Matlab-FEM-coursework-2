function[c] = rootmean_sq(actual,calculated)
err = (actual-calculated);
sq_Err = err.^2
mean_sq_Err = mean(sq_Err)
rootmean_sqa = sqrt(mean_sq_Err)
c = rootmean_sqa
end