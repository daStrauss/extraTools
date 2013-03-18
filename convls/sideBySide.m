
for mci  = 1
    for m = 1:11
        for n = 1:11
            a = randn(2^m,1);
            b = randn(2^n,1);
            
            %            853.766919 seconds for 2^30 random fft.