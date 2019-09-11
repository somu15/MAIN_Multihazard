function my_x = my_random(myECDF, xi, N)
    % Generate N uniformly distributed samples between 0 and 1.
    u = rand(N,1);
    [myECDF, index] = unique(myECDF);
    
    % Map these to the points on the empirical CDF.
    my_x = interp1(myECDF, xi(index), u, 'linear');
end