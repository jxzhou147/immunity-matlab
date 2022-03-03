function root = RealRootPromise( x, n )
    % avoid imagine root for x^n
    
    if (n ~= fix(n)) && (x < 0)
        root = 0;
    else
        root = x ^ n;
    end
    
end

