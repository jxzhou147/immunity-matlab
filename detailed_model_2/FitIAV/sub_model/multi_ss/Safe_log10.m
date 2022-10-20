function y = Safe_log10(x)
    y = log10(x);
    for i = 1:length(x)
        if x(i) <= 0
            y(i) = 0;
        end
    end
        
end