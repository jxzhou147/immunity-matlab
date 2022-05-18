function [r_s, g_s, b_s] = arr2rgb( s, bin_num, r, g, b )
    
    bin_len = max(s) / bin_num;
    s_sorted = sort(s);
    
    freq = zeros(bin_num, 1);
    
    j = fix(s_sorted(1) / bin_len) + 1;
    if (j < 1)
        j = 1;
    end
    
    for i = 1:length(s)
        if (s_sorted(i) <= j * bin_len)
            freq(j) = freq(j) + 1;
        else
            if j < bin_num
                j = j + 1;
            end
            freq(j) = freq(j) + 1;
        end
    end
    
    r_s = zeros(bin_num, 1);
    g_s = zeros(bin_num, 1);
    b_s = zeros(bin_num, 1);
    
    for i = 1:length(freq)
        r_s(i) = (r + (255 - r) * freq(i) / max(freq)) / 255;
        g_s(i) = (g + (255 - g) * freq(i) / max(freq)) / 255;
        b_s(i) = (b + (255 - b) * freq(i) / max(freq)) / 255;
    end
    
    r_s = flip(r_s);
    g_s = flip(g_s);
    b_s = flip(b_s);
    
end

