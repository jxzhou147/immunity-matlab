function x = FracNoInf(a, b)
    % Deal with inf when a / 0
    if b == 0
        x = 0;
    else
        x = a / b;
    end
end