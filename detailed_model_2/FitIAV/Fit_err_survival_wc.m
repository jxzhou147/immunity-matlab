function err = Fit_err_survival_wc( data, par )
    % Error function of survival and weight change
    % data(:, 1): list of weight change
    % data(:, 2): list for survival probability
    
    K = par(1);
    n = par(2);
        
    err = 0;
    for i = 1:length(data)
        err = err + ((1 - data(i, 1) .^ n ./ (K ^ n + data(i, 1) .^ n)) - data(i, 2)) .^ 2;
    end
    
end

