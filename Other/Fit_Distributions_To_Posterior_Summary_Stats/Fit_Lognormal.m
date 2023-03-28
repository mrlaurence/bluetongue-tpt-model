function param_info = Fit_Lognormal(parameter_name, parameter_mean, parameter_lower_95_interval, ...
parameter_upper_95_interval)
    
    variance = 0:0.00001:2;

    lower_intervals = logninv(0.025, log(parameter_mean - variance), sqrt(variance));

    best_variance = variance(find(abs(lower_intervals - parameter_lower_95_interval) == ...
    min(abs(lower_intervals - parameter_lower_95_interval))));

    best_mu = log(parameter_mean - best_variance);

    error_upper_95 = abs(parameter_upper_95_interval - logninv(0.975, best_mu, sqrt(best_variance)));

    per_error_upper_95 = 100 * error_upper_95 / parameter_upper_95_interval;

    fprintf("---- Lognormal fitting results for parameter %s ----\n", parameter_name);
    fprintf("mu: %.2f\n", best_mu);
    fprintf("sigma^2: %.7f\n", best_variance);
    fprintf("Absolute error in upper confidence interval: %.7f\n", error_upper_95);
    fprintf("Percentage error in upper confidence interval: %.7f\n\n", per_error_upper_95);

    if(best_variance == variance(end))
        fprintf("WARNING: Variance trial values too small. Increase and try again.\n\n");
    end

    param_info = {parameter_name, "lognormal", best_mu, sqrt(best_variance)};
end