function param_info = Fit_Normal(parameter_name, parameter_mean, parameter_lower_95_interval, ...
parameter_upper_95_interval)
    
    best_mu = parameter_mean;

    variance = 0:0.01:200000;

    lower_intervals = norminv(0.025, best_mu, sqrt(variance));

    best_variance = variance(find(abs(lower_intervals - parameter_lower_95_interval) == ...
    min(abs(lower_intervals - parameter_lower_95_interval))));

    error_upper_95 = abs(parameter_upper_95_interval - norminv(0.975, best_mu, sqrt(best_variance)));

    per_error_upper_95 = 100 * error_upper_95 / abs(parameter_upper_95_interval);

    fprintf("---- Normal fitting results for parameter %s ----\n", parameter_name);
    fprintf("mu: %.2f\n", best_mu);
    fprintf("sigma^2: %.7f\n", best_variance);
    fprintf("Absolute error in upper confidence interval: %.7f\n", error_upper_95);
    fprintf("Percentage error in upper confidence interval: %.7f\n\n", per_error_upper_95);

    if(best_variance == variance(end))
        fprintf("WARNING: Variance trial values too small. Increase and try again.\n\n");
    end

    param_info = {parameter_name, "normal", best_mu, sqrt(best_variance)};
end