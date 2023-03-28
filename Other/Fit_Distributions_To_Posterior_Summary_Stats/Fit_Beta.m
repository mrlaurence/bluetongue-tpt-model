function param_info = Fit_Beta(parameter_name, parameter_mean, parameter_lower_95_interval, ...
parameter_upper_95_interval)
    
    beta = 0.01:0.01:1000;

    lower_intervals = betainv(0.025, parameter_mean / (1 - parameter_mean) * beta, beta);

    best_beta = beta(find(abs(lower_intervals - parameter_lower_95_interval) == ...
    min(abs(lower_intervals - parameter_lower_95_interval))));

    best_alpha = parameter_mean / (1 - parameter_mean) * best_beta;

    error_upper_95 = abs(parameter_upper_95_interval - betainv(0.975, best_alpha, best_beta));

    per_error_upper_95 = 100 * error_upper_95 / parameter_upper_95_interval;

    fprintf("---- Beta fitting results for parameter %s ----\n", parameter_name);
    fprintf("Alpha: %.2f\n", best_alpha);
    fprintf("Beta: %.2f\n", best_beta);
    fprintf("Absolute error in upper confidence interval: %.4f\n", error_upper_95);
    fprintf("Percentage error in upper confidence interval: %.2f\n\n", per_error_upper_95);

    if(best_beta == beta(end))
        fprintf("WARNING: Beta trial values too small. Increase and try again.\n\n");
    end

    param_info = {parameter_name, "beta", best_alpha, best_beta};
end