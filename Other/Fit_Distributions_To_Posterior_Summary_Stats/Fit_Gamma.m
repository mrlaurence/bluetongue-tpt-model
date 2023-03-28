function param_info = Fit_Gamma(parameter_name, parameter_mean, parameter_lower_95_interval, ...
parameter_upper_95_interval)
    
    theta = 0.01:0.5:10000;

    lower_intervals = gaminv(0.025, parameter_mean ./ theta, theta);

    best_theta = theta(find(abs(lower_intervals - parameter_lower_95_interval) == ...
    min(abs(lower_intervals - parameter_lower_95_interval))));

    best_k = parameter_mean / best_theta;

    error_upper_95 = abs(parameter_upper_95_interval - gaminv(0.975, best_k, best_theta));

    per_error_upper_95 = 100 * error_upper_95 / parameter_upper_95_interval;

    fprintf("---- Gamma fitting results for parameter %s ----\n", parameter_name);
    fprintf("theta: %.2f\n", best_theta);
    fprintf("k: %.2f\n", best_k);
    fprintf("Absolute error in upper confidence interval: %.4f\n", error_upper_95);
    fprintf("Percentage error in upper confidence interval: %.2f\n\n", per_error_upper_95);

    if(best_theta == theta(end))
        fprintf("WARNING: Theta trial values too small. Increase and try again.\n\n");
    end

    param_info = {parameter_name, "gamma", best_k, best_theta};
end