function param_info = Fit_Exponential(parameter_name, parameter_mean, parameter_lower_95_interval, ...
parameter_upper_95_interval)
    
    best_lambda = parameter_mean;

    error_lower_95 = abs(parameter_lower_95_interval - expinv(0.025, best_lambda));
    per_error_lower_95 = 100 * error_lower_95 / parameter_lower_95_interval;
    error_upper_95 = abs(parameter_upper_95_interval - expinv(0.975, best_lambda));
    per_error_upper_95 = 100 * error_upper_95 / parameter_upper_95_interval;

    fprintf("---- Exponential fitting results for parameter %s ----\n", parameter_name);
    fprintf("Lambda: %.4f\n", best_lambda);
    fprintf("Absolute error in lower confidence interval: %.4f\n", error_lower_95);
    fprintf("Percentage error in lower confidence interval: %.2f\n\n", per_error_lower_95);
    fprintf("Absolute error in upper confidence interval: %.4f\n", error_upper_95);
    fprintf("Percentage error in upper confidence interval: %.2f\n\n", per_error_upper_95);

    param_info = {parameter_name, "exponential", best_lambda, NaN};
end