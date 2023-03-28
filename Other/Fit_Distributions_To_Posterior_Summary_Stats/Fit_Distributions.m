posterior_matrix = ["b",      "0.84",   "0.68",    "0.96",    "beta";
                    "beta",   "0.023",  "0.0073",  "0.042",   "beta";
                    "mu_V",   "1806.7", "688.2",   "3141.1",  "gamma";
                    "s_V",    "1.69",   "0.54",    "3.17",    "gamma";
                    "sigma",  "0.095",  "0.0029",  "0.40",    "exponential";
                    "1/r_C",  "20.53",  "18.76",   "22.26",   "normal";
                    "n_C",    "5.3",    "4",       "6",       "normal";
                    "d_C",    "0.0012", "0.00012", "0.0029",  "normal";
                    "1/r_S",  "16.10",  "14.12",   "18.20",   "normal";
                    "n_S",    "12.4",   "6",       "20",      "gamma";
                    "d_S",    "0.0068", "0.00059", "0.017",   "normal";
                    "alpha",  "0.020",  "0.016",   "0.024",   "normal";
                    "T_min",  "13.2",   "12.8",    "13.7",    "normal";
                    "k",      "11",     "3",       "22",      "gamma";
                    "b_11",  "-1.59",  "-1.80",   "-1.38",    "normal";
                    "b_21",  "-3.80",  "-4.37",   "-3.22",    "normal";
                    "b_12",  "-1.46",  "-1.60",   "-1.32",    "normal";
                    "b_22",  "-0.99",  "-1.38",   "-0.60",    "normal";
                   ];

param_info_matrix = {};

posterior_matrix_size = size(posterior_matrix);

param_count = posterior_matrix_size(1);

for i = 1 : param_count
    param_posterior = posterior_matrix(i,:);

    switch param_posterior(5)
        case "beta"
            param_info = Fit_Beta(param_posterior(1),str2double(param_posterior(2)),str2double(param_posterior(3)),str2double(param_posterior(4)));
        case "gamma"
            param_info = Fit_Gamma(param_posterior(1),str2double(param_posterior(2)),str2double(param_posterior(3)),str2double(param_posterior(4)));
        case "exponential"
            param_info = Fit_Exponential(param_posterior(1),str2double(param_posterior(2)),str2double(param_posterior(3)),str2double(param_posterior(4)));
        case "normal"
            param_info = Fit_Normal(param_posterior(1),str2double(param_posterior(2)),str2double(param_posterior(3)),str2double(param_posterior(4)));
        case "lognormal"
            param_info = Fit_Lognormal(param_posterior(1),str2double(param_posterior(2)),str2double(param_posterior(3)),str2double(param_posterior(4)));
    end
    
    param_info_matrix(end + 1, :) = param_info;
end

param_info_matrix