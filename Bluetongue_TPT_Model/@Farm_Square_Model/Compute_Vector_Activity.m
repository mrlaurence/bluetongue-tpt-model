function activity = Compute_Vector_Activity(par, t, activityNormaliser)
    
    % Computes the seasonal vector activity value on the farm square on a given day.
    %
    % par: The within-farm BTV parameters struct. Bluetongue_TPT_Model.Define_Within_Farm_
    % Parameters() generates a structure of the required form.
    % t: The current simulation time.
    % activityNormaliser: A value to multiply the "raw" vector activity value by. This may
    % be used, for example, to ensure that the maximum possible vector activity value is 1.
    %
    % Returns the seasonal vector activity at time t.
    %
    % AUTHOR: Laurence Dhonau.

    activity = activityNormaliser * exp(par.b_11 * sin(2*pi*t/365) + par.b_21 * cos(2*pi*t/365) + ...
    par.b_12 * sin(4*pi*t/365) + par.b_22 * cos(4*pi*t/365));
end