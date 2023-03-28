function tsplot(timeSeries, tau, noLabel, yearsScale)

    % Plots a time series matrix with an appropriate scale for the x-axis. It is
    % recommended to run Bluetongue_TPT_Model.Configure_Figure_Presentation for clear
    % display of MATLAB figures.
    %
    % timeSeries: A two-column matrix, where the first column specifies 'unscaled' time
    % values and the second specifies the value of the time series at each time value.
    % tau (optional): A scale factor used to compute the 'scaled' time values which are
    % plotted. If s,t are unscaled/scaled time values respectively, we have t = tau * (s -
    % 1).
    % noLabel (optional): A bool indicating whether a placeholder title and basic labels
    % for the axes should be added to the plot.
    % yearsScale (optional): A bool indicating whether scaled time values should be
    % further divided by 365 before plotting.
    %
    % AUTHOR: Laurence Dhonau.

    % Set default values for optional arguments.
    if nargin == 1
        tau = 1;
        noLabel = false;
        yearsScale = false;
    end

    if nargin == 2
        noLabel = false;
        yearsScale = false;
    end

    if nargin == 3
        yearsScale = false;
    end

    if ~yearsScale
        plot(tau * (timeSeries(:,1) - 1), timeSeries(:,2));
    else
        plot(tau * (timeSeries(:,1) - 1) / 365, timeSeries(:,2));
    end

    if ~noLabel

        if yearsScale
            xlabel("Time after BT introduced / years");
        else
            xlabel("Time after BT introduced / days");
        end
        
        ylabel("Dummy label");
        title("Dummy title");
    end
end