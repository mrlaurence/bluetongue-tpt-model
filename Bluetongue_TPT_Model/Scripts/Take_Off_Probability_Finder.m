close all;

for i = 1 : 200

    fprintf("[INFO] Running simulation %d.\n", i);

    generatedValidParams = false;

    while ~generatedValidParams
        % Construct a structure containing within-farm transmission parameters.
        farmSquareParameters = Bluetongue_TPT_Model.Define_Within_Farm_Parameters();

        if min(cell2mat(struct2cell(rmfield(farmSquareParameters, {'b_11', 'b_12', 'b_21', 'b_22'})))) >= 0
            generatedValidParams = true;
        else
            fprintf("[WARN] Invalid within-farm parameters drawn. Redrawing...\n");
        end
    end

    farm = Bluetongue_TPT_Model.Initialise_Farm_Square_Model(1,1,farmSquareParameters,140,[1500,2000,400]);

    infecVecCount = round(0.000005 * farm.vectorStates.tensor(1,1,1));
    
    farm.vectorStates.Set_Val("I", 1, infecVecCount);

    for j = 1 : 8000
        farm.Tau_Leap_Step([1 1 1 1 1 1], 0.05, 23);
    end

    figure("Visible","off");
    hold on;
    susceptible = farm.dairyCattleStates.Get_Time_Series("X_", true); infectious = farm.dairyCattleStates.Get_Time_Series("Y_", true); recovered = farm.dairyCattleStates.Get_Time_Series("Z_", true); total = susceptible(:,2) + infectious(:,2) + recovered(:,2); total = [susceptible(:,1) total]; tsplot(susceptible,0.05);tsplot(infectious,0.05);tsplot(recovered,0.05);tsplot(total,0.05);legend(["S","I","R","Total"]);
    saveas(gcf, sprintf("Figures/%d.png", i));
    close;
end