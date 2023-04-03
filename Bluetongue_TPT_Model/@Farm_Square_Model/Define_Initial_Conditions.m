function Define_Initial_Conditions(self)

    % Prepares the Farm_Square_Model object for simulation by setting the
    % dairyCattleStates, beefCattleStates, sheepStates and vectorStates Labelled_Tensor
    % properties based on the value of self.initialPopulations.
    %
    % It is assumed that the 4 Labelled_Tensors have each been constructed with a zero
    % matrix of the correct size (spanning just one time value). It is also assumed that
    % initialPopulations is a struct with three fields: "dairy", "beef", "sheep". The
    % value of each field should be a table with three columns, where each row specifies a
    % batch of values to be set. In each row, the first value specifies the query string
    % used to find the model compartments for which an initial value is to be set; the
    % second is a cell containing a vector and the third is a cell containing a boolean.
    % There are three valid cases for a row:
    % 
    % CASE 1 (vector with >1 element):
    % e.g. "X", {[3, 4, 5, 6, 7]}, {false}
    % The bool must be false. The number of elements in the vector must equal the number
    % of compartments which "X" matches to when the Labelled_Tensor is searched by query
    % string (see Labelled_Tensor.Set_By_Query_String()). For example, the Labelled_Tensor
    % could contain compartments "X_Y1.", "X_Y2.", "X_A.", "X_P1." and "X_V.". In this case,
    % the value of each of these compartments is set to the corresponding value in the
    % vector.
    %
    % CASE 2 (vector with 1 element, splitEqually = false):
    % e.g. "X", {[3]}, {false}
    % In this case, each of the compartments "X" matches to when the Labelled_Tensor is
    % searched by query string are set to the value in the vector (3).
    %
    % CASE 3 (vector with 1 element, splitEqually = true):
    % e.g. "X", {[25]}, {true}
    % In this case, each of the compartments "X" matches to when the Labelled_Tensor is
    % searched by query string are set to the value in the vector divided by the number of
    % such compartments (rounded to the nearest integer).
    %
    % The vectorStates Labelled_Tensor is also set to have an entirely suspectible vector
    % population. The size of this population is based on the total cattle and sheep
    % populations specified in self.initialPopulations, as well as the vector-to-host
    % ratios for cattle/sheep (found in self.parameters).
    %
    % AUTHOR: Laurence Dhonau.

    % Iterate over the rows in the dairy table.
    for i = 1:height(self.initialPopulations.dairy)

        % Retrieve the query string, vector of values and splitEqually bool for this row.
        queryString = table2array(self.initialPopulations.dairy(i,1));
        classVector = cell2mat(table2array(self.initialPopulations.dairy(i,2)));
        splitEqually = cell2mat(table2array(self.initialPopulations.dairy(i,3)));

        if splitEqually
            
            % Find the number of dairy compartments matching the query string.
            [labelRows, ~] = find(contains(self.dairyCattleStates.tensorLabelMatrix, ...
            queryString));
            compartmentCount = size(labelRows, 1);

            % We will set each compartment to the requested value divided by the
            % compartment count.
            valueToSet = round(classVector / compartmentCount);
        else

            % We will set each compartment to the requested value.
            valueToSet = classVector;
        end

        % Set each of the matched compartments to valueToSet without notifying the user of
        % which compartments were matched.
        self.dairyCattleStates.Set_By_Query_String(queryString, valueToSet, 1, true);
    end

    % Repeat the above process for the beef cattle population (the process is identical).
    for i = 1:height(self.initialPopulations.beef)

        queryString = table2array(self.initialPopulations.beef(i,1));
        classVector = cell2mat(table2array(self.initialPopulations.beef(i,2)));
        splitEqually = cell2mat(table2array(self.initialPopulations.beef(i,3)));

        if splitEqually

            [labelRows, ~] = find(contains(self.beefCattleStates.tensorLabelMatrix, ...
            queryString));
            compartmentCount = size(labelRows, 1);

            valueToSet = round(classVector / compartmentCount);
        else
            valueToSet = classVector;
        end

        self.beefCattleStates.Set_By_Query_String(queryString, valueToSet, 1, true);
    end

    % Repeat the above process for the sheep population (the process is identical).
    for i = 1:height(self.initialPopulations.sheep)

        queryString = table2array(self.initialPopulations.sheep(i,1));
        classVector = cell2mat(table2array(self.initialPopulations.sheep(i,2)));
        splitEqually = cell2mat(table2array(self.initialPopulations.sheep(i,3)));

        if splitEqually

            [labelRows, ~] = find(contains(self.sheepStates.tensorLabelMatrix, ...
            queryString));
            compartmentCount = size(labelRows, 1);

            valueToSet = round(classVector / compartmentCount);
        else
            valueToSet = classVector;
        end

        self.sheepStates.Set_By_Query_String(queryString, valueToSet, 1, true);
    end

    % Find the total number of dairy cattle, beef and sheep in the farm square.
    totalDairyCattle = sum(sum(self.dairyCattleStates.tensor(:,:,1)));
    totalBeefCattle = sum(sum(self.beefCattleStates.tensor(:,:,1)));
    totalSheep = sum(self.sheepStates.tensor(:,:,1));

    % Set the suspectible vector compartment based on the total populations and
    % vector-to-host ratios.
    self.vectorStates.Set_Val("S.", 1, round(self.parameters.m_C * (totalDairyCattle + ...
    totalBeefCattle) + self.parameters.m_S * totalSheep));
end