function cattleLabelMatrix = Build_Cattle_Label_Matrix(classCounts, classLabels)

    % Constructs a string array which contains labels for all the compartments for either
    % dairy cattle or beef cattle in the farm square.
    %
    % classCounts: A 3-element vector, where the elements are, in order: the number of
    % youth stages, the number of pregnancy stages, the number of viraemia stages.
    %
    % classLabels: A 4-element string array, where the elements are, in order: the
    % character(s) which denote the youth, adult, pregnancy and postpartum classes.
    % An example would be ["Y", "A", "P", "V"].
    %
    % Returns the constructed string array. For example, if classCounts = [2, 3, 2] and
    % classLabels = ["A", "B", "C", "D"] then the return value is:
    % ["X_A1."    "Y_1,A1."    "Y_2,A1."    "Z_A1.";
    %  "X_A2."    "Y_1,A2."    "Y_2,A2."    "Z_A2.";
    %  "X_B."     "Y_1,B."     "Y_2,B."     "Z_B." ;
    %  "X_C1."    "Y_1,C1."    "Y_2,C1."    "Z_C1.";
    %  "X_C2."    "Y_1,C2."    "Y_2,C2."    "Z_C2.";
    %  "X_C3."    "Y_1,C3."    "Y_2,C3."    "Z_C3.";
    %  "X_D."     "Y_1,D."     "Y_2,D."     "Z_D." ].
    %
    % AUTHOR: Laurence Dhonau.

    cattleLabelMatrix = strings(classCounts(1) + classCounts(2) + 2, classCounts(3) + 2);

    ageLabels = strings(classCounts(1) + classCounts(2) + 2, 1);
    ageLabels(classCounts(1) + 1) = classLabels(2);
    ageLabels(classCounts(1) + classCounts(2) + 2) = classLabels(4);

    for i = 1 : classCounts(1)
        ageLabels(i) = strcat(classLabels(1), num2str(i));
    end

    for i = 1 : classCounts(2)
        ageLabels(classCounts(1) + 1 + i) = strcat(classLabels(3), num2str(i));
    end

    cattleLabelMatrix(:, 1) = vecsprintf("X_%s.", ageLabels);
    cattleLabelMatrix(:, end) = vecsprintf("Z_%s.", ageLabels);

    for i = 1 : classCounts(3)
        cattleLabelMatrix(:, i+1) = vecsprintf(strcat(sprintf("Y_%s,", num2str(i)), "%s."), ageLabels);
    end
end