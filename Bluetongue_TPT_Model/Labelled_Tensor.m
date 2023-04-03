classdef Labelled_Tensor < handle

    % A 'wrapper' class for MATLAB's multidimensional arrays, which we call 'tensors'. A
    % Labelled_Tensor object encapsulates a tensor as the evolution of a matrix over the
    % third-dimension index and allows each matrix position to be labelled with a string
    % value. Elements in the tensor can then be accessed by their label and
    % third-dimensional 'time' index.

    % AUTHOR: Laurence Dhonau.

    properties
        tensor
        tensorLabelMatrix
    end

    methods
        function self = Labelled_Tensor(tensor, tensorLabelMatrix)

            % Constructs a Labelled_Tensor object with the passed values.
            %
            % tensor: A three-dimensional array.
            % tensorLabelMatrix: A string array matrix of the same size as the first two
            % dimensions of tensor.

            self.tensor = tensor;
            self.tensorLabelMatrix = tensorLabelMatrix;
        end

        function value = Get_Val(self, label, index)

            % Retrieves an element in the tensor given the label matching its position in
            % the matrix, and its third-dimension index.
            %
            % label: A string matching exactly one element in tensorLabelMatrix.
            % index: A positive integer (at most the size of the third-dimension of
            % tensor).
            %
            % Returns the retrieved value.

            % Find the position of the passed label in the label matrix.
            [labelRow, labelColumn] = find(self.tensorLabelMatrix == label);

            % Raise an error if the label does not exist or there are multiple matching
            % labels in the label matrix.
            if length(labelRow) + length(labelColumn) ~= 2
                error("Invalid lookup attempted on a Labelled_Tensor. Label does not exist " + ...
                      "or label matrix is invalid (duplicate entries).");
            end

            % Return the requested value.
            value = self.tensor(labelRow, labelColumn, index);
        end

        function Set_Val(self, label, index, new_value)

            % Updates an element in the tensor given the label matching its position in
            % the matrix, and its third-dimension index.
            %
            % label: A string matching exactly one element in tensorLabelMatrix.
            % index: A positive integer (at most the size of the third-dimension of
            % tensor).
            % new_value: The value to set the tensor element to.

            % Find the position of the passed label in the label matrix.
            [labelRow, labelColumn] = find(self.tensorLabelMatrix == label);

            % Raise an error if the label does not exist or there are multiple matching
            % labels in the label matrix.
            if length(labelRow) + length(labelColumn) ~= 2
                error("Invalid lookup attempted on a Labelled_Tensor. Label does not exist " + ...
                      "or label matrix is invalid (duplicate entries).");
            end

            % Update the value in the tensor.
            self.tensor(labelRow, labelColumn, index) = new_value;
        end


        function Set_By_Query_String(self, queryString, value, index, silent)

            % Updates multiple elements with a shared third-dimensional index in the
            % tensor based on a 'query string' matching the labels of each element.
            %
            % queryString: A string *contained by* the elements of tensorLabelMatrix which
            % occupy the positions in tensor to be updated.
            % value: Either a single value (if all elements are to be set to the same
            % value) or a vector (which must be the same size as the number of elements
            % matched by the query string)
            % index (optional): A positive integer (at most the size of the third-dimension
            % of tensor) in which to update elements of the tensor.
            % silent (optional): A bool indicating whether output informing the user which
            % labels in tensorLabelMatrix were matched should be surpressed.

            % Set default values for optional arguments.
            if nargin == 3
                index = 1;
                silent = false;
            end
            if nargin == 4
                silent = false;
            end

            % Find the position of the labels in the label matrix which contain the passed
            % query string.
            [labelRows, labelColumns] = find(contains(self.tensorLabelMatrix, queryString));

            % If no labels were matched.
            if size(labelRows, 1) == 0
                error("Invalid 'set by query string' operation attempted on a Labelled_Tensor." + ...
                " No labels match query string.");

            % If the number of values passed is not 1 but does not agree with the number
            % of labels matched.
            elseif size(value, 2) ~= 1 && size(value, 2) ~= size(labelRows, 1)
                error("Invalid 'set by query string' operation attempted on a Labelled_Tensor." + ...
                " Value vector is not of valid length.");

            elseif ~silent
                fprintf("Matched labels: ");
            end

            if size(labelRows, 1) == 1
                labelPosition = [labelRows' labelColumns'];
            else
                labelPosition = [labelRows labelColumns];
            end

            % Iterate over the matched labels.
            for i = 1 : max(size(labelRows))

                % If not in silent mode, output the matched label (and a comma if this is
                % not the final label).
                if ~silent
                    fprintf(self.tensorLabelMatrix(labelPosition(i,1), labelPosition(i,2)));

                     if i < max(size(labelRows))
                         fprintf(", ");
                     end
                end

                % Update the element in tensor which matches this label, either setting it
                % to value (if value is a scalar) or the ith element of value (if value is
                % a vector).
                if size(value, 2) == 1
                    self.Set_Val(self.tensorLabelMatrix(labelPosition(i,1), ...
                    labelPosition(i,2)), index, value);
                else
                    self.Set_Val(self.tensorLabelMatrix(labelPosition(i,1), ...
                    labelPosition(i,2)), index, value(i));
                end

            end
            if ~silent
                fprintf("\n");
            end
        end

        function timeSeries = Get_Time_Series(self, queryString, silent, timeIndices)

            % Compute a time series spanning the passed third-dimensional indices by
            % summing up values in the tensor, selected based on a 'query string' matching
            % the labels of each element.
            %
            % queryString: A string *contained by* the elements of tensorLabelMatrix which
            % occupy the positions in tensor to be included in the time series.
            % silent (optional): A bool indicating whether output informing the user which
            % labels in tensorLabelMatrix were matched should be surpressed.
            % timeIndices (optional): A vector of positive integers (each at most the size
            % of the third-dimension of tensor) for which a sum of the elements in the
            % positions matched by queryString should be summed and included in the time
            % series.
            %
            % Returns a vector giving the value of the sum at each third-dimensional index
            % specified in timeIndices.

            % Set default values for optional arguments.
            if nargin == 2
                timeIndices = 1:size(self.tensor, 3);
                silent = false;
            end

            if nargin == 3
                timeIndices = 1:size(self.tensor, 3);
            end

            % Initalise a vector to store the time series.
            timeSeries = zeros(size(timeIndices));

            % Find the position of the labels in the label matrix which contain the passed
            % query string.
            [labelRows, labelColumns] = find(contains(self.tensorLabelMatrix, queryString));

            % Raise an error if no labels were matched.
            if size(labelRows, 1) == 0
                error("Invalid time series construction attempted on a Labelled_Tensor." + ...
                " No labels match query string.");
            elseif ~silent
                fprintf("Matched labels: ");
            end

            if size(labelRows, 1) == 1
                labelPosition = [labelRows' labelColumns'];
            else
                labelPosition = [labelRows labelColumns];
            end

            % Iterate over the time indices to be included in the time series.
            for i = 1 : max(size(timeIndices))

                % Retrieve the third-dimensional index.
                index = timeIndices(i);

                % Iterate over the matched labels.
                for j = 1 : max(size(labelRows))

                    % If this is the first time index and not in silent mode, output the
                    % matched label (and a comma if this is not the final label).
                    if i == 1 && ~silent
                        fprintf(self.tensorLabelMatrix(labelPosition(j,1), labelPosition(j,2)));

                        if j < max(size(labelRows))
                            fprintf(", ");
                        end
                    end

                    % Add the value in the tensor matching the label position with
                    % third-dimensional index given by index to the value in the time
                    % series for this index.
                    timeSeries(i) = timeSeries(i) + self.tensor(labelPosition(j,1), ...
                    labelPosition(j,2), index);
                end
                
            end

            if ~silent
                fprintf("\n");
            end
        end
    end
end