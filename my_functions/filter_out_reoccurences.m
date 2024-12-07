function [filtered_matrix] = filter_out_reoccurences(matrix, target_column)
% Filters rows of a matrix based on unique values in a specific column.
% 
% Inputs:
%   matrix: The input matrix to filter.
%   target_column: The column index to use for filtering (e.g., 1 for the first column).
% Output:
%   filtered_matrix: The filtered matrix with rows corresponding to unique values.

    % Find unique values and their first occurrences in the target column
    [~, first_indices] = unique(matrix(:, target_column), 'stable');

    % Extract rows based on the first occurrences
    filtered_matrix = matrix(first_indices, :);
end
