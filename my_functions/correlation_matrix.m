function correlation_matrix_color_and_value(colorMatrix, eventMatrix, eventLabels, criteriaLabels, colorlimitsArray)
    % Input:
    % colorMatrix: Used for coloring
    % eventMatrix: Used for displaying values
    % eventLabels: X-axis labels
    % criteriaLabels: Y-axis labels
    % colorlimitsArray: Color limits
    
    % Create a custom colormap (red-white-blue)
    nColors = 128;
    nHalf = nColors / 2;
    redToWhite = [linspace(1, 1, nHalf)', linspace(0, 1, nHalf)', linspace(0, 1, nHalf)'];
    whiteToBlue = [linspace(1, 0, nHalf)', linspace(1, 1, nHalf)', linspace(1, 1, nHalf)'];
    customColormap = [redToWhite; whiteToBlue];

    % Normalize colorMatrix
    nRows = size(colorMatrix, 1);
    nCols = size(colorMatrix, 2);
    colSubsetSize = numel(colorlimitsArray);
    colsPerSubset = ceil(nCols / colSubsetSize);
    normMatrix = zeros(size(colorMatrix));
    
    for subsetIdx = 1:colSubsetSize
        colStart = (subsetIdx - 1) * colsPerSubset + 1;
        colEnd = min(subsetIdx * colsPerSubset, nCols);
        
        subsetData = colorMatrix(:, colStart:colEnd);
        limits = colorlimitsArray{subsetIdx};
        normMatrix(:, colStart:colEnd) = (subsetData - limits(1)) / (limits(2) - limits(1));
        normMatrix(:, colStart:colEnd) = max(0, min(1, normMatrix(:, colStart:colEnd)));
    end
    
    % Plot the heatmap
    figure;
    imagesc(normMatrix);
    colormap(customColormap);
    clim([0, 1]); % Set color axis limits between 0 and 1
    colorbar;
    
    % Set labels
    set(gca, 'YTick', 1:nRows, 'YTickLabel', eventLabels, 'TickLabelInterpreter', 'none');
    set(gca, 'XTick', 1:nCols, 'XTickLabel', criteriaLabels, 'TickLabelInterpreter', 'none');
    
    % Display eventMatrix on the heatmap
    for i = 1:nRows
        for j = 1:nCols
            text(j, i, sprintf('%.2f', eventMatrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', 'FontSize', 8);
        end
    end
    
    title('Heatmap with Values from Different Matrix');
end

