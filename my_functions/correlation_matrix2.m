function correlation_matrix2(eventMatrix, eventLabels, criteriaLabels, colorlimitsArray)
    nRows = size(eventMatrix, 1);
    nCols = size(eventMatrix, 2);
    nColors = 128;
    
    % Create colormap: red-white-blue
    nHalf = nColors / 2;
    redToWhite = [linspace(1, 1, nHalf)', linspace(0, 1, nHalf)', linspace(0, 1, nHalf)'];
    whiteToBlue = [linspace(1, 0, nHalf)', linspace(1, 1, nHalf)', linspace(1, 1, nHalf)'];
    customColormap = [redToWhite; whiteToBlue];
    
    % Normalize each subset separately
    colSubsetSize = numel(colorlimitsArray);
    colsPerSubset = ceil(nCols / colSubsetSize);
    normMatrix = zeros(size(eventMatrix));
    
    for subsetIdx = 1:colSubsetSize
        colStart = (subsetIdx - 1) * colsPerSubset + 1;
        colEnd = min(subsetIdx * colsPerSubset, nCols);
        
        % Normalize the matrix within the subset
        subsetData = eventMatrix(:, colStart:colEnd);
        limits = colorlimitsArray{subsetIdx};
        normMatrix(:, colStart:colEnd) = (subsetData - limits(1)) / (limits(2) - limits(1));
        normMatrix(:, colStart:colEnd) = max(0, min(1, normMatrix(:, colStart:colEnd))); % Clamp between 0 and 1
    end
    
    % Plot heatmap using imagesc
    figure;
    imagesc(normMatrix); 
    colormap(customColormap);
    colorbar;
    
    % Set labels
    set(gca, 'YTick', 1:nRows, 'YTickLabel', criteriaLabels, 'TickLabelInterpreter', 'none');
    set(gca, 'XTick', 1:nCols, 'XTickLabel', eventLabels, 'TickLabelInterpreter', 'none');
    
    % Annotate with values
    for i = 1:nRows
        for j = 1:nCols
            value = eventMatrix(i, j);
            text(j, i, sprintf('%.2f', value), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', 'FontSize', 8);
        end
    end
    
    title('Combined Heatmap with Different Color Limits');
end
