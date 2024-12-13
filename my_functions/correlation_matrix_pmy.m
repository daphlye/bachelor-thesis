function correlation_matrix_pmy(colorMatrix, eventMatrix, eventLabels, criteriaLabels, colorlimitsArray)
% correlation_matrix_pmy can plot a hatmapped table that has a heatmap
% based on a different array than the
% Input:
% colorMatrix: Used for coloring
% eventMatrix: Used for displaying values, has to be the same size as
% colorMatrix
% eventLabels: X-axis labels
% criteriaLabels: Y-axis labels
% colorlimitsArray: Color limits
    
    if isequal(size(colorMatrix), size(eventMatrix))
        disp('The matrices have the same size.');
    else
        disp('The matrices have different sizes.');
    end

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
        
        subsetData = table2array(colorMatrix(:, colStart:colEnd));
        limits = colorlimitsArray{subsetIdx};
        normMatrix(:, colStart:colEnd) = (subsetData - limits(1)) / (limits(2) - limits(1));
        normMatrix(:, colStart:colEnd) = max(0, min(1, normMatrix(:, colStart:colEnd)));
    end
    
    %% Plot the heatmap
    figure;
    imagesc(normMatrix);
    colormap(customColormap);
    clim([0, 1]); % Set color axis limits between 0 and 1
    cb = colorbar;
    cb.Label.String = 'correlation';
    cb.Label.FontSize = 15;
    
    % Set labels
    set(gca, 'YTick', 1:nRows, 'YTickLabel', eventLabels, 'TickLabelInterpreter', 'latex', 'FontAngle','normal','FontSize',15);
    set(gca, 'XTick', 1:nCols, 'XTickLabel', criteriaLabels, 'TickLabelInterpreter', 'latex', 'FontAngle','normal','FontSize',15);
    
    % Display eventMatrix on the heatmap
    for i = 1:nRows
        for j = 1:nCols
            text(j, i, sprintf('%.f', eventMatrix(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', 'FontSize', 15);
        end
    end
    
    title('$\mathbf{J}$ in predicted direction', 'Interpreter','latex');
end

