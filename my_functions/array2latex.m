function latex_code = array2latex(input)
    % Convert MATLAB arrays, symbolic arrays, or tables into a LaTeX matrix or table format.
%--------------------------------------------------------------------------
    
    if isa(input, 'table')
        % Handle MATLAB tables
        latex_code = '\begin{table}[h]\centering\begin{tabular}{';
        
        % Generate column alignment (one column per variable + row names if they exist)
        hasRowNames = ~isempty(input.Properties.RowNames);
        if hasRowNames
            latex_code = [latex_code, 'l'];  % Left-align for row names
        end
        latex_code = [latex_code, repmat('c', 1, width(input)), '}'];
        
        % Add header row with variable names
        if hasRowNames
            latex_code = [latex_code, '\hline ', '& '];
        end
        latex_code = [latex_code, strjoin(input.Properties.VariableNames, ' & '), ' \\ \hline'];
        
        % Add table rows
        for r = 1:height(input)
            row_str = '';
            
            % Add row names if they exist
            if hasRowNames
                row_str = [row_str, input.Properties.RowNames{r}, ' & '];
            end
            
            % Add data for the row
            row_data = input(r, :).Variables;
            row_str = [row_str, strjoin(arrayfun(@(x) num2str(x), row_data, 'UniformOutput', false), ' & ')];
            
            % Add to LaTeX code
            latex_code = [latex_code, row_str, ' \\'];
        end
        
        % Close the LaTeX tabular environment
        latex_code = [latex_code, ' \hline \end{tabular}\end{table}'];
    
    elseif isa(input, 'sym') || isnumeric(input)
        % Handle arrays (numeric or symbolic)
        if isa(input, 'sym')
            is_sym = true;  % Mark as symbolic
        else
            is_sym = false;
        end
        
        [rows, cols] = size(input);
        latex_code = '\begin{align*}\begin{pmatrix}';  % Start the LaTeX matrix
        
        % Loop through the rows
        for r = 1:rows
            row_str = '';  % Initialize row string
            
            % Loop through the columns
            for c = 1:cols
                if is_sym
                    % If symbolic, use char to handle symbolic expression
                    row_str = [row_str, char(input(r, c))];
                else
                    % If numeric, just use num2str
                    row_str = [row_str, num2str(input(r, c))];
                end
                
                % Add space or separator if not the last element
                if c < cols
                    row_str = [row_str, ' & '];  % LaTeX row separator
                end
            end
            
            % Add the current row to the LaTeX code
            latex_code = [latex_code, row_str, ' \\'];  % End the row with LaTeX line break
        end
        
        % Finish the LaTeX matrix code
        latex_code = [latex_code, '\end{pmatrix}\end{align*}'];
    
    else
        error('Input must be a table, numeric array, or symbolic array.');
    end
end


