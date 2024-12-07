function latex_code = table2latex(tbl)
% Convert MATLAB table to LaTeX table format
% tbl: MATLAB table with mixed data types

    % Initialize LaTeX table code
    latex_code = '\begin{table}[h]\centering\n\begin{tabular}{';
    col_format = repmat('c', 1, width(tbl)); % Center all columns
    latex_code = [latex_code, col_format, '}\n']; % Column formatting
    
    % Add column headers
    header_row = strjoin(tbl.Properties.VariableNames, ' & ');
    latex_code = [latex_code, header_row, ' \\ \hline\n'];
    
    % Iterate through rows of the table
    for r = 1:height(tbl)
        row_str = ''; % Initialize row string
        for c = 1:width(tbl)
            cell_val = tbl{r, c};
            if isdatetime(cell_val) % Format datetime to string
                cell_str = char(cell_val);
            else % For numeric or other types, convert to string
                cell_str = num2str(cell_val, '%.2f'); % Format to 2 decimal places
            end
            row_str = [row_str, cell_str]; % Append to row string
            if c < width(tbl) % Add column separator if not last column
                row_str = [row_str, ' & '];
            end
        end
        % Add current row to LaTeX table
        latex_code = [latex_code, row_str, ' \\ \n'];
    end
    
    % Finish the LaTeX table code
    latex_code = [latex_code, '\hline\n\end{tabular}\n\end{table}'];
end


