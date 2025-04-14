function to_latex_convergence_table(nstep, A, filename)
% Open file for writing
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file for writing.');
end

% Write LaTeX table header
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{c|%s}\n', repmat('c', 1, length(nstep)*2-1));  % Centered columns
fprintf(fid, '\\hline\n');
%write header
fprintf(fid, ['&' ]);
for i=1:length(nstep)-1
    fprintf(fid, ['$E\\left( \\frac{1}{',num2str(nstep(i)),'}\\right)$ &  $\\mathcal{O}$   &' ]);
end

fprintf(fid, ['$E\\left( \\frac{1}{',num2str(nstep(end)),'}\\right)$ \\\\ \n' ]);

fprintf(fid, '\\hline\n');
% Write matrix data
for i=1:size(A,1)
     fprintf(fid, ['&' ]);
    for j = 1:size(A,2)-1
        fprintf(fid, [num2str(A(i,j),'%.2e'), '&', num2str(log(A(i,j)/A(i,j+1))/log(2),'%.2f'),'&' ]);
    end

    fprintf(fid, [num2str(A(i,end),'%.2e'),'\\\\ \n' ]);
end

% Write LaTeX table footer
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Generated LaTeX table}\n');
fprintf(fid, '\\label{tab:matrix}\n');
fprintf(fid, '\\end{table}\n');

% Close the file
fclose(fid);
fprintf('LaTeX table saved to %s\n', filename);
end
