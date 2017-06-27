function saveTable(tblName,x0,iter1,iter2,evals1,evals2,T1,T2)
T1 = 1000*T1; T2 = 1000*T2;

file = fopen(['../tables/' tblName],'w');
fprintf(file,'\\begin{tabular}{l|ccc|ccc|ccc} \\hline \\hline \n');
fprintf(file,'& \\multicolumn{3}{c}{$x_0 = (%d,%d)$} & \\multicolumn{3}{c}{$x_0 = (%d,%d)$} & \\multicolumn{3}{c}{$x_0 = (%d,%d)$} \\\\ \n',x0(1,1),x0(1,2),x0(2,1),x0(2,2),x0(3,1),x0(3,2));
fprintf(file,'Method & Iter & Evals & Time & Iter & Evals & Time & Iter & Evals & Time \\\\ \\hline \n');

fprintf(file,'Backtracking & %d & %d & %.1f & %d & %d & %.1f & %d & %d & %.1f \\\\ \n',iter1(1),evals1(1),T1(1),iter1(2),evals1(2),T1(2),iter1(3),evals1(3),T1(3));
fprintf(file,'Soft line search & %d & %d & %.1f & %d & %d & %.1f & %d & %d & %.1f \\\\ \n',iter2(1),evals2(1),T2(1),iter2(2),evals2(2),T2(2),iter2(3),evals2(3),T2(3));

fprintf(file,'\\hline \\hline \n');
fprintf(file,'\\end{tabular} \n');
fclose(file);