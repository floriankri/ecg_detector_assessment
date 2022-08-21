function pqrst = simPAF_correct_baseline(pqrstIn)
% 
% [] = simPAF_correct_baseline() corrects PQRST baseline.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
    pqrst(1,:) = pqrstIn;
    x = [1 length(pqrst)];
    xi = 1:1:length(pqrst);
    y = [pqrst(1,1) pqrst(1,end)];
    baseline = interp1(x,y,xi,'linear');
    pqrst = pqrst - baseline;
end