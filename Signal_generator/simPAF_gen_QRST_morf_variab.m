function par = simPAF_gen_QRST_morf_variab(par, f, phase, rangeL, rangeH, n)
% 
% par = simPAF_gen_QRST_morf_variab() returns parameter values that are used for
% altering the simulated QRST morphology in time.
% 
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
    Td = 0.05;
    par = par + (rangeH - rangeL)*sin(2*pi*f*Td*n+phase)/2 + (rangeH - rangeL)*rand(1,1)/2;
end