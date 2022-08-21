function par = simPAF_gen_P_morf_variab(f, phase, par, rangeLow, rangeHigh, n)
% 
% par = simPAF_gen_P_morf_variab() returns parameter values that are used for
% altering simulated P wave morphology in time.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
    Td = 0.05;
    par = par + (2*par*(rand(1,1)-0.5)/10) + ((rangeHigh-rangeLow)*(1+sin(2*pi*f*Td*n+phase))/2+rangeLow)/10;
end