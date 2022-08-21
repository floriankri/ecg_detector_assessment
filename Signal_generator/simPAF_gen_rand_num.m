function number = simPAF_gen_rand_num(rangeLow, rangeHigh)
%
% number = simPAF_gen_rand_num() returns a random number from the interval
% [rangeLow, rangeHigh].
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
number = (rangeHigh-rangeLow)*rand(1,1) + rangeLow;
number = round(number*100)/100;
end