function [fs , check] = get_sampling_frequency

% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 06/05/2019   JM      Initial version
%
% Copyright (c) 2019,  Jonathan Moeyersons, KULeuven-ESAT-STADIUS 
%
% This software is made available for non commercial research purposes only 
% under the GNU General Public License. However, notwithstanding any 
% provision of the GNU General Public License, this software may not be 
% used for commercial purposes without explicit written permission after 
% contacting jonathan.moeyersons@esat.kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Pre-allocate
fs = [];
check = 0;

% Set dialog properties
prompt = {'Enter sampling frequency (Hz):'};
dlg_title = 'Sampling frequency';
num_lines = 1;
defaultans = {'250'};

% Ask to enter the sampling frequency
while check == 0
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

    % Perform a check of the sampling frequency
    if isempty(answer) || sum(~isnan(str2double(answer{1})) && ~sum(contains(answer{1},[",","."])))
        check = 1;
    end
end

% Define sampling frequency, this will fail if cancel is pressed or the dialog box is closed 
try
    fs = str2double(answer{1}); 
catch
    check = 0;
end