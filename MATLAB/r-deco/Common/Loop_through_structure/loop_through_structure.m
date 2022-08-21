function [data,check] = loop_through_structure(data)
% Input
% data       - Structure that contains the data
% 
% Output
% data       - The data that the variable contains
% check      - 1 if the dialog is excited normally, 0 if closed
%
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


% Disable the back button
enableback = 'off';

% Pre-allocate depth in structure
depth = 1;

% Get the fieldnames
vars{depth} = fieldnames(data);

% Define the temporary data
temp_data = data;

% Loop until there is no more structure
while isa(temp_data,'struct')

    % Show data structure
    [temp_selection,check] = listdlg_JM('PromptString','Select the ECG file:',...
        'SelectionMode','single',...
        'ListString',vars{depth},...
        'EnableBack',enableback);

    if check == 1 % OK button is pressed   
        % Adjust selection
        selection(depth) = temp_selection; %#ok
        
        % Change the data variable
        temp_data = temp_data.(vars{depth}{selection(depth)});

        % Change depth
        depth = depth+1;

        % Enable the back button
        enableback = 'on';

        % Get the fieldnames, if this is not possible, we have reached the
        % end
        try
            vars{depth} = fieldnames(temp_data);
        catch
            break
        end
    elseif check == 2 % Back button is pressed
        % Change depth
        depth = depth-1;
        
        % Change the data variable
        if depth == 1 % Back to the start
            temp_data = data;

            % Disable the back button
            enableback = 'off';
        else % One step back
            temp_data = data;

            for ii = 1:length(depth-1)
                temp_data = temp_data.(vars{ii}{selection(ii)});
            end
        end
    else
        break
    end
end

if check == 1
    
    % Change the original data variable
    data = temp_data;
end
