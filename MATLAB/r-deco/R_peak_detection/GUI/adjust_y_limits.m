function handles = adjust_y_limits(handles,ax)

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

if nargin == 1
    ax = {'axes_ECG','axes_tachogram'};
end

% Get the fix y-limits checkbox value
checkbox_value = get(handles.checkbox_fix_ylimits,'Value');

% Set up the axes
if checkbox_value && length(ax) == 2
    ax = {'axes_ECG'};
elseif checkbox_value && strcmp(ax,'axes_tachogram')
    ax = [];
end

% Loop over the axes
for ii = 1:length(ax)
    % Get the line children
    children = findall(handles.(ax{ii}),'Visible','on','Type','Line',...
        'LineStyle','-','-not','Tag','');

    % Get the ammount of visible children
    L = numel(children);

    % Only continue if children are present
    if L >= 1
        % Get the x-limits
        Xlim = xlim(handles.(ax{ii}));
        
        % Pre-allocate
        Min = zeros(1,L);
        Max = zeros(1,L);
        Std = zeros(1,L);

        % Loop over the children
        for iii = 1:L
            % Get the start and finish
            start = find(children(iii).XData >= Xlim(1),1,'first');
            finish = find(children(iii).XData <= Xlim(2),1,'last');

            if isempty(start)
                start = finish;
            elseif isempty(finish)
                finish = start;
            elseif start > finish
                temp = start;
                start = finish;
                finish = temp;
            end

            % Get the minimum, maximum and standard deviation
            Min(iii) = min(children(iii).YData(start:finish));
            Max(iii) = max(children(iii).YData(start:finish));
            Std(iii) = std(children(iii).YData(start:finish));
        end

        % Adjust the y-limits
        if ~(min(Min) == max(Max))
            ylim(handles.(ax{ii}),[min(Min)-max(Std) max(Max)+max(Std)])
        end
    end
end