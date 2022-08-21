function handles = update_plots(handles)
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
% Get the listbox value
listbox_value = get(handles.listbox_channels,'value')-1; % 0=All, 1=Channel 1,...

% Define the looping variable
if listbox_value == 0
    loop_var = 1:handles.data.channels;
else % Only the selection
    loop_var = listbox_value;
end

% Define sampling frequency
fs = handles.data.fs;

% Get the start time
start = max([1 round(fs*seconds(handles.data.start_analysis))]);

% Get the stop time
stop = min([size(handles.data.signal.filtered,1) round(fs*seconds(handles.data.start_analysis+handles.data.duration_analysis))]);

% Get the time
time = seconds((1:size(handles.data.signal.filtered,1))/fs); 
time.Format = 'hh:mm:ss';

% Clear the axes
cla(handles.axes_ECG)
cla(handles.axes_tachogram)

% Adjust axis
for ii = 1:handles.data.channels
    % Loaded signal
    handles.graph.signal.original.(handles.labels{ii}) = plot(handles.axes_ECG,time(start:stop),handles.data.signal.original(start:stop,ii),...
        'color',[handles.color(ii,:),0.2],...
        'visible','off');

    % Filtered signal
    handles.graph.signal.filtered.(handles.labels{ii}) = plot(handles.axes_ECG,time(start:stop),handles.data.signal.filtered(start:stop,ii),...
        'color',handles.color(ii,:),...
        'visible','off',...
        'tag',handles.labels{ii},...
        'UIContextMenu', handles.context_lead_options);

    if ~isempty(handles.data.R)
        % Define the R-peaks that remain
        temp = handles.data.R{ii} > time(start);
        temp = temp + (handles.data.R{ii} < time(stop));
        R_peak = handles.data.R{ii}(temp == 2);
        
        if ~isempty(R_peak)
            % Define the RR-intervals
            RR_int = seconds(diff(R_peak))*1000;
            
            % Adjust the tachogram based on the units selection
            contents = get(handles.popupmenu_units_tachogram,'string');

            switch contents{get(handles.popupmenu_units_tachogram,'Value')}
                case 'HR (bpm)'
                    % Change the data
                    RR_int = 60000./RR_int;
            end

            % Plot the R-peaks on the ECG-graph
            handles.graph.R.(handles.labels{ii}) = plot(handles.axes_ECG,...
                R_peak,...
                handles.data.signal.filtered(round(seconds(R_peak)*fs),ii),...
                'o',...
                'color',handles.color(ii,:),...
                'markerfacecolor',handles.color(ii,:),...
                'visible','off',...
                'tag',handles.labels{ii},...
                'UIContextMenu', handles.context_R_peak_options);

            % Create the tachogram
            handles.graph.RR.(handles.labels{ii}) = plot(handles.axes_tachogram,...
                R_peak(2:end),...
                RR_int,...
                '-o',...
                'color',handles.color(ii,:),...
                'markerfacecolor',handles.color(ii,:),...
                'visible','off',...
                'tag',handles.labels{ii});
        end
    end
end

% Make the signals visible
for ii = loop_var
    set(handles.graph.signal.filtered.(handles.labels{ii}),'visible','on')

    if handles.checkbox_show_original.Value
        set(handles.graph.signal.original.(handles.labels{ii}),'visible','on')
    end

    if ~isempty(handles.data.R)
        set(handles.graph.R.(handles.labels{ii}),'visible','on')
        set(handles.graph.RR.(handles.labels{ii}),'visible','on')
    end
end

