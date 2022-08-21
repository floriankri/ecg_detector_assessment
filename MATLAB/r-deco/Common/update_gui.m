function handles = update_gui(handles)
% Updates the look of the GUI
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

%% Set menu settings
if ~isempty(handles.data.signal.original)
    set(handles.menu_save_session,'enable','on')
end

if ~isempty(handles.data.R)
    set(handles.menu_save_results,'enable','on')
    set(handles.menu_export_results,'enable','on')
    set(handles.menu_to_excel,'enable','on')
    set(handles.menu_to_workspace,'enable','on')
    
    % Enable the add context menu
    set(handles.context_add,'enable','on')
end

%% Set data panel
set(handles.text_filepath,'string',handles.data.filepath)
set(handles.text_variable,'string',handles.data.variable)
set(handles.text_sampling_frequency,'string',handles.data.fs)
set(handles.text_signal_duration,'string',char(handles.data.duration_recording))
set(handles.text_channels,'string',handles.data.channels)

set(handles.pushbutton_reset,'enable','on')

%% Set filter panel
set(findall(handles.uipanel_filter,'style','text'),'enable','on')
set(handles.edit_high_pass,'enable','on')
set(handles.edit_low_pass,'enable','on')
set(handles.edit_stop,'enable','on')
set(handles.pushbutton_check_psd,'enable','on')
set(handles.pushbutton_filter,'enable','on')

if ~isempty(handles.data.hp) || ~isempty(handles.data.lp) || ~isempty(handles.data.stop)
    set(handles.pushbutton_remove_filter,'enable','on')
    set(handles.checkbox_show_original,'enable','on')
end

%% Set analysis period panel
set(handles.edit_start,'string',char(handles.data.start_analysis));
set(handles.edit_duration,'string',char(handles.data.duration_analysis));

set(handles.togglebutton_analysis_period,'enable','on')


%% Set R-peak detection panel
set(handles.pushbutton_detect_peaks,'enable','on')
set(handles.pushbutton_load_peaks,'enable','on')


%% Set R-peak correction panel
if ~isempty(handles.data.R)
    set(findall(handles.uipanel_R_peak_correction,'-property','enable'),'enable','on')
end


%% Set range
set(handles.edit_range,'string',char(handles.data.range),...
    'enable','on')

set(handles.pushbutton_plus,'enable','on')
set(handles.pushbutton_minus,'enable','on')
set(handles.text_range,'enable','on')


%% Set slider
if handles.data.range == handles.data.duration_analysis
    % Adjust the slider
    set(handles.slider_range,'enable','off',...
        'max',round(handles.data.fs*seconds(handles.data.range)),...
        'value',0)

else
    % Get a new maximum for the slider
    new_max = round(handles.data.fs*seconds(handles.data.duration_analysis-handles.data.range));

    % Adjust the value if it is outside the new range
    if handles.slider_range.Value > new_max
        set(handles.slider_range,'Value',new_max)
    end

    % Get the new step size
    n = round(handles.data.fs*seconds(handles.data.range))/new_max;

    % Adjust the slider
    set(handles.slider_range,'max',new_max,...
        'sliderstep',[min([1 n]) min([1 2*n])],...
        'enable','on')
end

% Get the slider value
slider_value = get(handles.slider_range,'Value');

% Adjust the x-limits
stop = min([handles.data.start_analysis+seconds(slider_value/handles.data.fs)+handles.data.range handles.data.start_analysis+handles.data.duration_analysis]); 

if stop == handles.data.start_analysis+handles.data.duration_analysis
    if strcmp(handles.slider_range.Enable,'on')
        start = stop-range;
    else
        start = handles.data.start_analysis;
    end
else
    start = handles.data.start_analysis+seconds(slider_value/handles.data.fs);
end

% Adjust the x-limits
xlim(handles.axes_ECG,[start stop])
xlim(handles.axes_tachogram,[start stop])


%% Set channel listbox
str = cell(1,handles.data.channels+1);
str{1} = 'All';

for ii = 1:handles.data.channels
    color = handles.color(ii,:);
    colorStr = sprintf('%d,',int16(255*color)); 
    str{ii+1} = ['<HTML><FONT color="rgb(' colorStr ')">Channel_' num2str(ii) '</Font></html>'];
end

set(handles.listbox_channels,'string',str,...
    'value',1,...
    'enable','on')

%% Set the tachogram units
set(handles.popupmenu_units_tachogram,'enable','on')

%% Set the y-limits checkbox
set(handles.checkbox_fix_ylimits,'enable','on')

%% Set the axes
% Get the entire time
time = seconds((1:size(handles.data.signal.filtered,1))/handles.data.fs); 
time.Format = 'hh:mm:ss';

% Get the start and stop of the analysis window
start = max([1 round(handles.data.fs*seconds(handles.data.start_analysis))]);
stop = start + min([length(handles.data.signal.filtered) round(handles.data.fs*seconds(handles.data.duration_analysis))])-1;      

for ii = 1:handles.data.channels
    % Loaded signal
    handles.graph.signal.original.(handles.labels{ii}) = plot(handles.axes_ECG,time(start:stop),handles.data.signal.original(start:stop,ii),...
        'color',[handles.color(ii,:),0.2],...
        'visible','off');
    
    % Filtered signal
    handles.graph.signal.filtered.(handles.labels{ii}) = plot(handles.axes_ECG,time(start:stop),handles.data.signal.filtered(start:stop,ii),...
        'color',handles.color(ii,:),...
        'tag',handles.labels{ii},...
        'UIContextMenu', handles.context_lead_options);
    
    if ~isempty(handles.data.R)
        % Check if the input is in duration values. Adjust if necessary
        if ~isduration(handles.data.R{ii})
            handles.data.R{ii} = seconds(handles.data.R{ii})/handles.data.fs;
        end
        
        % Define R-peaks
        R_peak = handles.data.R{ii};
        
        % Define the RR-intervals
        RR_int = seconds(diff(R_peak))*1000;
        
        % R-peaks
        handles.graph.R.(handles.labels{ii}) = plot(handles.axes_ECG,...
                R_peak,...
                handles.data.signal.filtered(round(seconds(R_peak)*handles.data.fs),ii),...
                'o',...
                'color',handles.color(ii,:),...
                'tag',handles.labels{ii},...
                'markerfacecolor',handles.color(ii,:),...
                'UIContextMenu', handles.context_R_peak_options);
            
        % Tachogram
        handles.graph.RR.(handles.labels{ii}) = plot(handles.axes_tachogram,...
                R_peak(2:end),...
                RR_int,...
                '-o',...
                'color',handles.color(ii,:),...
                'markerfacecolor',handles.color(ii,:),...
                'tag',handles.labels{ii});
    end
end


