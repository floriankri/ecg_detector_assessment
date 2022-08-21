function create_excel_output(handles)

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

% Get all names
all_names = handles.listbox_channels.String;

% Get the units of the current tachogram
contents = cellstr(get(handles.popupmenu_units_tachogram,'String'));
value = get(handles.popupmenu_units_tachogram,'Value');
units = contents{value};      

% Activate excell
e = actxserver('Excel.application');

% Open a workbook
eWorkbook = e.Workbooks.Add; %#ok % New workbook

% Activate a worksheet
eAllSheets = e.ActiveWorkbook.Sheets;

% Loop over the ammount of channels
for ii = handles.data.channels:-1:1
    % Add a sheet if necessary
    if ii < handles.data.channels
        e.Sheets.Add;
    end
    eSheet = eAllSheets.Item(1);

    % Activate the sheet
    eSheet.Activate;

    % Change the sheet name
    L = length(all_names{ii+1});
    string = all_names{ii+1};
    name = string(find(string=='"',1,'last')+2:L-14);
    eSheet.Name = name;

    % Make the table
    eActivesheetRange = e.Activesheet.Range('A1:A6');
    eActivesheetRange.Value = {'Results' ; 'Nr beats:'; 'Mean HR (bpm):' ; 'STD HR (bpm):' ; 'Min HR (bpm):' ; 'Max HR (bpm):'};
    eActivesheetRange.Font.Bold = 1;
    eActivesheetRange.EntireColumn.AutoFit; % Autofit the cell width
    
    % Adjust the borders (all outer borders)
    eActivesheetRange = e.Activesheet.Range('A1:B6');
    for iii = 7:10
        eActivesheetRange.Borders.Item(iii).Linestyle = 1;
    end

    % Adjust the title
    eActivesheetRange = e.Activesheet.Range('A1:B1');
    eActivesheetRange.MergeCells = 1;
    eActivesheetRange.HorizontalAlignment = 3; % Set it to the middle
    for iii = 7:10
        eActivesheetRange.Borders.Item(iii).Linestyle = 1;
    end
    
    % Adjust the values
    if strcmp(units,'HR (bpm)')
        handles.graph.RR.(handles.labels{ii}).YData = 60*1000./handles.graph.RR.(handles.labels{ii}).YData;
    end
    
    % Fill in the values
    eActivesheetRange = e.Activesheet.Range('B2:B6');
    eActivesheetRange.Value = [length(handles.graph.R.(handles.labels{ii}).XData);...
        mean(60*1000./handles.graph.RR.(handles.labels{ii}).YData);...
        std(60*1000./handles.graph.RR.(handles.labels{ii}).YData);...
        min(60*1000./handles.graph.RR.(handles.labels{ii}).YData);...
        max(60*1000./handles.graph.RR.(handles.labels{ii}).YData)];
    
    % Create the value names
    eActivesheetRange = e.Activesheet.Range('D1:E1');
    eActivesheetRange.Value = {'R-peak location (hh:mm:ss:SSSS)','RR-interval (ms)'};
    eActivesheetRange.Font.Bold = 1;
    eActivesheetRange.EntireColumn.AutoFit; % Autofit the cell width
    
    % Get the ammount of R-peaks in the selected lead
    N = length(handles.graph.R.(handles.labels{ii}).XData);
    
    % Set the R-peak locations in the correct format
    temp = handles.graph.R.(handles.labels{ii}).XData';
    temp.Format = 'hh:mm:ss.SSSS';
    R_peaks = arrayfun(@char,temp,'UniformOutput',false);
    
    % Fill in the R-peak locations
    eActivesheetRange = e.Activesheet.Range(strcat('D2:D',num2str(N+1)));
    eActivesheetRange.Value = R_peaks;
    
    % Fill in the RR-intervals
    eActivesheetRange = e.Activesheet.Range(strcat('E3:E',num2str(N+1)));
    eActivesheetRange.Value = handles.graph.RR.(handles.labels{ii}).YData';
end

% Create the first sheet
e.Sheets.Add;
eSheet = eAllSheets.Item(1);

% Activate the sheet
eSheet.Activate;

% Change the sheet name
eSheet.Name = 'File information';

% Make the table
eActivesheetRange = e.Activesheet.Range('A1:A5');
eActivesheetRange.Value = {'File information' ; 'Nr channels:'; 'Sampling frequency (Hz):' ; 'Total duration recording (hh:mm:ss):' ; 'Total duration analysis period (hh:mm:ss):'};
eActivesheetRange.Font.Bold = 1;
eActivesheetRange.EntireColumn.AutoFit; % Autofit the cell width

% Adjust the borders (all outer borders)
eActivesheetRange = e.Activesheet.Range('A1:B5');
for ii = 7:10
    eActivesheetRange.Borders.Item(ii).Linestyle = 1;
end

% Adjust the title
eActivesheetRange = e.Activesheet.Range('A1:B1');
eActivesheetRange.MergeCells = 1;
eActivesheetRange.HorizontalAlignment = 3; % Set it to the middle
for ii = 7:10
    eActivesheetRange.Borders.Item(ii).Linestyle = 1;
end

% Fill in the values
eActivesheetRange = e.Activesheet.Range('B2:B5');
eActivesheetRange.Value = {handles.data.channels;...
    handles.data.fs;...
    char(handles.data.duration_recording);...
    char(handles.data.duration_analysis)};

% Make excell visible
e.visible = 1;
    
% Save, close and quit
% e.WorkBooks.Item('myResults').Save;
% e.WorkBooks.Item('myResults').Close;
% e.Quit;