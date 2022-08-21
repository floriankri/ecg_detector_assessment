function varargout = QRS_selection(varargin)
% QRS_SELECTION MATLAB code for QRS_selection.fig
%      QRS_SELECTION, by itself, creates a new QRS_SELECTION or raises the existing
%      singleton*.
%
%      H = QRS_SELECTION returns the handle to a new QRS_SELECTION or the handle to
%      the existing singleton*.
%
%      QRS_SELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QRS_SELECTION.M with the given input arguments.
%
%      QRS_SELECTION('Property','Value',...) creates a new QRS_SELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QRS_selection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QRS_selection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Edit the above text to modify the response to help QRS_selection

% Last Modified by GUIDE v2.5 11-Jun-2018 10:47:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QRS_selection_OpeningFcn, ...
                   'gui_OutputFcn',  @QRS_selection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT




% --- Executes just before QRS_selection is made visible.
function QRS_selection_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QRS_selection (see VARARGIN)

% Choose default command line output for QRS_selection
handles.output = hObject;

% Get sampling frequency
handles.fs = varargin{1,5};

% Get the beat length
handles.beat_length = length(varargin{1,1});

% Get the window length
handles.window = varargin{1,3};

% Set the selection based on the sign of peaks that was mostly present
if varargin{1,4} > 0
    handles.selection = 1;
    
    % Plot the QRS
    handles.QRS_pos = line(handles.axes_QRS,...
        1:handles.beat_length,varargin{1,1},...
        'Color','b',...
        'Linewidth',2,...
        'Tag','Positive',...
        'ButtonDownFcn', @LineCallback);
    handles.QRS_neg = line(handles.axes_QRS,...
        1:handles.beat_length,varargin{1,2},...
        'Color','r',...
        'Tag','Negative',...
        'ButtonDownFcn', @LineCallback);
else
    handles.selection = -1;
    
    % Plot the QRS
    handles.QRS_pos = line(handles.axes_QRS,...
        1:handles.beat_length,varargin{1,1},...
        'Color','b',...
        'Tag','Positive',...
        'ButtonDownFcn', @LineCallback);
    handles.QRS_neg = line(handles.axes_QRS,...
        1:handles.beat_length,varargin{1,2},...
        'Color','r',...
        'Linewidth',2,...
        'Tag','Negative',...
        'ButtonDownFcn', @LineCallback);
end

% Set the x- and y-limits
xlim(handles.axes_QRS,[0 handles.beat_length])
try
    ylim(handles.axes_QRS,[min(min([handles.QRS_pos.YData handles.QRS_neg.YData])) - min(std([handles.QRS_pos.YData handles.QRS_neg.YData]))...
        max(max([handles.QRS_pos.YData handles.QRS_neg.YData])) + min(std([handles.QRS_pos.YData handles.QRS_neg.YData]))])
catch
    try
        ylim(handles.axes_QRS,[min(handles.QRS_pos.YData) - std(handles.QRS_pos.YData)...
            max(handles.QRS_pos.YData) + std(handles.QRS_pos.YData)])
    catch
        try
            ylim(handles.axes_QRS,[min(handles.QRS_neg.YData) - std(handles.QRS_neg.YData)...
                max(handles.QRS_neg.YData) + std(handles.QRS_neg.YData)])
        catch
        end
    end
end

% Plot the R-peak line
handles.R_peak = line(handles.axes_QRS,...
        [round(handles.window*handles.fs)+1 round(handles.window*handles.fs)+1], ylim(handles.axes_QRS),...
        'Color','k',...
        'Linewidth',1,...
        'Tag','R_peak',...
        'ButtonDownFcn', @LineCallback);
    
% Create the legend    
legend(handles.axes_QRS,'Pos. R-peak','Neg. "R-peak"','R-peak location',...
    'Location','eastoutside')    

% Define the labels
ylabel(handles.axes_QRS,'ECG (a.u.)')
xlabel(handles.axes_QRS,'# samples')

% Enable the pointer manager and adjust the pointer
iptPointerManager(handles.figure1,'enable');

% Change the pointer behavior for the lines
pointerBehavior.enterFcn    = [];
pointerBehavior.exitFcn     = [];
pointerBehavior.traverseFcn = @overline;

iptSetPointerBehavior(handles.R_peak, pointerBehavior);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QRS_selection wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = QRS_selection_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.selection;
varargout{2} = handles.R_peak.XData(1) - (round(handles.window*handles.fs)+1);

% Delete the figure
delete(handles.figure1);


% --- Executes on button press in pushbutton_accept.
function pushbutton_accept_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.figure1)


% Helper functions
function LineCallback(hObject, eventdata) 
% hObject    handle to selected line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

% Get the handles
handles = guidata(eventdata.Source);

switch hObject.Tag
    case 'Positive'
        set(handles.QRS_pos,'Linewidth',2)
        set(handles.QRS_neg,'Linewidth',0.5)
        
        handles.selection = 1;
    case 'Negative'
        set(handles.QRS_pos,'Linewidth',0.5)
        set(handles.QRS_neg,'Linewidth',2)
        
        handles.selection = -1;   
    case 'R_peak'
        set(handles.QRS_pos,'Linewidth',0.5)
        set(handles.QRS_neg,'Linewidth',0.5)
        set(handles.R_peak,'Linewidth',2)
        
        % Disable the pointer manager
        iptPointerManager(handles.figure1,'disable');

        % Fix the pointer to left
        set(handles.figure1,'Pointer','left')
        
        set(handles.figure1,'WindowButtonMotionFcn',{@move_R_peak,handles},...
            'WindowButtonUpFcn',{@define_R_peak,handles})
        
end

% Update handles structure
guidata(hObject, handles);


function move_R_peak(hObject, eventdata, handles) %#ok

% Get the mouse position
mousepos = get(handles.axes_QRS,'CurrentPoint');

if mousepos(1) <= 1 % Past the left limit
    idx = 1;
elseif mousepos(1) >= handles.beat_length % Past the right limit
    idx = handles.beat_length;
else % Perfect
    idx = round(mousepos(1));
end

% Adjust the temporary R-peak
set(handles.R_peak,'XData',[idx idx])

% Update structure
guidata(hObject,handles);


function define_R_peak(hObject, eventdata, handles) %#ok

% Remove the window button motion and up function
set(handles.figure1,'WindowButtonMotionFcn','',...
    'WindowButtonUpFcn','');

% Change the linewidth back to previous
switch handles.selection
    case 1
        set(handles.QRS_pos,'Linewidth',2)
        set(handles.QRS_neg,'Linewidth',0.5)
        
    case -1
        set(handles.QRS_pos,'Linewidth',0.5)
        set(handles.QRS_neg,'Linewidth',2)
end

set(handles.R_peak,'Linewidth',1)

% Change the pointer back to it's original state
set(handles.figure1,'Pointer','arrow')

% Enable the pointer manager and adjust the pointer
iptPointerManager(handles.figure1,'enable');

% Change the pointer behavior for the lines
pointerBehavior.enterFcn    = [];
pointerBehavior.exitFcn     = [];
pointerBehavior.traverseFcn = @overline;

iptSetPointerBehavior(handles.R_peak, pointerBehavior);

% Update structure
guidata(hObject,handles);

