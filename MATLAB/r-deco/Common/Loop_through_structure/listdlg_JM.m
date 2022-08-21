function [selection,ok] = listdlg_JM(varargin)
%LISTDLG  List selection dialog box.
%   [SELECTION,OK] = LISTDLG_JM('ListString',S) creates a modal dialog box
%   which allows you to select a string or multiple strings from a list.
%   SELECTION is a vector of indices of the selected strings (length 1 in
%   the single selection mode).  This will be [] when OK is 0.  OK is 1 if
%   you push the OK button, 2 if you push the Back button, or 0 if you push
%   the Cancel button or close the figure.
%
%   Double-clicking on an item or pressing <CR> when multiple items are
%   selected has the same effect as clicking the OK button.  Pressing <CR>
%   is the same as clicking the OK button. Pressing <ESC> is the same as
%   clicking the Cancel button.
%
%   Inputs are in parameter,value pairs:
%
%   Parameter       Description
%   'ListString'    cell array of strings for the list box.
%   'SelectionMode' string; can be 'single' or 'multiple'; defaults to
%                   'multiple'.
%   'ListSize'      [width height] of listbox in pixels; defaults
%                   to [160 300].
%   'InitialValue'  vector of indices of which items of the list box
%                   are initially selected; defaults to the first item.
%   'Name'          String for the figure's title; defaults to ''.
%   'PromptString'  string matrix or cell array of strings which appears
%                   as text above the list box; defaults to {}.
%   'OKString'      string for the OK button; defaults to 'OK'.
%   'BackString'    string for the Back button; defaults to 'Back'.
%   'CancelString'  string for the Cancel button; defaults to 'Cancel'.
%   'EnableBack'    string to enable or disable the Back button
%
%   A 'Select all' button is provided in the multiple selection case.
%
%   Example:
%     d = dir;
%     str = {d.name};
%     [s,v] = listdlg_JM('PromptString','Select a file:',...
%                     'SelectionMode','single',...
%                     'ListString',str)
%
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG,
%    MSGBOX, QUESTDLG, WARNDLG.

%   Copyright 1984-2014 The MathWorks, Inc.

% simple test:
%
% d = dir; [s,v] = listdlg_JM('PromptString','Select a file:','ListString',{d.name});
%

% Check the ammount of inputs
narginchk(1,inf)

% Set the defaults
figname = '';
smode = 2;   % (multiple)
promptstring = {};
liststring = [];
listsize = [160 300];
initialvalue = [];
okstring = getString(message('MATLAB:uistring:popupdialogs:OK'));
backstring = 'Back';
cancelstring = getString(message('MATLAB:uistring:popupdialogs:Cancel'));
enableback = 'off';
fus = 8; % frame/uicontrol spacing, in pixels; default = 8.
ffs = 8; % frame/figure spacing, in pixels; default = 8.
uh = 22; % uicontrol button height, in pixels; default = 22.

if mod(length(varargin),2) ~= 0
    % input arguments have not come in pairs, create an error message
    error(message('MATLAB:listdlg:InvalidArgument'))
end

% Go one input at a time
for i=1:2:length(varargin)
    switch lower(varargin{i}) % Make capitals of the text small
        case 'name'
            figname = varargin{i+1};
        case 'promptstring'
            promptstring = varargin{i+1};
        case 'selectionmode'
            switch lower(varargin{i+1})
                case 'single'
                    smode = 1;
                case 'multiple'
                    smode = 2;
            end
        case 'listsize'
            listsize = varargin{i+1};
        case 'liststring'
            liststring = varargin{i+1};
        case 'initialvalue'
            initialvalue = varargin{i+1};
        case 'uh'
            uh = varargin{i+1};
        case 'fus'
            fus = varargin{i+1};
        case 'ffs'
            ffs = varargin{i+1};
        case 'okstring'
            okstring = varargin{i+1};
        case 'backstring'
            backstring = varargin{i+1};
        case 'cancelstring'
            cancelstring = varargin{i+1};
        case 'enableback'
            enableback = varargin{i+1};
        otherwise
            error(message('MATLAB:listdlg:UnknownParameter', varargin{ i }))
    end
end

if ischar(promptstring)
    promptstring = cellstr(promptstring);
end

if isempty(initialvalue)
    initialvalue = 1;
end

if isempty(liststring)
    error(message('MATLAB:listdlg:NeedParameter'))
end

ex = get(0,'DefaultUicontrolFontSize')*1.7;  % height extent per line of uicontrol text (approx)

fp = get(0,'DefaultFigurePosition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

fig_props = { ...
    'windowstyle'            'modal' ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'visible'                'off' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
    };

liststring=cellstr(liststring);

fig = figure(fig_props{:});

if ~isempty(promptstring)
    prompt_text = uicontrol('Style','text','String',promptstring,...
        'HorizontalAlignment','left',...
        'Position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring)) ...
        listsize(1) ex*length(promptstring)]); %#ok
end

% Define button width
btn_wid = (fp(3)-2*(ffs+fus)-2*fus)/3;

% Define the uicontrols
listbox = uicontrol('Style','listbox',...
    'Position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
    'String',liststring,...
    'BackgroundColor','w',...
    'Max',smode,...
    'Tag','listbox',...
    'Value',initialvalue, ...
    'Callback', {@doListboxClick});

ok_btn = uicontrol('Style','pushbutton',...
    'String',okstring,...
    'Position',[ffs+fus ffs+fus btn_wid uh],...
    'Tag','ok_btn',...
    'Callback',{@doOK,listbox});

back_btn = uicontrol('Style','pushbutton',...
    'String',backstring,...
    'Position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
    'Tag','ok_btn',...
    'Enable',enableback,...
    'Callback',{@doBack,listbox});

cancel_btn = uicontrol('Style','pushbutton',...
    'String',cancelstring,...
    'Position',[ffs+3*fus+2*btn_wid ffs+fus btn_wid uh],...
    'Tag','cancel_btn',...
    'Callback',{@doCancel,listbox});

% If a select all button is required
if smode == 2
    selectall_btn = uicontrol('Style','pushbutton',...
        'String',getString(message('MATLAB:uistring:popupdialogs:SelectAll')),...
        'Position',[ffs+fus 4*fus+ffs+uh listsize(1) uh],...
        'Tag','selectall_btn',...
        'Callback',{@doSelectAll, listbox});
    
    % If only one value is present
    if length(initialvalue) == length(liststring)
        set(selectall_btn,'Enable','off')
    end
    set(listbox,'Callback',{@doListboxClick, selectall_btn})
end

% Define the keypress functions
set([fig, ok_btn, cancel_btn, back_btn, listbox], 'KeyPressFcn', {@doKeypress, listbox});

% Define a nice dialog location
set(fig,'Position',getnicedialoglocation(fp, get(fig,'Units')));

% Make ok_btn the default button.
fig.setDefaultButton(ok_btn);

% make sure we are on screen
movegui(fig)
set(fig, 'Visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(listbox);
    c = matlab.ui.internal.dialog.DialogUtils.disableAllWindowsSafely();
    uiwait(fig);
    delete(c);
catch
    if ishghandle(fig)
        delete(fig)
    end
end

% Define the output
if isappdata(0,'ListDialogAppData__')
    ad = getappdata(0,'ListDialogAppData__');
    selection = ad.selection;
    ok = ad.value;
    rmappdata(0,'ListDialogAppData__')
else
    % figure was deleted
    selection = [];
    ok = 0;
end
drawnow; % Update the view to remove the closed figure (g1031998)

%% figure, OK, Back and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
    case 'escape'
        doCancel([],[],listbox);
end

%% OK callback
function doOK(ok_btn, evd, listbox) %#ok
if (~isappdata(0, 'ListDialogAppData__'))
    ad.value = 1;
    ad.selection = get(listbox,'Value');
    setappdata(0,'ListDialogAppData__',ad);
    delete(gcbf);
end

%% Back callback
function doBack(back_btn, evd, listbox) %#ok
ad.value = 2;
ad.selection = [];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);

%% Cancel callback
function doCancel(cancel_btn, evd, listbox) %#ok
ad.value = 0;
ad.selection = [];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox) %#ok
set(selectall_btn,'Enable','off')
set(listbox,'Value',1:length(get(listbox,'String')));

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);
elseif nargin == 3
    if length(get(listbox,'String'))==length(get(listbox,'Value'))
        set(selectall_btn,'Enable','off')
    else
        set(selectall_btn,'Enable','on')
    end
end
