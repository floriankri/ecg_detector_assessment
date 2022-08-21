% Adjust the paper position mode
handles.figure_main.PaperPositionMode = 'auto';

% Save R-DECO figure
print(handles.figure_main, 'test', '-depsc ','-painters');

