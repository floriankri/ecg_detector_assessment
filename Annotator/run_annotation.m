function run_annotation(data)
% run_mimic_imp_annotation is designed to annotate the quality and breaths
% of signals in the impedance pneumography signal quality index (SQI) dataset
% (extracted from the MIMIC III database). It was used to evaluate the
% impedance pneumography SQI (using only the breath annotation functionality).
%
%               run_mimic_imp_annotation
%
%	Inputs:
%       mimic_imp_sqi_data.mat - a file containing curated data from the MIMIC database.
%        - Specify the path to this file in the "setup_universal_params" function below.
%
%	Outputs:
%       Matlab files containing the annotations for each subject, called
%          "ImP_SQI#_##_an.mat" (where #_## is the subject number, followed
%          by the initials of the annotator)
%           
%   Further Information:
%       This version of the run_mimic_imp_annotation is provided to facilitate
%       replication of the analysis reported in:
%           Charlton P. H. et al., "An impedance pneumography signal
%           quality index for respiratory rate monitoring: design,
%           assessment and application", [under review]
%       Further information on this study can be obtained at:
%           http://peterhcharlton.github.io/RRest/imp_sqi.html
%       In addition, further information on RRest, a toolbox of respiratory
%       algorithms which can be used with this dataset, can be obtained at:
%           http://peterhcharlton.github.io/RRest/index.html
%
%   Comments, Questions, Criticisms, Feedback, Contributions:
%       See: http://peterhcharlton.github.io/RRest/contributions.html
%
%   Version:
%       v.0.1 - accompanying peer review, 6th Aug 2020 by Peter H Charlton
%
%   Source:
%       This script has been adapted from 'MIMICII_data_importer.m', which
%       is one of the scripts in the RRest toolbox, available at:
%           https://github.com/peterhcharlton/RRest
%
%   Licence:
%       This program is available under the GNU public license. This
%       program is free software: you can redistribute it and/or modify 
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version. This program is distributed in
%       the hope that it will be useful, but WITHOUT ANY WARRANTY; without
%       even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%       PARTICULAR PURPOSE.  See the GNU General Public License for more
%       details: <http://www.gnu.org/licenses/>.
%

%% Setup
% ~~~ This function needs adjusting by the individual user as it contains paths
% specific to their computer ~~~
up = setup_universal_params;

%% Initialisation

% Load data
load(up.paths.data_file)
subjects = 1:length(data);   % create a list of subjects which are to be analysed.

% mode
up.mode = 2;        % select 1 if annotating high or low quality, 2 if annotating breaths
if up.mode == 1, fprintf('\n *** Mode selected: quality annotation *** '), elseif up.mode == 2, fprintf('\n *** Mode selected: Breath annotation *** '), end

fprintf(['\n\n\n\n\n~~~~~~~  Starting annotations  ~~~~~~~\n' ...
    'Instructions: Use the left mouse key to\n' ...
    'annotate a peak, and use the right mouse\n' ...
    'key to annotate a trough. To delete the\n' ...
    'most recent annotation, hold down SHIFT\n' ...
    'and press a mouse key.\n' ...
    'Please note that any previous annotations\n' ...
    'will be kept as long as the results files\n' ...
    'are in the correct location.\n\n']);

operator=input('(To enter review mode, type  review  )\n\nOtherwise, please type your first and last initial,\nand press enter:\n','s');
%operator = 'PC';
if strcmp(operator, 'review') == 1
    paw_verification_script_checker(subjects, data, up);
    return
end

question=input('\nWould you like to start from the \nfirst subject? (y or n, followed by enter)\n', 's');
%question = 'y';
if strcmp(question, 'n') == 1
    start_subject = str2num(input('\nWhich subject?\n', 's'));
else
    start_subject = 1;
end

%% Cycle through different subjects and periods
start_subject_el = find(subjects == start_subject);
if isempty(start_subject_el)
    fprintf('\nThis subject doesn''t exist. Try all over again.');
else
    subjects = subjects(start_subject_el:end);
    for s = subjects(1:end)
        
        %% Load and process ECG signal
        ekg = load_and_process_ecg(data, s, up);
        
        %% Set up file for keeping annotations in:
        % see if this subject has already been annotated:
        if up.mode == 1
            qual.t = []; qual.v = [];
            save(up.paths.annotations_file, 'qual', 'up')
        else
            pk_anns.t=[]; pk_anns.v=[];
            save(up.paths.annotations_file, 'pk_anns', 'up')
        end 
        
        %% plot waveforms and annotations (for further annotation)
        duration_of_signal = ekg.t(end) - ekg.t(1);
        NUMWINS = floor(duration_of_signal / up.win_length); clear duration_of_signal
        edge_time = up.edge_prop*up.win_length;
        for win = 1:NUMWINS
            
            % Find window timings
            starttime = ekg.t(1)+(win-1)*up.win_length;
            endtime = starttime + up.win_length;
            
            % Skip if this window has already been annotated
            savename = [up.paths.annotations_folder, 'ImP_SQI' num2str(s), '_' operator '_an.mat'];
            if exist(savename, 'file')
                temp = load(savename);
                rel_els = temp.pk_anns.t>=starttime & temp.pk_anns.t<=endtime;
                min_diff = min(diff(temp.pk_anns.t(rel_els)));
                if min_diff < 1
                    error(['Check window ', num2str(win), ', subj ', num2str(s)])
                end
                if sum(rel_els)>0
                    continue
                end
            end
            
            % Setup plot
            rel_els = find(ekg.t >= starttime & ekg.t < endtime);
            grey_els = find(ekg.t >= (starttime - edge_time) & ekg.t < (endtime + edge_time));
            scrsz = get(0,'ScreenSize');
            load(up.paths.annotations_file)
            
            % Plot
            h=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
            %save('plot_vars.mat', 'h', 'ekg', 'Re', 'grey_paw_els', 'rel_paw_els', 'rel_Re_els', 'grey_Re_els', 'period', 's', 'subjects')
            [axis_h, h2] = annotation_plotter2(up.paths.annotations_file, h, ekg, grey_els, rel_els, s, subjects, up, win, NUMWINS,starttime);
            
            % use pointer to identify breaths
            set(gcf, 'Pointer', 'cross');
            if up.mode == 1
                % Construct a questdlg with three options
                choice = MFquestdlg([0.1,0.1],'High quality?', ...
                    '', ...
                    'Yes','No','No');
                % Handle response
                switch choice
                    case 'Yes'
                        load(up.paths.annotations_file);
                        [qual.t,inds]=unique([qual.t; starttime]);
                        qual.v=[qual.v; 1]; qual.v = qual.v(inds);
                        [qual.t, inds] = sort(qual.t);
                        qual.v = qual.v(inds);
                        save(up.paths.annotations_file, 'qual');
                    case 'No'
                        load(up.paths.annotations_file);
                        [qual.t,inds]=unique([qual.t; starttime]);
                        qual.v=[qual.v; 0]; qual.v = qual.v(inds);
                        [qual.t, inds] = sort(qual.t);
                        qual.v = qual.v(inds);
                        save(up.paths.annotations_file, 'qual');
                end
            else
                set(h2,'ButtonDownFcn',{@Click_CallBack2 gca axis_h h2 up up.paths.annotations_file starttime});
                pause
            end
            close(h);
            
        end
        
        if up.DOSAVE
            fprintf(['~~~~~~ Saving subject ' num2str(s) '  ~~~~~~\n']);
            if exist(savename, 'file')
                old_data = load(savename);
            end
            new_data.up = up;
            temp = load(up.paths.annotations_file);
            % If there are annotations in both the old file and the current annotations file
            if exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'pk_anns')) && ~isempty(old_data.pk_anns.t) ...
                    && sum(strcmp(fieldnames(temp), 'pk_anns')) && ~isempty(temp.pk_anns.t)
                all_t = [temp.pk_anns.t; old_data.pk_anns.t];
                all_v = [temp.pk_anns.v; old_data.pk_anns.v];
                [~, order] = sort(all_t);
                new_data.pk_anns.t = all_t(order);
                new_data.pk_anns.v = all_v(order);
                clear order all_t all_v
            % If there are annotations in the old file
            elseif exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'pk_anns')) && ~isempty(old_data.pk_anns.t)
                new_data.pk_anns = old_data.pk_anns;
            % If there are annotations in the current annotations file
            elseif sum(strcmp(fieldnames(temp), 'pk_anns')) && ~isempty(temp.pk_anns.t)
                new_data.pk_anns = temp.pk_anns;
            else
                [new_data.pk_anns.t, new_data.pk_anns.v] =deal([]);
            end
            if exist('qual', 'var') && ~isempty(qual.t)
                new_data.qual = qual;
            elseif sum(strcmp(fieldnames(temp), 'qual'))  && ~isempty(temp.qual.t)
                new_data.qual = temp.qual;
            elseif exist('old_data', 'var') && sum(strcmp(fieldnames(old_data), 'qual'))  && ~isempty(old_data.qual.t)
                new_data.qual = old_data.qual;
            else
                new_data.qual.t = []; qual.v = [];
            end
            pk_anns = new_data.pk_anns;
            qual = new_data.qual;
            up = new_data.up;
            save(savename, 'pk_anns', 'qual', 'up')
            clear old_data new_data temp pk_anns qual
        end
        
        clear pk_anns qual NUMWINS pawmean STUDYID
        
    end
    
end

end

function up = setup_universal_params

close all

%% folderpaths
up.paths.root_data_folder = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\Annotator\data\';
up.paths.data_file = [up.paths.root_data_folder, 'MIMIC_PERform_truncated_train_all_data.mat'];
up.paths.resultspath = [up.paths.root_data_folder, 'SQI_development', filesep];
up.paths.annotations_folder = [up.paths.root_data_folder, '2019_annotations', filesep];
up.paths.annotations_file = [up.paths.root_data_folder, '2019_annotations', filesep, 'temp.mat'];
up.paths.dlg = 'C:\Users\flori\OneDrive\Dokumente\TU\Bachelor Thesis\Code\Annotator\required\';
addpath(up.paths.dlg)
if ~exist(up.paths.annotations_folder, 'dir')
    mkdir(up.paths.annotations_folder)
end

%% Operational specifications
up.DOSAVE = 1;                  % Save plots of filter characteristics?
up.edge_prop = 0.3;             % proportion of window to show either side of the one currently being viewed.
up.win_length = 10;            % window length for analysis in secs

%% Filtering specifications
% - eliminating low freqs
up.elim_low_freqs.Fpass = 1.02;  % in Hz
up.elim_low_freqs.Fstop = 0.6;  % in Hz      % 0.6 - 1.02 gives -3dB cutoff of 40 bpm (i.e. 0.67 Hz)
up.elim_low_freqs.Dpass = 0.05;
up.elim_low_freqs.Dstop = 0.01;

end

function [axis_h, h2] = annotation_plotter2(annotations_filepath, h, imp, grey_els, rel_els, s, subjects, up, win, NUMWINS, starttime)

%% ------------    annotation_plotter  ------------------------------
%
% Created by: Peter Charlton (17/10/2012)
%
% Purpose: To plot the Paw and impedance waveforms, and the annotated
%           breaths.
%
% Context: Some safety nets have been added to make it run ok even if the
% impedance signal isn't recorded.
%
% Inputs:
%
% Outputs:
%
% Files required:
%
% Called by:            Run_paw_verification
%
% ------------------------------------------------------------------------

ftsize = 16; lwidth = 2;
init_time = imp.t(rel_els(1));

load(annotations_filepath)
mean_imp = mean(imp.v(rel_els));
h3 = plot(imp.t(rel_els) - init_time, zeros(1,length(rel_els)), 'k', 'LineWidth', lwidth); hold on    % plot zero line on Paw plot
set(h3, 'HandleVisibility', 'off')
h1 = plot(imp.t(grey_els) - init_time, imp.v(grey_els) - mean_imp, 'Color', [0.5 0.5 0.5], 'LineWidth', lwidth); hold on,    % plot Paw (grey)
set(h1, 'HandleVisibility', 'off')
h2 = plot(imp.t(rel_els) - init_time, imp.v(rel_els) - mean_imp, 'LineWidth', lwidth);    % plot Paw
set(h2, 'HandleVisibility', 'off')
xlim([imp.t(grey_els(1)) - init_time, imp.t(grey_els(end)) - init_time])

%% Calculate and plot mix signal
if up.mode == 2
    % plot annotations which are relevant to this plot:
    %breathels = find(pk_anns.t>=imp.t(grey_els(1)) & pk_anns.t<=imp.t(grey_els(end)));
    %plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels)-mean_imp, 'ro','LineWidth',4)
    breathels = find(pk_anns.t>=imp.t(rel_els(1)) & pk_anns.t<=imp.t(rel_els(end)));
    plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels), 'ro','LineWidth',4)
    
    title(['Please annotate the peaks (red) then press a button  -  ' num2str(win) ' of ' num2str(NUMWINS) ' windows for subject ' num2str(s) ' of ' num2str(subjects(end))], 'FontSize', ftsize)
    ylabel('ECG', 'FontSize', ftsize), xlabel('Time [s]', 'FontSize', ftsize);
    set(gca,'XTick',((imp.t(rel_els(1)) - init_time):5:ceil((imp.t(rel_els(end)) - init_time))), 'YTick', [])
    set(gca, 'FontSize', ftsize)
else
    title(['Please annotate high quality (left) or low quality (right)  -  ' num2str(win) ' of ' num2str(NUMWINS) ' windows for subject ' num2str(s) ' of ' num2str(subjects(end))], 'FontSize', ftsize)
    ylabel('Impedance', 'FontSize', ftsize), xlabel('Time [s]', 'FontSize', ftsize);
    ylabel('Imepdance (above line = inhalation)', 'FontSize', ftsize)
    set(gca,'XTick',((imp.t(rel_els(1)) - init_time):5:ceil((imp.t(rel_els(end)) - init_time))), 'YTick', [])
    set(gca, 'FontSize', ftsize)
    if find(qual.t == starttime)
        rel_qual = qual.v(qual.t == starttime);
        if rel_qual
            annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','High Quality', 'Color', 'b', 'FontSize', ftsize*3, 'LineStyle', 'None')
        else
            annotation('textbox',[0.75, 0.2, 0.1,0.1],'String','Low Quality', 'Color', 'r', 'FontSize', ftsize*3, 'LineStyle', 'None')
        end
    end
    if exist('pk_anns', 'var')
        % plot annotations which are relevant to this plot:
        breathels = find(pk_anns.t>=imp.t(rel_els(1)) & pk_anns.t<=imp.t(rel_els(end)));
        plot(pk_anns.t(breathels) - init_time,pk_anns.v(breathels), 'ro','LineWidth',4)
    end
end

axis_h = gca;

end

function [pk_anns, TRS]= Click_CallBack2(h,e, a, axis_h, h2, up, annotations_filepath, start_time)

switch get(ancestor(a,'figure'),'SelectionType')
    
    case 'normal' %left click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)+start_time]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'alt'  % right click
        point = get(a,'CurrentPoint');
        load(annotations_filepath);
        [pk_anns.t,inds]=unique([pk_anns.t; point(1,1)]);
        pk_anns.v=[pk_anns.v; point(1,2)]; pk_anns.v = pk_anns.v(inds);
        [pk_anns.t, inds] = sort(pk_anns.t);
        pk_anns.v = pk_anns.v(inds);
        save(annotations_filepath, 'pk_anns');
        
    case 'extend' % right click whilst holding down shift
        load(annotations_filepath)
        pk_anns.t = pk_anns.t(1:(end-1));
        pk_anns.v = pk_anns.v(1:(end-1));
        save(annotations_filepath, 'pk_anns');
        
end

cla(axis_h)
plot(pk_anns.t-start_time,pk_anns.v, 'ro','LineWidth',4)

end

function paw_verification_script_checker(subjects, data, up)
% this function is used to check that the results can indeed be plotted
% back on the original waveform:
%% initialisation
fprintf('\n\n\n\n\n~~~~~~~  Reviewing Annotations  ~~~~~~~\n');

%% find out the file details
operator=input('Please type the first and last initial,\nof the operator whose annotations\nyou''d like to review, and press enter:\n','s');
question=input('\nWould you like to start from the \nfirst subject? (y or n, followed by enter)\n', 's');
if strcmp(question, 'n') == 1
    start_subject = input('\nWhich subject?\n', 's');
    start_subject = str2double(start_subject);
else
    start_subject = 1;
end

start_subject_el = find(subjects == start_subject);
if isempty(start_subject_el)
    fprintf('\nThis subject doesn''t exist. Try all over again.');
else
    subjects = subjects(start_subject_el:end);
    for s = subjects(1:end)

        %% Load and process ECG signal
        ekg = load_and_process_ecg(data, s, up);
        
        %% Load annotations
        up_copy = up;
        load(up.paths.annotations_file);
        up = up_copy; clear up_copy
        
        %% plot waveforms and annotations
        duration_of_signal = ekg.t(end) - ekg.t(1);
        NUMWINS = floor(duration_of_signal / up.win_length); clear duration_of_signal
        edge_time = up.edge_prop*up.win_length;
        for win = 1:NUMWINS
            % Setup
            starttime = ekg.t(1)+edge_time+(win-1)*up.win_length;
            endtime = starttime + up.win_length;
            rel_els = find(ekg.t >= starttime & ekg.t < endtime);
            grey_els = find(ekg.t >= (starttime - edge_time) & ekg.t < (endtime + edge_time));
            scrsz = get(0,'ScreenSize');
            
            % Plot
            h=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
            [axis_h, h2] = annotation_plotter2(up.paths.annotations_file, h, ekg, grey_els, rel_els, s, subjects, up, win, NUMWINS, starttime);
            
            % close plot after a pause
            pause
            close all
        end
    end
end

end

function s_filt = elim_vlfs(s, filt_characteristics)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 4.435e-3;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

if length(s.v) > (length(AMfilter.numerator)-1)*3
    % - Tukey Window to avoid edge effects
    win_edge_durn = 0.5; % in secs
    prop_win = win_edge_durn*2/((length(s.v)-1)/s.fs);
    tw = tukeywin(length(s.v),prop_win); 
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v.*tw);
    s_filt.v = s.v-s_filt.v;
else
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
s_filt.t = s.t;
end

function ekg = load_and_process_ecg(data, s, up)

%% Load ECG
ekg = data(s).ekg;
ekg.t = [0:(length(ekg.v)-1)]/ekg.fs;

%% Process ECG
filt_characteristics = up.elim_low_freqs;
ekg = elim_vlfs(ekg, filt_characteristics);

end
