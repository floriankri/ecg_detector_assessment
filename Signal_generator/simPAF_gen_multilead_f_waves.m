function multilead_f_waves = simPAF_gen_multilead_f_waves(fibFreqz, fWaveLength)
% multilead_f_waves = simPAF_gen_multilead_f_waves() returns multilead (15 lead) 
% simulated f-waves generated using an extended sawtooth f wave model. For
% details see a paper by Petrenas et al. An echo state neural network for
% QRST cancellation during atrial fibrillation. Transactions on Biomedical
% Engineering, 59(10), 2012. f-waves model is further extended by adding low
% frequency and high frequency noise components.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% Generated leads:
% multilead_f_waves(1,:) - I        multilead_f_waves(7,:) - V1           
% multilead_f_waves(2,:) - II       multilead_f_waves(8,:) - V2        
% multilead_f_waves(3,:) - III      multilead_f_waves(9,:) - V3        
% multilead_f_waves(4,:) - aVR      multilead_f_waves(10,:) - V4 
% multilead_f_waves(5,:) - aVL      multilead_f_waves(11,:) - V5 
% multilead_f_waves(6,:) - aVF      multilead_f_waves(12,:) - V6 
%
% multilead_f_waves(13,:) - X
% multilead_f_waves(14,:) - Y
% multilead_f_waves(15,:) - Z

% Guillem matrix (optimized for reconstructing 12 lead ECG from vcg)
Z = [-0.213276421450882 0.210438180414305 -1.32199603270407;
    0.527759188353754 0.0863232920392896 -1.16070747844017;
    0.752551447000432 0.315418742662643 -1.19110821365617;
    0.926632574756099 0.418704634911852 -1.00757071040099;
    0.968764605349974 0.401332344956578 -0.555083558926398;
    0.937743573516388 0.408711021211355 -0.317655185251036;
    0.987645105751296 0.0158866121670719 0.277151634261775;
    0.0628391410856677 1.29202360737584 0.200506800589521];

% Diagonal matrix 
coef = [1.5 0 0 0 0 0 0 0;
        0 1.2 0 0 0 0 0 0;
        0 0 0.8 0 0 0 0 0;
        0 0 0 0.5 0 0 0 0;
        0 0 0 0 0.4 0 0 0;
        0 0 0 0 0 0.3 0 0;
        0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 1];

    
fWavesRMS = simPAF_gen_rand_num(15, 45);
vcg_f_waves(1,:) = simPAF_gen_single_lead_f_waves(fWaveLength, fibFreqz-0.05*fibFreqz, fWavesRMS);

fWavesRMS = simPAF_gen_rand_num(15, 40);
vcg_f_waves(2,:) = -1*simPAF_gen_single_lead_f_waves(fWaveLength, fibFreqz, fWavesRMS);

fWavesRMS = simPAF_gen_rand_num(25, 70);
vcg_f_waves(3,:) = simPAF_gen_single_lead_f_waves(fWaveLength, fibFreqz+0.05*fibFreqz, fWavesRMS);

multilead_f_waves= coef*Z*vcg_f_waves;
multilead_f_waves(9,:) = multilead_f_waves(8,:) - multilead_f_waves(7,:);     % III
multilead_f_waves(10,:) = -(multilead_f_waves(7,:)+multilead_f_waves(8,:))/2; % aVR
multilead_f_waves(11,:) =  multilead_f_waves(7,:)-multilead_f_waves(8,:)/2;   % aVL
multilead_f_waves(12,:) =  multilead_f_waves(8,:)-multilead_f_waves(7,:)/2;   % aVF

multilead_f_waves(13,:) = vcg_f_waves(1,:);  % X
multilead_f_waves(14,:) = vcg_f_waves(2,:);  % Y
multilead_f_waves(15,:) = vcg_f_waves(3,:);  % Z

multilead_f_waves_temp = multilead_f_waves;
multilead_f_waves(1:6,:,:) = multilead_f_waves(7:12,:,:);
multilead_f_waves(7:12,:,:) = multilead_f_waves_temp(1:6,:,:);
  
end