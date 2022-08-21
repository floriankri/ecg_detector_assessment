function multileadAA = simPAF_gen_multilead_AA(targets_SR_AF, targets_APB, QRSindex, fibFreqz, realAAon, ecgLength, AFburden)
% multileadAA = simPAF_gen_multilead_AA() returns multilead (15 lead) atrial
% activity. 
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% Generated leads:
% multileadAA(1,:) - I      multileadAA(7,:) - V1    multileadAA(13,:) - X     
% multileadAA(2,:) - II     multileadAA(8,:) - V2    multileadAA(14,:) - Y  
% multileadAA(3,:) - III    multileadAA(9,:) - V3    multileadAA(15,:) - Z  
% multileadAA(4,:) - aVR    multileadAA(10,:) - V4 
% multileadAA(5,:) - aVL    multileadAA(11,:) - V5 
% multileadAA(6,:) - aVF    multileadAA(12,:) - V6 

% Decide which rhythm type to generate
if AFburden == 0
    rhythmType = 0; % SR
elseif AFburden == 1
    rhythmType = 1; % AF
elseif AFburden > 0 && AFburden < 1
    rhythmType = 2; % PAF
else
    error('AF burden must be a value between 0 and 1')
end

multileadAA = zeros(15, ecgLength);
% Generate atrial activity
switch rhythmType
    case 0 % Entire signal is SR
       if realAAon == 0
           Nrr = length(targets_SR_AF);
           P_waves = simPAF_gen_multilead_P_waves(Nrr); 
           
           % Replace P waves during ectopic beats
           numABPs = length(find(targets_APB == 1));
           if numABPs > 0
               P_waves_APBs = simPAF_gen_multilead_P_waves(numABPs);
               iAPB = 1;
               for n = 1:Nrr 
                   if targets_APB(n) == 1
                       P_waves(:,n+1,:) = P_waves_APBs(:,iAPB,:);
                       iAPB = iAPB + 1;
                   end
               end
           end
           % Insert P waves
           for p_num = 2:Nrr
               multileadAA(:, QRSindex(p_num)-249:QRSindex(p_num)-100) = P_waves(:,p_num,:);             
           end
       end      
    case 1 % Entire signal is AF
        if realAAon == 0
            multileadAA = simPAF_gen_multilead_f_waves(fibFreqz, ecgLength);
        else
            sigNum = randi([1 20]);
            load('DATA_f_waves_real')
            f_waves = DATAfWaves(sigNum).f_waves;
            
            clear DATA_f_waves_real
            
            if ecgLength > length(f_waves)   % Concatenate if RR series is shorter than the desired length
                nCycles = ceil(ecgLength/length(f_waves));
                multileadAA = [];
                for i = 1:nCycles
                    multileadAA = [multileadAA f_waves];
                end
                multileadAA = multileadAA(:,1:ecgLength);  
            else
                fStart = randi([1 (length(f_waves)-ecgLength)]); % Randomly select start position
                multileadAA = f_waves(:,fStart+1:fStart+ecgLength);
            end   
        end
                
    case 2 % PAF
        Nrr = length(targets_SR_AF);
        if realAAon == 0  % Synthetic atrial activity is prefered
           P_waves = simPAF_gen_multilead_P_waves(Nrr); 
           % Replace P waves during ectopic beats
           numABPs = length(find(targets_APB == 1));
           if numABPs > 0
               P_waves_APBs = simPAF_gen_multilead_P_waves(numABPs);
               iAPB = 1;
               for n = 1:Nrr 
                   if targets_APB(n) == 1
                       P_waves(:,n+1,:) = P_waves_APBs(:,iAPB,:);
                       iAPB = iAPB + 1;
                   end
               end
           end
           % Generate f waves
           f_waves = simPAF_gen_multilead_f_waves(fibFreqz, ecgLength);
           % Insert P waves
           for p_num = 2:Nrr
               multileadAA(:, QRSindex(p_num)-249:QRSindex(p_num)-100) = P_waves(:,p_num,:);  
           end
           % Insert f waves
           for p_num = 2:Nrr
               if targets_SR_AF(p_num) == 1
                   praIndex = QRSindex(p_num);
                   while targets_SR_AF(p_num) > 0 && p_num < Nrr
                      p_num = p_num + 1;
                   end 
                   pabIndex = QRSindex(p_num);
                   multileadAA(:, praIndex:pabIndex) = f_waves(:,praIndex:pabIndex);   
               end
           end
               
        else % Real atrial activity is prefered
            sigNum = randi([1 20]);
            load('DATA_f_waves_real')
            f_waves = DATAfWaves(sigNum).f_waves;
            
            clear DATA_f_waves_real
            
            if ecgLength > length(f_waves)   % Concatenate if RR series is shorter than the desired length
                nCycles = ceil(ecgLength/length(f_waves));
                multileadfWaves = [];
                for i = 1:nCycles
                    multileadfWaves = [multileadfWaves f_waves];
                end
                multileadfWaves = multileadfWaves(:,1:ecgLength);  
            else
                fStart = randi([1 (length(f_waves)-ecgLength)]); % Randomly select start position
                multileadfWaves = f_waves(:,fStart+1:fStart+ecgLength);
            end 
            for p_num = 2:Nrr
                if targets_SR_AF(p_num) == 1
                    praIndex = QRSindex(p_num);
                    while targets_SR_AF(p_num) > 0 && p_num < Nrr
                        p_num = p_num + 1;
                    end 
                    pabIndex = QRSindex(p_num);
                    multileadAA(:, praIndex:pabIndex) = multileadfWaves(:,praIndex:pabIndex); 
                end
            end
            
        end
             
end

