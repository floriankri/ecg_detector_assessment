function [rr, targets_SR_AF, targets_APB, pafBoundaries, pafEpisodeLength] = Adapted_simPAF_gen_rr_intervals(rrLength, fibFreqz, APBrate, realRRon, stayInAF, AFburden, DATArrSR, DATArrAF)
%
% [] = simPAF_gen_rr_intervals() returns RR intervals during PAF. 
% The dominant length of PAF episodes, and the percentage of the total time AF is
% present (AF burden) is determined by a first-order two-state Markov chain.
%
% An atrioventricular node model in which the ventricles are assumed to be 
% activated by the arriving atrial impulses according to a Poisson process
% is used to generate realistic RR series during an episode of AF (Corino et al 2011)
%
% Ventricular rhythm during SR is simulated using RR interval generator 
% proposed by McSharry et al(2003) in which both the impact of 
% parasympathetic stimulation (respiratory sinus arrhythmia) and baroreflex 
% regulation (Mayer waves) is modeled by a bimodal power spectrum
%
% In case of real data the entire MIT-BIH Normal Sinus Rhythm database, 
% consisting of 18 long-term ECG recordings, was taken as a basis for the set
% of 18 RR interval series of sinus rhythm.
% 
% The Long Term Atrial Fibrillation database was used to compose a set
% of AF rhythm. In total, 69 RR interval series were extracted out of 
% 84 long-term ECG recordings. The remaining records were excluded since 
% they did not meet the criteria of <5000 RR intervals of AF.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%

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

switch rhythmType  % 0 - sinus rhythm, 1 - AF, 2 - PAF
    case 0 % The entire rhythm is SR
        if realRRon == 1 % Use real RR series
            rr = Adapted_simPAF_construct_real_RR(0, rrLength, DATArrSR, DATArrAF);
        else % Use simulated RR series
            rr = simPAF_gen_SR_RR_intervals(rrLength);
            rr = rr(1:rrLength);
        end
        targets_SR_AF = zeros(1,rrLength);
        pafBoundaries(1,1) = 0; 
        pafBoundaries(1,2) = 0;
        pafEpisodeLength = 0;
        
    case 1 % The entire rhythm is AF
        if realRRon == 1 % Use real RR series
            rr = Adapted_simPAF_construct_real_RR(1, rrLength, DATArrSR, DATArrAF);
        else
            rr = simPAF_gen_AF_RR_intervals(fibFreqz,rrLength);
            rr = rr(1:rrLength); 
        end
        targets_SR_AF = ones(1,rrLength);
        pafBoundaries(1,1) = 1; 
        pafBoundaries(1,2) = rrLength;
        pafEpisodeLength = rrLength;
        
    case 2 % PAF
        % Generate alternating SR and AF episodes
        goToNS = 1-stayInAF;
        goToAF = (goToNS*AFburden)/(1-AFburden);
        stayInNS = 1 - goToAF;

        trans = [stayInNS,1-stayInNS;1-stayInAF, stayInAF];
        [~,targets_SR_AF] = hmmgenerate(rrLength,trans,[0; 1]);

        % Find the number of AF episodes and the length of each episode
        k = 1;
        pafEpisodeLength = [];
        for p = 1:length(targets_SR_AF)-1
            if targets_SR_AF(p) == 2
                if targets_SR_AF(p+1) == 2
                    k = k + 1;
                    if p == (length(targets_SR_AF)-1)
                        pafEpisodeLength = [pafEpisodeLength k]; 
                    end
                else
                pafEpisodeLength = [pafEpisodeLength k];
                k = 1;
                end         
            end
        end
        % Find boundaries of each PAF episode
        if rrLength == length(find(targets_SR_AF == 1))% The entire signal is SR   
            pafBoundaries(1,1) = 0; 
            pafBoundaries(1,2) = 0;
        elseif rrLength == length(find(targets_SR_AF == 2))% The entire signal is AF
            pafBoundaries(1,1) = 1; 
            pafBoundaries(1,2) = rrLength;
        elseif rrLength ~= length(find(targets_SR_AF == 2))% The signal with PAF
            diffTar = diff(targets_SR_AF);
            j = 1;
            k = 1;
            flag = 1;
            for i = 1:length(diffTar)
                if diffTar(i) == 1
                    pafBoundaries(j,1) = i + 1;
                    j = j + 1;
                    flag = 1;
                end
                if diffTar(i) == -1
                    pafBoundaries(k,2) = i;
                    k = k + 1;
                    flag = 2;
                end   
                if i == length(diffTar)
                    if flag == 1
                        pafBoundaries(k,2) = i+1;
                    end
                end
            end
        end
        
        % Generate RR intervals during AF
        rrLengthAF = length(find(targets_SR_AF == 2));
        if realRRon == 1 % Use real RR series
            rrAF = simPAF_construct_real_RR(1, rrLengthAF);
        else
            rrAF = simPAF_gen_AF_RR_intervals(fibFreqz,rrLengthAF);
            rrAF = rrAF(1:rrLengthAF); 
        end
        
        % Generate RR intervals during SR
        rrLengthSR = length(find(targets_SR_AF == 1));
        if realRRon == 1 % Use real RR series
            rrSR = simPAF_construct_real_RR(0, rrLengthSR);
        else % Use simulated RR series
            rrSR = simPAF_gen_SR_RR_intervals(rrLengthSR);
            rrSR = rrSR(1:rrLengthSR);
        end
        % Make sure that the average RR value during SR is larger than that in AF
        if mean(rrSR) < mean(rrAF)
            rrSR = rrSR + (mean(rrAF) - mean(rrSR));
        end
        
        % Construct PAF RR series
        j = 1;
        k = 1;
        for i = 1:rrLength
            if targets_SR_AF(i) == 1
                rr(i) = rrSR(j); 
                j = j+1;
            else
                rr(i) = rrAF(k); 
                k = k+1;
            end
        end
        targets_SR_AF = targets_SR_AF - 1;
end

%% Insert atrial premature beats
if APBrate > 0.5
    error('APB rate must be less or equal to 0.5')
elseif APBrate > 0
    targets_APB = zeros(1,rrLength);
    transAPB = [1-APBrate,APBrate;0.99999,0.00001];
    [~,statesAPB] = hmmgenerate(rrLength,transAPB, [0; 1]);
    % Remove APBs during AF
    statesAPB = statesAPB.*(abs(targets_SR_AF-1));
    rrAPB=rr;
    
    apbType = randi([0 3],1,1);

    switch apbType 
        case 0 % APBs with sinus reset
            for i=1:length(statesAPB)-1
                if statesAPB (i) == 2
                    rrAPB(i)=0.8*rrAPB(i);   
                    targets_APB(i) = 1;
                end
            end
        case 1 % Interpolated PACs
            for i=1:length(statesAPB)-1
                if statesAPB (i) == 2
                    rrAPBtemp = rrAPB(i);
                    rrAPB(i)=0.6*rrAPB(i);
                    rrAPB(i+1)=rrAPBtemp - rrAPB(i);   
                    targets_APB(i) = 1;
                end
            end   
        case 2 % PACs with delayed sinus reset
            for i=1:length(statesAPB)-1
                if statesAPB (i) == 2
                    rrAPBtemp = rrAPB(i);
                    rrAPB(i)= 0.8*rrAPB(i);
                    rrAPB(i+1)= 1.2*rrAPB(i+1);   
                    targets_APB(i) = 1;
                end
            end
            case 3 %PACs with full compensatory pause
                for i=1:length(statesAPB)-1
                    if statesAPB (i) == 2
                        rrAPBtemp = rrAPB(i);
                        rrAPB(i)= 0.8*rrAPB(i);
                        rrAPB(i+1)= 2*rrAPBtemp-rrAPB(i);   
                        targets_APB(i) = 1;
                    end
                end            
    end
    rr = rrAPB;
else
    targets_APB = zeros(1,rrLength);
end


% %% Switching between AF an SR
% goToNS = 1-stayInAF;
% goToAF = (goToNS*AFburden)/(1-AFburden);
% stayInNS = 1 - goToAF;
% 
% trans = [stayInNS,1-stayInNS;1-stayInAF, stayInAF];
% [~,states_SR_AF] = hmmgenerate(ecgLength,trans,[0; 1]);
% 
% %% Calculate the number of AF episodes generate by Markov chain
% k = 1;
% episodeLength = [];
% for p = 1:length(states_SR_AF)-1
%     if states_SR_AF(p) == 2
%         if states_SR_AF(p+1) == 2
%             k = k + 1;
%            if p == (length(states_SR_AF)-1)
%               episodeLength = [episodeLength k]; 
%            end
%         else
%             episodeLength = [episodeLength k]; 
%             k = 1;
%         end         
%     end
% end
% 
% %% Find boundaries of each PAF episode
% if ecgLength == length(find(states_SR_AF == 1))
%     % The entire signal is SR
%     beginOfAF(1,1) = 0; 
%     beginOfAF(1,2) = 0;
% elseif ecgLength == length(find(states_SR_AF == 2))
%     % The entire signal is AF
%     beginOfAF(1,1) = 1; 
%     beginOfAF(1,2) = ecgLength;
% elseif ecgLength ~= length(find(states_SR_AF == 2))
%     diffTar = diff(states_SR_AF);
%     j = 1;
%     k = 1;
%     flag = 1;
%     for i = 1:length(diffTar)
%         if diffTar(i) == 1
%             beginOfAF(j,1) = i + 1;
%             j = j + 1;
%             flag = 1;
%         end
%         if diffTar(i) == -1
%             beginOfAF(k,2) = i;
%             k = k + 1;
%             flag = 2;
%         end
%     
%         if i == length(diffTar)
%             if flag == 1
%             beginOfAF(k,2) = i+1;
%             end
%         end
%     end
% end
% 
% %% 
%    
% 
end