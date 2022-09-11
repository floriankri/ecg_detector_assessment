function [QRSindex, rr, multileadVA, ecgLength] = Adapted_simPAF_gen_multilead_VA(ecgLength, targets_SR_AF, rr, AFburden, realVAon, realAAon, DATApqrst)
% 
% [] = simPAF_gen_multilead_VA() returns multilead (15 lead) ventricular
% activity. A set of 100 15-lead ECGs with SR selected from the PTB Diagnostic
% ECG Database is used as a basis for modeling ventricular activity. The ECGs 
% of the PTB database are first subjected to baseline removal and QRST delineation. 
% The original T-waves are then resampled to a fixed width and, depending on 
% the type of rhythm, width-adjusted to match prevailing heart rate. Since 
% the original ECGs last just for about 2 min, QRST complexes are subjected
% to repeated concatenation until desired length of ECG is obtained. The TQ 
% interval is interpolated using a cubic spline interpolation.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% Generated leads:
% multileadVA(1,:) - I        multileadVA(7,:) - V1      multileadVA(13,:) - X     
% multileadVA(2,:) - II       multileadVA(8,:) - V2      multileadVA(14,:) - Y  
% multileadVA(3,:) - III      multileadVA(9,:) - V3      multileadVA(15,:) - Z  
% multileadVA(4,:) - aVR      multileadVA(10,:) - V4 
% multileadVA(5,:) - aVL      multileadVA(11,:) - V5 
% multileadVA(6,:) - aVF      multileadVA(12,:) - V6 

rr = round(rr*1000);    % rr intervals in miliseconds

% Generate PQRST complexes
switch realVAon 
    case 0 % Generate synthetic PQRST complexes
        [data.pqrst, data.Qind, data.Rind, data.Sind, data.Tind] = simPAF_gen_syn_VA(100);  
    case 1 % Load random patch of PQRST complexes extracted from PTB database
        sigNum = randi([1 100]);
        %load('DATA_PQRST_real')
        data.pqrst = DATApqrst(sigNum).pqrst;
        data.Qind = DATApqrst(sigNum).Qind;
        data.Rind = DATApqrst(sigNum).Rind;
        data.Sind = DATApqrst(sigNum).Sind;
        data.Tind = DATApqrst(sigNum).Tind;
end

% Original PQRST is subjected to repeated concatenation until required
% number of beats is obtained
arraySize = length(data.pqrst(1,:,1));
pqrstLength = length(data.pqrst(1,1,:));
pqrst = zeros(15, ecgLength,pqrstLength);

j = 1;
for i = 1:ecgLength+1 % Additonal is required to meat exact number of RR intervals
    pqrst(:,i,:) = data.pqrst(:,j,:);
    j = j + 1;
    if j > arraySize
        j = 1; 
    end
end

%Load QRST complexes for all leads
cnt = 0;
for lead = 1:15
    ecgSig = [];
    AFend = 0;
    cnt = cnt + 1;
    switch realAAon
        case 0 % Synthetic atrial activity
            rIndex = data.Rind - data.Qind;
            Qind = data.Qind; Rind = data.Rind - data.Qind; 
            Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
            for beatNr = 1:ecgLength
                QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                QRSTtempNext = pqrst(lead,beatNr+1,Qind+1:end);
                QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                % Divide PQRST into two parts 
                QRST_PS = QRSTtemp(1,1:Sind); % The PQRS part 
                QRST_ST = QRSTtemp(1,Sind+1:Tind); % The ST part 
                % Resample T wave according to current RR interval
                % Signal is prolonged in order to avoid boundary effects of resampling
                QRST_ST  = [QRST_ST(1,1)*ones(1,20) QRST_ST QRST_ST(1,end)*ones(1,20)];
        
                % Gradual prolongation of T wave after AF is terminated
                if targets_SR_AF(beatNr) == 0
                    if AFend == 1
                        ST_length = round(length(QRST_ST)*sqrt((0.3+m*0.1)*rr(beatNr)/1000));
                        m = m + 1;
                        if m == 7
                            AFend = 0;
                        end
                    else
                        ST_length = round(length(QRST_ST)*sqrt(rr(beatNr)/1000));
                    end
                else
                    ST_length = round(length(QRST_ST)*sqrt(0.35));
                    AFend = 1;
                    m = 1;
                end
                
                QRST_ST = resample(QRST_ST, ST_length+40, length(QRST_ST));  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR)
                QRSTc = [QRST_PS QRST_ST(21:end-20)]; % Corrected PQRST
                % Find TQ length
                TQlength = rr(beatNr) - length(QRSTc);
                % Protect against the error of to low heart rate
                if  TQlength < 2
                    TQlength = 2;
                end
                % Interpolate TQ interval
                x = [1 TQlength];
                xi = 1:1:TQlength;
                y = [QRSTc(1,end) QRSTtempNext(1, 1)];
                TQ = interp1(x,y,xi,'linear');
                ecgSig = [ecgSig QRSTc TQ];
                if lead == 15
                    rIndex = [rIndex (rIndex(1,end) + TQlength + length(QRSTc))];
                end
            end
        case 1 % Real atrial activity
            if AFburden == 1 % Entire signal is AF
                rIndex = data.Rind - data.Qind;
            else
                rIndex = data.Rind;
            end
            for beatNr = 1:ecgLength
                if beatNr == 1
                    if AFburden == 1
                        Qind = data.Qind; Rind = data.Rind - data.Qind; 
                        Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;
                        QRSTtemp = pqrst(lead,beatNr,Qind + 1:end);
                        QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                    else
                        Qind = data.Qind; Rind = data.Rind;
                        Sind = data.Sind; Tind = data.Tind;
                        QRSTtemp = pqrst(lead,beatNr,1:end);
                        QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                    end
                else
                    if targets_SR_AF(beatNr) == 0
                        if targets_SR_AF(beatNr-1) == 0
                            QRSTtemp = pqrst(lead,beatNr,1:end);
                            QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind; 
                            Sind = data.Sind; Tind = data.Tind;
                        else % targets_SR_AF(beatNr-1) == 1
                            QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                            QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind - data.Qind; 
                            Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;  
                        end
                    else % targets_SR_AF(beatNr) == 1
                        if targets_SR_AF(beatNr-1) == 1
                            QRSTtemp = pqrst(lead,beatNr,Qind+1:end);
                            QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind - data.Qind; 
                            Sind = data.Sind - data.Qind; Tind = data.Tind - data.Qind;      
                        else % targets_SR_AF(beatNr-1) == 0
                            QRSTtemp = pqrst(lead,beatNr,1:end);
                            QRSTtemp = simPAF_correct_baseline(QRSTtemp);
                            Qind = data.Qind; Rind = data.Rind; 
                            Sind = data.Sind; Tind = data.Tind;  
                        end
                    end
                end
                
                if beatNr == length(targets_SR_AF)
                    if targets_SR_AF(beatNr) == 0
                        RindNext = data.Rind; 
                        QRSTtempNext = pqrst(lead,beatNr,1:end);
                        QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                    else
                        RindNext = data.Rind-data.Qind;
                        QRSTtempNext = pqrst(lead,beatNr,Qind+1:end);
                        QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                    end
                else
                    if targets_SR_AF(beatNr+1) == 0
                        if targets_SR_AF(beatNr) == 1
                            RindNext = data.Rind-data.Qind;
                            QRSTtempNext = pqrst(lead,beatNr+1,Qind+1:end);
                            QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                        else
                            RindNext = data.Rind;
                            QRSTtempNext = pqrst(lead,beatNr+1,1:end);
                            QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                        end
                    else
                        if targets_SR_AF(beatNr) == 0
                            RindNext = data.Rind;
                            QRSTtempNext = pqrst(lead,beatNr+1,1:end);
                            QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                        else
                            RindNext = data.Rind-data.Qind;
                            QRSTtempNext = pqrst(lead,beatNr+1,Qind+1:end);
                            QRSTtempNext = simPAF_correct_baseline(QRSTtempNext);
                        end
                    end
                end
                % Divide PQRST into two parts
                QRST_PS = QRSTtemp(1,1:Sind); % The PQRS part 
                QRST_ST = QRSTtemp(1,Sind+1:Tind); % The ST part 
                % Resample T wave according to current RR interval
                % Signal is prolonged in order to avoid boundary effects of resampling            
                QRST_ST = [QRST_ST(1,1)*ones(1,20) QRST_ST QRST_ST(1,end)*ones(1,20)];
                % Gradual prolongation of T wave after AF is terminated
                if targets_SR_AF(beatNr) == 0
                    if AFend == 1
                        ST_length = round(length(QRST_ST)*sqrt((0.3+m*0.1)*rr(beatNr)/1000));
                        m = m + 1;
                        if m == 7
                            AFend = 0;
                        end
                    else
                        ST_length = round(length(QRST_ST)*sqrt(rr(beatNr)/1000));
                    end
                else
                    ST_length = round(length(QRST_ST)*sqrt(0.35));
                    AFend = 1;
                    m = 1;
                end               
                QRST_ST = resample(QRST_ST, ST_length+40, length(QRST_ST));  % T wave length correction QT = QTc*sqrt(RR) ~ T = Tc*sqrt(RR)
                QRSTc = [QRST_PS QRST_ST(21:end-20)]; % Corrected PQRST
                % Find TQ length
                TQlength = rr(beatNr) - length(QRSTc) + Rind - RindNext;
                % Protect against the error of to low heart rate
                if  TQlength < 2
                    TQlength = 2;
                end
                % Interpolate TQ interval
                x = [1 TQlength];
                xi = 1:1:TQlength;
                y = [QRSTc(1,end) QRSTtempNext(1, 1)];
                TQ = interp1(x,y,xi,'linear');
                ecgSig = [ecgSig QRSTc TQ];             
                if lead == 15
                    rIndex = [rIndex (rIndex(1,end) + TQlength + length(QRSTc) - Rind + RindNext)];
                end
            end
    if cnt == 1        
    multileadVA = zeros(15, length(ecgSig));
    end
    end
    if lead == 15
        rr = diff(rIndex);
        QRSindex = rIndex(1:end-1);  
    end
    multileadVA(lead, :) = ecgSig;
end
ecgLength = length(ecgSig);
end