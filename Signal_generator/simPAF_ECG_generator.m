function [simPAFdata, initialParameters] = simPAF_ECG_generator(rrLength, realRRon, realVAon, realAAon, AFburden, stayInAF, APBrate, noiseType, noiseRMS, onlyRR)
% [] = simPAF_ECG_generator() returns a 15-by-N matrix containing 15 lead
% ECGs. Three types of ECG signals can be generated: sinus rhythm (SR) (AF burden set to 0, 
% AF (AF burden set to 1) or PAF (AF burden any value from the interval 
% (0, 1)). Standard leads I, II, III, aVR, aVL, aVF, V1, V2, V3, V4, V5, V6
% and Frank leads X, Y, Z are generated (sampling frequency 1000 Hz).
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% Generated leads:
% multileadVA(1,:) - I      multileadVA(7,:) - V1    multileadVA(13,:) - X     
% multileadVA(2,:) - II     multileadVA(8,:) - V2    multileadVA(14,:) - Y  
% multileadVA(3,:) - III    multileadVA(9,:) - V3    multileadVA(15,:) - Z  
% multileadVA(4,:) - aVR    multileadVA(10,:) - V4 
% multileadVA(5,:) - aVL    multileadVA(11,:) - V5 
% multileadVA(6,:) - aVF    multileadVA(12,:) - V6 
%
% Input arguments:
% rrLength indicates the length of the desired ECG signal (in RR intervals)
%
% realRRon 1 indicates that real RR intervals are used, 0 - synthetic
% realVAon 1 indicates that real ventricular activity is used, 0 - synthetic
% realAAon 1 indicates that real atrial activity is used, 0 - synthetic
%
% onlyRR 1 - only RR intervals are generated, 0 - multilead ECG is generated
%
% AFburden is a value between 0 and 1. 0 - the entire signal is SR, 
% 1 - the entire signal is AF.
%
% stayInAF denotes the probability to stay in AF. 
%
% APBrate is any value from 0 to 0.5, where 0 implies that no APBs occur, 
% whereas 0.5 means that 50 % of all beats are APBs. Note that P wave
% morphology is altered only when synthetic atrial activiti is used 
% (realAAon = 0).
%
% noiseType: a number from 0 to 4
% 0 - no noise added (noise RMS = 0 mV)
% 1 - motion artefacts
% 2 - electrode movement artefacts
% 3 - baseline wander
% 4 - mixture of type 1, type 2 and type 3 noises
%
% noiseLevel - noise level in milivolts, i.e. 0.02 corresponds to 20 uV
%
% Output arguments:
% simPAFdata returns generated data (multilead ECG, multilead ventricular
% activity, multilead atrial activity, QRS index, etc.). initialParameters 
% returns initial parameter values used to generated ECG signals. 
%
% Known problems:
% *AV node model used for generating RR intervals during AF is relatively slow.
% *Synthetic P wave amplitude is nearly 1.5 lower in several leads than that
% observed in reality (at least for healthy patients). Parameters for 
% simulating Type 2 P waves are taken from the paper by Havmoller et al. 
% Age-related changes in P wave morphology in healthy subjects. 
% BMC Cardiovascular Disorders, 7(1), 22, 2007.
% *Interpolated TQ intervals (using a cubic spline interpolation) sometimes 
% do not look realistic.
%

switch onlyRR
    case 1 % only RR intervals are generated
        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simPAF_gen_initial_parameters;
        % Generate RR intervals
        [rr, targets_SR_AF, targets_APB, pafBoundaries, pafEpisLength] = simPAF_gen_rr_intervals(rrLength, fibFreqz, APBrate, realRRon, stayInAF, AFburden);
        simPAFdata.rr = rr;
        simPAFdata.multileadECG = [];
        simPAFdata.multileadVA = [];
        simPAFdata.multileadAA = [];
        simPAFdata.multileadNoise = [];
        simPAFdata.QRSindex = [];
        simPAFdata.targets_SR_AF = targets_SR_AF;
        simPAFdata.targets_APB = targets_APB;
        simPAFdata.pafBoundaries = pafBoundaries;
        simPAFdata.pafEpisLength = pafEpisLength;
        simPAFdata.ecgLength = [];

        initialParameters.fibFreqz = fibFreqz;
        initialParameters.rrLength = rrLength;
        initialParameters.APBrate = APBrate;
        initialParameters.realRRon = realRRon;
        initialParameters.realVAon = realVAon;
        initialParameters.realVAon = realAAon;
        initialParameters.stayInAF = stayInAF;
        initialParameters.AFburden = AFburden;
        initialParameters.noiseType = noiseType;
        initialParameters.noiseRMS = noiseRMS;
    case 0 % multilead ECG is generated
        % Check for errors:
        if (realVAon == 0) && (realAAon == 1)
            msg = ('Selection of synthetic ventricular activity and real atrial activity is not allowed');
            error('MyComponent:incorrectType', msg);
        end

        % Generate initial parameters (fibrillatory frequency)
        fibFreqz = simPAF_gen_initial_parameters;  
        % Generate RR intervals
        [rrIn, targets_SR_AF, targets_APB, pafBoundaries, pafEpisLength] = simPAF_gen_rr_intervals(rrLength, fibFreqz, APBrate, realRRon, stayInAF, AFburden);
        % Generate multilead ventricular activity
        [QRSindex, rr, multileadVA, ecgLength] = simPAF_gen_multilead_VA(rrLength, targets_SR_AF, rrIn, AFburden, realVAon, realAAon);
        % Generate multilead atrial activity
        multileadAA = simPAF_gen_multilead_AA(targets_SR_AF, targets_APB, QRSindex, fibFreqz, realAAon, ecgLength, AFburden);
        % Generate multilead noise
        multileadNoise = simPAF_gen_noise(ecgLength, noiseType, noiseRMS);
        % Generate multilead noise
        multileadECG = multileadVA + multileadAA + multileadNoise;

        simPAFdata.rr = rr;
        simPAFdata.multileadECG = multileadECG;
        simPAFdata.multileadVA = multileadVA;
        simPAFdata.multileadAA = multileadAA;
        simPAFdata.multileadNoise = multileadNoise;
        simPAFdata.QRSindex = QRSindex;
        simPAFdata.targets_SR_AF = targets_SR_AF;
        simPAFdata.targets_APB = targets_APB;
        simPAFdata.pafBoundaries = pafBoundaries;
        simPAFdata.pafEpisLength = pafEpisLength;
        simPAFdata.ecgLength = ecgLength;

        initialParameters.fibFreqz = fibFreqz;
        initialParameters.rrLength = rrLength;
        initialParameters.APBrate = APBrate;
        initialParameters.realRRon = realRRon;
        initialParameters.realVAon = realVAon;
        initialParameters.realVAon = realAAon;
        initialParameters.stayInAF = stayInAF;
        initialParameters.AFburden = AFburden;
        initialParameters.noiseType = noiseType;
        initialParameters.noiseRMS = noiseRMS;          
end

end