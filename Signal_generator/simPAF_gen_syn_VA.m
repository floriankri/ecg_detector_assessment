function [pqrst, Qind, Rind, Sind, Tind] = simPAF_gen_syn_VA(num)
%
% [] = simPAF_gen_syn_VA() returns simulated QRST complexes. 
% Ventricular activity is simulated using a dynamic ECG model proposed in 
% Sameni et al. 2007, which is an extended version of a single channel
% ECG simulator originally proposed in McSharry et al. 2003. By using this
% model, three orthogonal Frank leads of N samples length are generated.
% Ventricular activity is generated in 3 Frank leads and then 12-lead signals
% are obtained by applying a transformation matrix.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
N = 1000;                       % Signal length in samples 
fs = 1000;                      % sampling rate
% num - number of PQRST compelxes
pqrst = zeros(15,num,560);
pqrstFrank = zeros(3,560);

% Dipole parameters 
F = 1;                          % heart rate 
R0 = simPAF_Rotate3D(0,0,0);    % dipole rotation matrices (tetax,tetay,tetaz) 
Lambda = eye(3); 
teta0 = -pi/2;                  % initial phase of the ECG 
 
%% QRS
Qw = simPAF_gen_rand_num(0.05, 0.08);
Rw = simPAF_gen_rand_num(0.05, 0.08);
Sw = simPAF_gen_rand_num(0.05, 0.08);

alphaiQRS.x = [ simPAF_gen_rand_num(-0.4, -0.05)  simPAF_gen_rand_num(0.4, 1.5)     0]; 
biQRS.x     = [Qw    Rw    Sw]; 
tetaiQRS.x  = [-0.1    0    0.1]; 

alphaiQRS.y = [ 0    simPAF_gen_rand_num(0.1, 0.7)   simPAF_gen_rand_num(-0.3, -0.05)]; 
biQRS.y     = [Qw   Rw   Sw]; 
tetaiQRS.y  = [-0.1   0   0.1]; 
  
alphaiQRS.z = [simPAF_gen_rand_num(-0.4, -0.05)   0  simPAF_gen_rand_num(0.1, 1)]; 
biQRS.z     = [Qw   Rw  Sw]; 
tetaiQRS.z  = [ -0.1   0   0.1];

%% T wave
Tw = simPAF_gen_rand_num(0.5, 0.7);

Txa = simPAF_gen_rand_num(0.02, 0.08);
Tya = simPAF_gen_rand_num(0.01, 0.03);
Tza = simPAF_gen_rand_num(-0.02, -0.06);

alphaiT.x = [Txa   2*Txa   3*Txa]; 
biT.x     = [Tw     Tw/2    Tw/4]; 
tetaiT.x  = [1.1     1.4    1.6];

alphaiT.y = [Tya   2*Tya   3*Tya]; 
biT.y     = [Tw     Tw/2    Tw/4]; 
tetaiT.y  = [1.1     1.4    1.6];

alphaiT.z = [Tza   2*Tza   3*Tza]; 
biT.z     = [Tw     Tw/2    Tw/4]; 
tetaiT.z  = [1.1     1.2    1.4];
%%
alphai.x = [alphaiQRS.x alphaiT.x];
alphai.y = [alphaiQRS.y alphaiT.y];
alphai.z = [alphaiQRS.z alphaiT.z];

bi.x = [biQRS.x biT.x];
bi.y = [biQRS.y biT.y];
bi.z = [biQRS.z biT.z];

tetai.x = [tetaiQRS.x tetaiT.x];
tetai.y = [tetaiQRS.y tetaiT.y];
tetai.z = [tetaiQRS.z tetaiT.z];

f = randi([10 30])/10;
% Generate phase for each variability related frequency
if rand(1,1) > 0.5
    phaseX = 0;
else 
    phaseX = pi;
end

if rand(1,1) > 0.5
    phaseY = 0;
else 
    phaseY = pi;
end

if rand(1,1) > 0.5
    phaseZ = 0;
else 
    phaseZ = pi;
end

%
for n = 1:num
    for m = 1:6   
        rangeLow = alphai.x(1,m) - alphai.x(1,m)/10;
        rangeHigh = alphai.x(1,m) + alphai.x(1,m)/10;
        alphaiTemp.x(1,m) = simPAF_gen_QRST_morf_variab(alphai.x(1,m), f, phaseX, rangeLow, rangeHigh, n);    
        
        rangeLow = alphai.y(1,m) - alphai.y(1,m)/10;
        rangeHigh = alphai.y(1,m) + alphai.y(1,m)/10;
        alphaiTemp.y(1,m) = simPAF_gen_QRST_morf_variab(alphai.y(1,m), f, phaseY, rangeLow, rangeHigh, n);
        
        rangeLow = alphai.z(1,m) - alphai.z(1,m)/10;
        rangeHigh = alphai.z(1,m) + alphai.z(1,m)/10;       
        alphaiTemp.z(1,m) = simPAF_gen_QRST_morf_variab(alphai.z(1,m), f, phaseZ, rangeLow, rangeHigh, n);
    end
    
       
[DIP, ~] = simPAF_DipoleGenerator(N,fs,F,alphaiTemp,bi,tetai,teta0); 
pqrstFrank = R0*Lambda*[DIP.x(51:610) ; DIP.y(51:610) ; DIP.z(51:610)]; 

% Dower transform
D = [  
-0.515  0.157   -0.917;
 0.044  0.164   -1.387;
 0.882  0.098   -1.277;
 1.213  0.127   -0.601;
 1.125  0.127   -0.086;
 0.831  0.076    0.230;
 0.632 -0.235    0.059;
 0.235  1.066   -0.132];
  
 pqrst_v1_v6_I_II = D*pqrstFrank; 

 pqrst(1:8,n,:) = pqrst_v1_v6_I_II(1:8,:); % V1, V2, V3, V4, V5, V6, I, II, 
 pqrst(9,n,:) = pqrst_v1_v6_I_II(8,:) - pqrst_v1_v6_I_II(7,:);     % III
 pqrst(10,n,:) = -(pqrst_v1_v6_I_II(7,:) + pqrst_v1_v6_I_II(8,:))/2; % aVR
 pqrst(11,n,:) =  pqrst_v1_v6_I_II(7,:) - pqrst_v1_v6_I_II(8,:)/2;   % aVL
 pqrst(12,n,:) =  pqrst_v1_v6_I_II(8,:) - pqrst_v1_v6_I_II(7,:)/2;   % aVF
    
 pqrst(13,n,:) = pqrstFrank(1,:); % X
 pqrst(14,n,:) = pqrstFrank(2,:); % Y
 pqrst(15,n,:) = pqrstFrank(3,:); % Z
end
 
 [~,qrsMaxInd] = max(pqrstFrank(1,1:300));
 
Qind = qrsMaxInd-50;
Rind = qrsMaxInd;
Sind = qrsMaxInd+50;
Tind = length(pqrstFrank(1,:));

pqrst_temp = pqrst;
pqrst(1:6,:,:) = pqrst(7:12,:,:);
pqrst(7:12,:,:) = pqrst_temp(1:6,:,:);
end

