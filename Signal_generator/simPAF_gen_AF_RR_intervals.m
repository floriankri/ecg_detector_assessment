function rr = simPAF_gen_AF_RR_intervals(lamba,nRR)
% rr = simPAF_gen_AF_RR_intervals() returns AF RR series modelled according
% to the paper by Corino et al. An atrioventricular node model for analysis
% of the ventricular response during atrial fibrillation. IEEE Transactions
% on Biomedical Engineering. 2011, 58(12), 3386-3395.
%
% Copyright (C) 2011 Valentina D. A. Corino*, Frida Sandberg**, 
% Luca T. Mainardi*, Leif Sornmo**
% *Department of Bioengineering, Politecnico di Milano;
% **Department of Biomedical Engineering and 
% Center of Integrative Electrocardiology, Lund University;
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
if nRR < 70       
    time = 1;
else
    time = nRR/70;
end

prob_alpha = 0.6;
P = [0.15 0.25];
difftau = 0.2;
slope_beta = 10;

tau1_in=polyval(P,1);
tau1=polyval(P,0:.001:3);
P2=[P(1) P(2)+difftau];
tau2=polyval(P2,0:.001:3);
tau2_in=polyval(P2,1);

Vt=-40;
Vr=-90;
dVdt=0;

aa=exprnd(1/lamba,1000000,1);
t_a=cumsum(aa);
aa=aa(t_a<time*60);
atr_t=0;
avj_tmA=0;
vtr_tmA=0;
phase='phase4';
Vm=Vr;
Ts=0.001; 
t=0;
nextAA=aa(1);
i_R=0;
vtr_t=0;
vtr_nA=0;
avj_t=0;avj_nA=0;
R=zeros(length(aa),1);

nAA=0;
prob=zeros(length(aa),1);
prob1=round(length(aa)*prob_alpha); 
prob(1:prob1)=1;
prob=prob(randperm(length(prob)));
prob_vera=zeros(size(prob));
tautau=zeros(size(aa));

tempo_x=0:.001:3;

beta=-slope_beta*tempo_x+1;
beta(beta<0)=0;

prob2=round(length(aa)*beta); 
Nbeta=length(find(prob2>0));
prob_beta=zeros(length(aa),length(Nbeta));

for ii=1:Nbeta 
    prob_beta(1:prob2(ii),ii)=1;
    prob_beta(:,ii)=prob_beta(randperm(length(aa)),ii);
end
 
blocked_beta=zeros(length(aa),1);
blocked_beta_t=zeros(length(tempo_x),1);
non_blocked_beta_t=zeros(length(tempo_x),1);
tempo_tot=zeros(length(aa),1);

while nAA<length(aa)
    
%%% function UpdateAtTs
% At every msec this is what happens: Vm is linearly increased
% by dVdt*Ts if in phase4; or the avj_t is increased

    t=t+Ts; % update time
    atr_t=atr_t+Ts;
    vtr_t=vtr_t+Ts;
    
    avj_tmA=avj_tmA+Ts;
    vtr_tmA=vtr_tmA+Ts;
    
    % update timers
    if strcmp(phase,'phase4') % otherwise it's in refractory period
        Vm=Vm+dVdt*Ts;         
    else
        avj_t=avj_t+Ts;
    end
    
%%% function AnteHitAVJ
% an AA arrives at the AVJ: Vm is increased of deltaV if in phase4 or 
% tau (=RP) is prolonged

    if atr_t>=nextAA
        
%         t1=atr_t;
        atr_t=0;
        nAA=nAA+1; % increments AF counter
        if nAA<length(aa)
            nextAA=aa(nAA+1); % next AA
        end

        if strcmp(phase,'phase4') % otherwise it's in refractory period 
            
            beats=find(R>0);
            if beats>0
                tempo=t-R(beats(end))-tau;
                [m, pos]=min(abs(tempo_x-tempo));
                if pos<=Nbeta
                    if prob_beta(nAA,pos)==0
                        deltaV=(Vt-Vr)+1; 
                        Vm=Vm+deltaV;
                        non_blocked_beta_t(pos)=non_blocked_beta_t(pos)+1;
                    else blocked_beta(nAA)=1;
                        blocked_beta_t(pos)=blocked_beta_t(pos)+1;
                        tempo_tot(nAA)=tempo;
                    end 
                else 
                    deltaV=(Vt-Vr)+1; 
                    Vm=Vm+deltaV;
                    non_blocked_beta_t(pos)=non_blocked_beta_t(pos)+1;
                end
            else 
                deltaV=(Vt-Vr)+1; 
                Vm=Vm+deltaV;
            end            
        end       
    end  
%%% function VtrSense
% there is a wave in the ventricle (vtr_nA>0) and the ventricle 
% is not in refractory period (vtr_tmA>=AntDly)
 
    if vtr_nA>0 
        i_R=i_R+1;
        R(i_R)=t;
        vtr_tmA=0;
        vtr_t=0;
        vtr_nA=0;
    end
      
%%% function ActivateAVJ + StartAVJref if in phase4:
% when Vm>=Vt AVJ is activated (avj_nA=avj_nA+1;) and then the
% RP starts in AVJ

% OR StartAVJph4 if in phase0
% when avj_t>tau, i.e. the RP of AVJ finishes, phase4 again

    if strcmp(phase,'phase4') % otherwise it's in refractory period
        if Vm>=Vt       
           %%% ActivateAVJ 
            avj_nA=avj_nA+1; 
            vtr_nA=vtr_nA+1;     
            phase='phase0';
            avj_t=0;
            
            if prob(nAA)==0 
                
                beats=find(R>0);
                if beats>0
                    new_RR=t-R(beats(end));
                    [m xx]=min(abs(tempo_x-new_RR));
                    tau=tau2(xx);
                    prob_vera(nAA)=2;
                else
                    tau=tau2_in;
                    prob_vera(nAA)=2;
                end               
            else 
                beats=find(R>0);
                if beats>0
                    new_RR=t-R(beats(end));
                    [m, xx]=min(abs(tempo_x-new_RR));
                    tau=tau1(xx);
                    prob_vera(nAA)=1;
                else
                    tau=tau1_in;
                    prob_vera(nAA)=1;
                end
            end
            tautau(nAA)=tau;
        end 
        
    elseif exist('tau1','var')==1                  
        if avj_t>tau
            %%% StartAVJph4
            phase='phase4';
            Vm=Vr;
        end

    end
   
end
R=R(R>0);
rr=diff(R);
