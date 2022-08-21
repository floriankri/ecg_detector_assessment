function rr = simPAF_rrprocess(flo, fhi, flostd, fhistd, lfhfratio, hrmean, hrstd, sfrr, n)
%
% [] = simPAF_rrprocess() returns sinus rhythm RR intervals. 
%
% Original title: [] = rrprocess()
%
% RR intervals are generated according to the principle reported in the paper by 
% McSharry P, Clifford G, Tarassenko L & Smith L 2003 A dynamical model for generating 
% syntheticelectrocardiogram signals IEEE Transactions on Biomedical Engineering,
% 50(3), 289–294.
%
% Copyright (c) 2003 by Patrick McSharry & Gari Clifford, All Rights Reserved  
% Contact P. McSharry (patrick@mcsharry.net) or G. Clifford (gari@mit.edu)
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public License for more details. You should have received a copy of the 
% GNU General Public License along with this program; if not, write to the 
% Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
% MA 02111-1307  USA
%
% The function is also freely availble from Physionet at 
% https://physionet.org/content/ecgsyn/1.0.0/
%
w1 = 2*pi*flo;
w2 = 2*pi*fhi;
c1 = 2*pi*flostd;
c2 = 2*pi*fhistd;
sig2 = 1;
sig1 = lfhfratio;
rrmean = 60/hrmean;
rrstd = 60*hrstd/(hrmean*hrmean);

df = sfrr/n;
w = (0:n-1)'*2*pi*df;
dw1 = w-w1;
dw2 = w-w2;

Hw1 = sig1*exp(-0.5*(dw1/c1).^2)/sqrt(2*pi*c1^2);
Hw2 = sig2*exp(-0.5*(dw2/c2).^2)/sqrt(2*pi*c2^2);
Hw = Hw1 + Hw2;
Hw0 = [Hw(1:n/2); Hw(n/2:-1:1)];
Sw = (sfrr/2)*sqrt(Hw0);

ph0 = 2*pi*rand(n/2-1,1);
ph = [0; ph0; 0; -flipud(ph0) ]; 
SwC = Sw .* exp(1i*ph);
x = (1/n)*real(ifft(SwC));

xstd = std(x);
ratio = rrstd/xstd;
rr = rrmean + x*ratio;
end