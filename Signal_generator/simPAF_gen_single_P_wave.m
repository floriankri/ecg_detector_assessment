function [P, fi0, fi1, fi2] = simPAF_gen_single_P_wave(t, k0,k1,k2,b0,b1,b2)
%
% [] = simPAF_gen_single_P_wave() returns a single P wave.
% A linear combination of Hermite functions is used to model P-waves in
% orthogonal leads as first proposed in Sörnmo et al. 1981.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%

% Hermite functions
    fi0 = k0*(1/(sqrt(b0*sqrt(pi))))*exp((-t.^2)/(2*b0.^2));
    fi1 = -1*k1*((sqrt(2))/(sqrt(b1*(sqrt(pi)))))*(t/b1).*exp((-t.^2)/(2*(b1.^2)));
    fi2 = k2*(1/(sqrt(2*b2*sqrt(pi))))*((2*((t.^2)/(b2.^2)))-1).*exp((-t.^2)/(2*(b2.^2)));
    % Resulting P wave
    P = fi0 + fi1 + fi2;
end