function multileadNoise = simPAF_gen_noise(ecgLength, noiseType, noiseRMS)
%
% multileadNoise = simPAF_gen_noise() returns physiological noises obtained
% from the MIT-BIH Noise Stress Test Database. The signals of the
% the MIT–BIH Noise Stress Test database is used as a noise source for 
% orthogonal leads X and Y, whereas noise for lead Z is obtained as the sum 
% of squares of these two noises of leads X and Y. Then, an inverse of 
% Dowers’s transformation matrix is applied to generate noise in the remaining 
% 12 ECG leads. Ultimately, the noise is rescaled to the desired RMS value.
%
% Copyright (C) 2017  Andrius Petrenas
% Biomedical Engineering Institute, Kaunas University of Technology
%
% Available under the GNU General Public License version 3
% - please see the accompanying file named "LICENSE"
%
% noiseType: value from 0 to 4
% 0 - no noise added (noise RMS = 0 mV)
% 1 - motion artefacts
% 2 - electrode movement
% 3 - baseline wander
% 4 - mixture of noises
%
noiseLength = 1805556; % Noise length
if ecgLength < noiseLength
    noiseStart = randi([1 (noiseLength-ecgLength)]); % take starting point randomly
    noiseStart = 1;
else
    noiseStart = randi([1 noiseLength/2]); %1805556 - cut half of length noise segments
    noiseStart = 1;
    cycles = ceil(ecgLength/(noiseLength/2));
    noiseTemp = [];
end

switch noiseType 
    case 0     % no noise added 
       multileadNoise = zeros(15, ecgLength);  
    case 1     % motion artefacts
       data = load('DATA_noises_real'); 
       if ecgLength < noiseLength
          multileadNoise = data.motion_artefacts(:, noiseStart+1:noiseStart+ecgLength);
       else
          multileadNoise = data.motion_artefacts(:, noiseStart+1:noiseStart+(noiseLength/2));
          for i = 1:cycles
              noiseTemp = [noiseTemp multileadNoise]; 
          end
          multileadNoise = noiseTemp(:,1:ecgLength);
       end
       
    case 2     % electrode movement
       data = load('DATA_noises_real'); 
       if ecgLength < noiseLength
          multileadNoise = data.electrode_movement(:, noiseStart+1:noiseStart+ecgLength);
       else
          multileadNoise = data.electrode_movement(:, noiseStart+1:noiseStart+(noiseLength/2));
          for i = 1:cycles
              noiseTemp = [noiseTemp multileadNoise]; 
          end
          multileadNoise = noiseTemp(:,1:ecgLength);
       end

    case 3     % baseline wander
       data = load('DATA_noises_real'); 
       if ecgLength < noiseLength
          multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+ecgLength);
       else
          multileadNoise = data.baseline_wander(:, noiseStart+1:noiseStart+(noiseLength/2));
          for i = 1:cycles
              noiseTemp = [noiseTemp multileadNoise]; 
          end
          multileadNoise = noiseTemp(:,1:ecgLength);
       end

    case 4     % mixture of noises
       data = load('DATA_noises_real'); 
       if ecgLength < noiseLength
          multileadNoise = data.mixture_of_noises(:, noiseStart+1:noiseStart+ecgLength);
       else
          multileadNoise = data.mixture_of_noises(:, noiseStart+1:noiseStart+(noiseLength/2));
          cycles
          for i = 1:cycles
              noiseTemp = [noiseTemp multileadNoise]; 
          end
          multileadNoise = noiseTemp(:,1:ecgLength);
       end   
end

if noiseType > 0
    % Adjust to desired noise RMS value
    for i = 1:15
        multileadNoise(i,:) = noiseRMS*(multileadNoise(i,:)/std(multileadNoise(i,:)));   
    end
end
end