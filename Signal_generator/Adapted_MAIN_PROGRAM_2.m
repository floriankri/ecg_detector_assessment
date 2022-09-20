close all
clear 
clc

rrLength_array = [40];
APBrate_array = [0];
onlyRR_array = [0]; 

medEpis_array = [15];
stayInAF_array = [1-log(2)/15];
AFburden_array = [0 1];

noiseType_array = [0 1 2 3];
noiseRMS_array = [0 0.5];

realRRon_array = [1];
realVAon_array = [1];
realAAon_array = [1];
repeats = 2;

newdata = Adapted_simPAF_ECG_generator_iterator(rrLength_array,APBrate_array,onlyRR_array,medEpis_array,stayInAF_array,AFburden_array,noiseType_array,noiseRMS_array,realRRon_array,realVAon_array,realAAon_array, repeats);