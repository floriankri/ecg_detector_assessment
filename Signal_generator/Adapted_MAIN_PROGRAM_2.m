close all
clear 
clc

rrLength_array = [50];
APBrate_array = [0.10];
onlyRR_array = [0]; 

medEpis_array = [15];
stayInAF_array = [1-log(2)/15];
AFburden_array = [0.8];

noiseType_array = [3 4];
noiseRMS_array = [0.02 0.03];

realRRon_array = [1];
realVAon_array = [0];
realAAon_array = [0];

newdata = Adapted_simPAF_ECG_generator_iterator(rrLength_array,APBrate_array,onlyRR_array,medEpis_array,stayInAF_array,AFburden_array,noiseType_array,noiseRMS_array,realRRon_array,realVAon_array,realAAon_array);