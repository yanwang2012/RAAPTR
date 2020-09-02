% Read a searchParams.json file and fluctuate the band boundary,
% then generate new search parameter set.

% QYQ 08/2020

%%
clear

input_file_name = '/Users/qyq/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/searchParams1.json'; % file defines the band edge to be shifted
s1 = jsondecode(fileread(input_file_name));
NumBands = 2; % number of bands
% create searchParams struct
searchParams.alpha = s1.alpha;
searchParams.delta = s1.delta;
searchParams.angular_velocity = s1.angular_velocity;
searchParams.phi0 = s1.phi0;
searchParams.inclination = s1.inclination;
searchParams.amplitude = s1.amplitude;
searchParams.polarization = s1.polarization;
searchParams.band_num = s1.band_num;


full_range_file_name = 'Nyquist.json'; % file defines the whole search range
fullRange = jsondecode(fileread(full_range_file_name));
FreqRange = fullRange.angular_velocity;

edge = s1.angular_velocity(1,1);
fluctuation = 10 * randn; % fluctuation range
new_edge = edge + fluctuation;

% update searchParams for 1st band with new band edge
s1.angular_velocity(1,1) = new_edge;
searchParams.angular_velocity = s1.angular_velocity;

[filePath,newFile,~] = fileparts(input_file_name);
newFile = regexp(newFile,'[a-zA-Z]*','match');
save([filePath,filesep,char(newFile),'Rand1','.mat'],'searchParams','NumBands','FreqRange');

fullRange.angular_velocity(2,1) = new_edge;
fullRange.band_num = 2;

% update searchParams for 2nd band
searchParams.angular_velocity = fullRange.angular_velocity;
searchParams.band_num = fullRange.band_num;

save([filePath,filesep,char(newFile),'Rand2','.mat'],'searchParams','NumBands','FreqRange')