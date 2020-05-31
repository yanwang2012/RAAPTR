% convert search_params file from .json to .mat, more flexible
% Generate specific band width search parameter files

clear;
Dir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test11/searchParams/2bands/superNarrow';
Filelist = dir([Dir,filesep,'*.json']);
n = length(Filelist);
NyqFreRan = [81.8345;1.0]; %upper; lower
NumBands = 2; % total number of bands
for i = 1:n
searchParams = jsondecode(fileread([Dir,filesep,Filelist(i).name]));
FreqRange = NyqFreRan;
%searchParams.band_num = i; % band number
[~,filename,~]=fileparts(Filelist(i).name);
save([Dir,filesep,filename,'.mat'],'searchParams','NumBands','FreqRange');
end