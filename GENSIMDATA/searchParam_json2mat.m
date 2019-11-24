% convert search_params file from .json to .mat, more flexible
% Generate specific band width search parameter files

clear;
Dir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test10/searchParams/';
Filelist = dir([Dir,filesep,'*.json']);
n = length(Filelist);
NyqFreRan = [81.8345;1.0]; %upper; lower
NumBands = 5; % total number of bands
for i = 1:n
searchParams = jsondecode(fileread([Dir,filesep,Filelist(i).name]));
FreqRange = NyqFreRan;
%searchParams.band_num = i; % band number
save([Dir,filesep,'searchParams_MBLT',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
end