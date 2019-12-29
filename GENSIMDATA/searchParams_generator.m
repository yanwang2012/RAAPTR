% load the entire xmaxmin from a json file then split the frequency into multiple bands
%%
clear;
filename = 'Nyquist.json';
outDir = '~/Research/PulsarTiming/SimDATA/MultiSource/Investigation/Test8/searchParams_Nyquist';
mkdir(outDir);
searchParams = jsondecode(fileread(filename));
NumBands = 5; %10; % total number of bands
bandwidth = searchParams.angular_velocity(1)/NumBands;% max frequency divided into 10 bands
FreqRange = searchParams.angular_velocity;
for i = 1:NumBands
    searchParams.angular_velocity(1) = bandwidth*i;
    if i ==1
        searchParams.angular_velocity(2) = 1;
    else
        searchParams.angular_velocity(2) = bandwidth*(i-1);
    end
    searchParams.band_num = i;
    save([outDir,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
end
