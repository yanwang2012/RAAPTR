% load the entire xmaxmin width from a json file then split the frequency into 10 bands
%%
clear;
filename = 'Nyquist.json';
outDir = 'searchParams_Nyquist';
mkdir(outDir);
searchParams = jsondecode(fileread(filename));
NumBands = 10; % total number of bands
bandwidth = searchParams.angular_velocity(1)/NumBands;% max frequency divided into 10 bands
FreqRange = searchParams.angular_velocity;
for i = 1:10
    searchParams.angular_velocity(1) = bandwidth*i;
    if i ==1
        searchParams.angular_velocity(2) = 1;
    else
        searchParams.angular_velocity(2) = bandwidth*(i-1);
    end
    searchParams.band_num = i;
    save([outDir,filesep,'searchParams_Nyquist',num2str(i),'.mat'],'searchParams','NumBands','FreqRange');
end
