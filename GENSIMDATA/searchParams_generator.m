% load the entire xmaxmin width from a txt file then split the frequency into 10 bands
%%
% load('searchParams_GWBsimDataSKA.mat')
% writematrix(xmaxmin,'xmaxmin.txt','delimiter','tab');
xmaxmin = readmatrix('xmaxmin.txt');
bandwidth = xmaxmin(3,1)/10;
for i = 1:10
    xmaxmin(3,1) = bandwidth*i;
    if i ==1
        xmaxmin(3,2) = 1;
    else
        xmaxmin(3,2) = bandwidth*(i-1);
    end
    nband = i;
    save(['searchParams_GWBsimDataSKA',num2str(i),'.mat'],'xmaxmin','nband');
end
