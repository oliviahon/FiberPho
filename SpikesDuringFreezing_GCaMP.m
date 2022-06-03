%to determine whether each spike occurs during freezing or mobility using
%spikes detected using SpikeThresholding_GCaMP
%use with Machine Results output from SimBA

%mice and correspnding trial must be in same order
mice=[1464 1454 1431 1432 1442 1444 1512 1513 1514 1521 1522 1531 1532 1533 1534 1535 1536 1537] ;
trials=[139 143 144 145 150 152 194 195 196 197 198 200 201 202 203 204 205 206];

exp='Learning'


for j=1:6
m=mice(j) 
trial=trials(j)

    %load in spike data   
    peakdatafile=sprintf('learning\\whole trace matlab output\\FearLearn_peakdataMAD_%d.xls',m)
    peakdata = readmatrix(peakdatafile,'Sheet','spike details');
    
    %load in freezing data from simba output
   freezedatafile= sprintf('VgatMachineResults\\Learning\\Trial%d.csv',trial)
   freezedata(:,1)=readmatrix(freezedatafile,'Range','A:A');%freezing data in col 1;
   freezedata(:,2)=readmatrix(freezedatafile,'Range','IN:IN');%freezing data in col 1;
   freezedata(1,:)=[];
   
    %converts frames to seconds by dividing frame # by 30 (recordings are
    %30fps)
    freezedata(:,3)= freezedata(:,1)/30;
    freezedata(:,4)=freezedata(:,1)/1800; %converts to minutes)
    
    %load in photometry traces 
    photomdatafile=sprintf('learning\\whole trace matlab output\\PagCeAvGATgcamp_%dwholetrace.xls',m);
    photomdata=readmatrix(photomdatafile);
    photomdata(1,:)=[];
    photomdata(:,1)=photomdata(:,1)-1; %fix time because these start at 2 min instead of 1, is an error, verified by looking at length of actual videos
    
   figname=sprintf('learning\\matlab figs\\spikeFreeze%d.fig',m)
    writefile=peakdatafile
    
    
%for each spike, find the frame closest in time in the simba output data,
%then write what the freezing status is during that frame (0 or 1)
for peak=1:length(peakdata)
    peaktime=find(freezedata(:,4)>=peakdata(peak,2));
    ind=peaktime(1);
    isfreeze(peak,1)=freezedata(ind,2);
end

peakdata=[peakdata isfreeze];

%make arrays of peaks during freezing and not during freezing, make second
%column for plotting tick marks
fyes=find(peakdata(:,6)==1);
fno=find(peakdata(:,6)==0);
     
freezeyes(:,1)=peakdata(fyes,2)
freezeyes(:,2)=5
freezeno(:,1)=peakdata(fno,2)
freezeno(:,2)=5

%downsample freezing data and convert so that it can be plotted as a
%scatter plot
freezedown=downsample(freezedata(:,:),3);
freezeplot=freezedown(:,2)*6;
freezeplot(freezeplot==0)=NaN;
freezedown=[freezedown freezeplot];

%plot everything all together for validation
figure
plot(photomdata(:,1), photomdata(:,4))
hold on
plot(freezedown(:,4), freezedown(:,5),'LineWidth',16)
hold on
plot(freezeyes(:,1), freezeyes(:,2),'|', 'Color', 'magenta','MarkerSize', 16)
hold on
plot(freezeno(:,1), freezeno(:,2),'|','Color', 'black', 'MarkerSize', 16)
savefig(figname)

%count number of spikes that occur during freezing and not
fyesn=length(fyes);
fnon=length(fno);

%calculate total time spent freezing by counting numer of freeze-pos frames
%and dividing by total frames to get % freezing and then mult by total
%seconds
timefreezep=(sum(freezedata(:,2)))/length(freezedata);
timemobilep=1-timefreezep;
totaltime=freezedata(length(freezedata),3);
timefreeze=timefreezep*totaltime;
timemobile=timemobilep*totaltime;

freqFreeze=fyesn/timefreeze;
freqMobile=fnon/timemobile;

spikequant(1,1)=timefreezep
spikequant(1,2)=timemobilep
spikequant(2,:)=totaltime
spikequant(3,1)=timefreeze
spikequant(3,2)=timemobile
spikequant(4,1)=freqFreeze
spikequant(4,2)=freqMobile

writematrix(isfreeze,writefile,'Sheet','is freeze')
writematrix(spikequant, writefile, 'Sheet', 'quant')

clearvars -except mice trials exp
close all
end

