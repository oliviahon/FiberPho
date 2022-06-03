% to calculate the average z-score of fiber pho signal during freezing
% bouts. appropriate for GCaMP terminal recordings or neurotransmitter
% biosensors

%using photometry code processed with FearLearn3ROIregression_GCaMP


mice=[1464 1454 1431 1432 1442 1444 1512 1513 1514 1521 1522 1531 1532 1533 1534 1535 1536 1537] ;
trials=[163 164 165 166 168 170 222 223 224 225 226 228 229 230 231 232 233 234];


exp='CuedFearCeA'

for j=1:18
m=mice(j) 
trial=trials(j)


    %load in freezing data from simba output
   freezedatafile= sprintf('VgatMachineResults\\Learning\\Trial%d.csv',trial)
   freezedata(:,1)=readmatrix(freezedatafile,'Range','A:A');%freezing data in col 1;
   freezedata(:,2)=readmatrix(freezedatafile,'Range','IN:IN');%freezing data in col 1;
   freezedata(1,:)=[];
   
    %converts frames to seconds by dividing frame # by 30 (recordings are
    %30fps)
    freezedata(:,3)= freezedata(:,1)/30;
    freezedata(:,4)=freezedata(:,1)/1800; %converts to minutes)
    
    % %downsample freezing data and convert so that it can be plotted as a
% %scatter plot
freezedown=downsample(freezedata(:,:),3);
freezeplot=freezedown(:,2)*5;
freezeplot(freezeplot==0)=NaN;
freezedown=[freezedown freezeplot];

    %load in photometry traces 
    photomdatafile=sprintf('learning\\whole trace matlab output\\VgatFearLearn_%dtrace.xls',m);
    photomdata=readmatrix(photomdatafile);
    photomdata(1,:)=[];
  
            %load in bout data
boutfile=sprintf('VgatMachineResults\\Learning\\FreezingBouts\\Trial%dbouts.xls',trial);
bouts=readmatrix(boutfile);

   figname=sprintf('learning\\matlab figs\\CeAFreeze%d.fig',m)
    writefile=boutfile
    
%we cut off the first minute of BL freezing data because we trim the first
%minute of photometry data
firstmin=find(bouts(:,2)<1)
tocut=length(firstmin)
bouts(1:tocut,:)=[]
[boutRows,boutCols] = size(bouts)
bouts(boutRows,:)=[]
boutRows=boutRows-1
    

%calculate avg during freezing bouts
for i=1:boutRows
    start=find(photomdata(:,5)>=bouts(i,1));
    stop=find(photomdata(:,5)>=bouts(i,2));
    bouts(i,5)=sum(photomdata(start(1):stop(1),2)); %sum during bout
     bouts(i,6)=sum(photomdata(start(1):stop(1),3));
     bouts(i,7)=length(photomdata(start(1):stop(1),2));%length of bout in data points
     
end
clear start stop


%create matrix of mobile epochs based on bouts
%triallength=freezedata(length(tracesL),4);

mobileend=boutRows+1
mobile1(1,1)=1
mobile1(2:mobileend,1)=bouts(:,2)
mobile2=bouts(:,1);
mobile2(length(mobile1),1)=mobile1(length(mobile1),1)+.01

mobile=[mobile1 mobile2];

for i=1:length(mobile)
    start=find(photomdata(:,5)>=mobile(i,1));
    stop=find(photomdata(:,5)>=mobile(i,2));
    mobile(i,3)=sum(photomdata(start(1):stop(1),2));
     mobile(i,4)=sum(photomdata(start(1):stop(1),3));
     mobile(i,5)=length(photomdata(start(1):stop(1),2));
end

avgs(1,1)=(sum(bouts(:,5)))/(sum(bouts(:,7)))%l freeze avg
avgs(1,2)=(sum(bouts(:,6)))/(sum(bouts(:,7)))%r freeze acg
avgs(2,1)=(sum(mobile(:,3)))/(sum(mobile(:,5)))%l mobile
avgs(2,2)=(sum(mobile(:,4)))/(sum(mobile(:,5)))%r mobile


figure
subplot(2,1,1)
plot(photomdata(2:end,5), photomdata(2:end,2));
hold on
plot(freezedown(:,4), freezedown(:,5),'LineWidth',16)
xlabel('time min')
ylabel('zscore')
title('Left')

subplot(2,1,2)
plot(photomdata(2:end,5), photomdata(2:end,3));
hold on
plot(freezedown(:,4), freezedown(:,5),'LineWidth',16)
xlabel('time min')
ylabel('zscore')
title('Right')
savefig(figname)

writematrix(bouts, boutfile,'Sheet','bouts w freez')
writematrix(mobile, boutfile,'Sheet','mobile')
writematrix(freezedown, boutfile,'Sheet','freeze 4 plot')
writematrix(avgs, boutfile,'Sheet','avgs weighted')

clearvars -except mice trials exp
close all
end


