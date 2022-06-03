%written for cued fear learning with TTL timestamps
%data recorded with Neurophotometrics FP3001, Trig 1
%packages needed= curve fitting toolbox, filter builder

%cohort info
mice=[1513 1521 1532 1533 1534 1536]  
shockmice=[1521 1532 1534]
ctrlmice=[1513 1533 1536]
ledstart=[1513 470;1521 415;1532 470;1533 470;1534 470;1536 470] 

CEALsh=[]
CEARsh=[]
PAGsh=[]
CEALc=[]
CEARc=[]
PAGc=[]


for j=1:length(mice)
    %tell it the file names to load
    m=mice(j) 
       
    datafile=sprintf('PagCeAVgatLearn_%d_1.csv',m)
    timefile=sprintf('PagCeAVgatLearn_%d_ts0.csv',m)
    FP.rawdata = readmatrix(datafile);
    FP.timestamps =readmatrix(timefile);
    FP.timestamps(:,2)=[] %delete second column of ts
  
    FP.background1=3075
    FP.background2=3080
    FP.background3=3085
    
%filter specs
Fstop=5;
Fsample=20;
order=10;

%for naming pictures
exp='VgatFear2cLearn'

%name all figures that we will save below
picz=sprintf('%s_%dzscore', exp,m)
piccorrection =sprintf( '%s_%dcorrection',exp,m)
picchopped=sprintf('%s_%dchop',exp,m)
motioncorrected= sprintf('%s_%dmotioncorrect',exp,m)
regression= sprintf('%s_%dregression',exp,m)

%cuts off last data point if total length is an odd number
trim=find(FP.rawdata(:,1)>=(FP.timestamps(end,1)+120000));
if rem(trim(1,1),2) == 0
    stop=trim
else
stop=(trim(1,1)-1);
end

% plots raw data
f1 = figure;
for i = 1:4
    subplot(4,1,i)
    plot(FP.rawdata(:,i))
end


%deinterleave data and adjust time to start based on first timestamp
if ledstart(j,2)==415

% % starts wtih 415
FP.calcium_dependent(:,1) = downsample(FP.rawdata(2:stop,1),2); % time
FP.calcium_dependent(:,2) = downsample(FP.rawdata(2:stop,2),2); % CeA L
FP.calcium_dependent(:,3) = downsample(FP.rawdata(2:stop,3),2); % CeA R
FP.calcium_dependent(:,4) = downsample(FP.rawdata(2:stop,4),2); % PAG
FP.calcium_dependent(:,5) = (FP.calcium_dependent(:,1)-(FP.timestamps(1,1)-120000))./60000 %convert time to min based on first TS

FP.isos(:,1) = downsample(FP.rawdata(1:stop,1),2); % time
FP.isos(:,2) = downsample(FP.rawdata(1:stop,2),2); % CeA L
FP.isos(:,3) = downsample(FP.rawdata(1:stop,3),2); % CeA R
FP.isos(:,4) = downsample(FP.rawdata(1:stop,4),2); % PAG
FP.isos(:,5) = (FP.calcium_dependent(:,1)-(FP.timestamps(1,1)-120000))./60000

else
%starts wtih 470
FP.calcium_dependent(:,1) = downsample(FP.rawdata(1:stop,1),2); 
FP.calcium_dependent(:,2) = downsample(FP.rawdata(1:stop,2),2); % l
FP.calcium_dependent(:,3) = downsample(FP.rawdata(1:stop,3),2); % r
FP.calcium_dependent(:,4) = downsample(FP.rawdata(1:stop,4),2); % pag
FP.calcium_dependent(:,5) = (FP.calcium_dependent(:,1)-(FP.timestamps(1,1)-120000))./60000

FP.isos(:,1) = downsample(FP.rawdata(2:stop,1),2); 
FP.isos(:,2) = downsample(FP.rawdata(2:stop,2),2); 
FP.isos(:,3) = downsample(FP.rawdata(2:stop,3),2);
FP.isos(:,4) = downsample(FP.rawdata(2:stop,4),2); 
FP.isos(:,5) = (FP.calcium_dependent(:,1)-(FP.timestamps(1,1)-120000))./60000
end


%subtract the background for each column from the raw data
FP.backgroundsubtracted(:,1)= FP.calcium_dependent(:,1); 
FP.backgroundsubtracted(:,2)= (FP.calcium_dependent(:,2)-FP.background1);
FP.backgroundsubtracted(:,3)= (FP.calcium_dependent(:,3)-FP.background2);
FP.backgroundsubtracted(:,4)= (FP.calcium_dependent(:,4)-FP.background3);
FP.backgroundsubtracted(:,5)= FP.calcium_dependent(:,5);

FP.backgroundsubtractedisos(:,1)= FP.isos(:,1); 
FP.backgroundsubtractedisos(:,2)= (FP.isos(:,2)-FP.background1);
FP.backgroundsubtractedisos(:,3)= (FP.isos(:,3)-FP.background2);
FP.backgroundsubtractedisos(:,4)= (FP.isos(:,4)-FP.background3);
FP.backgroundsubtractedisos(:,5)= FP.isos(:,5);

%design filter based on parameters specified above
butterfive = designfilt('lowpassiir', 'FilterOrder', order, 'HalfPowerFrequency', Fstop, 'SampleRate', Fsample, 'DesignMethod', 'butter');
filt=butterfive;

FP.filtered(:,1)= FP.backgroundsubtracted(:,1);
FP.filtered(:,2)= filtfilt(filt,FP.backgroundsubtracted(:,2));
FP.filtered(:,3)= filtfilt(filt,FP.backgroundsubtracted(:,3));
FP.filtered(:,4)= filtfilt(filt,FP.backgroundsubtracted(:,4));
FP.filtered(:,5)= FP.backgroundsubtracted(:,5);

FP.filteredisos(:,1)= FP.backgroundsubtractedisos(:,1);
FP.filteredisos(:,2)= filtfilt(filt,FP.backgroundsubtractedisos(:,2));
FP.filteredisos(:,3)= filtfilt(filt,FP.backgroundsubtractedisos(:,3));
FP.filteredisos(:,4)= filtfilt(filt,FP.backgroundsubtractedisos(:,4));
FP.filteredisos(:,5)= FP.backgroundsubtracted(:,5);

%cut off first minute of BL then curve fit
start=find(FP.calcium_dependent(:,5)>=1);

f4 = figure; %fit signal to bioexp and subtract and normalize 470 data
subplot(4,3,1);
temp_fit = fit(FP.filtered(start:end,5),FP.filtered(start:end,2),'exp2'); 
plot(temp_fit,FP.filtered(start:end,5),FP.filtered(start:end,2))
title('Uncorrected CeA L 470','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,2);
temp_fit2 = fit(FP.filtered(start:end,5),FP.filtered(start:end,3),'exp2'); 
plot(temp_fit2,FP.filtered(start:end,5),FP.filtered(start:end,3))
title('Uncorrected CeA R 470','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,3);
temp_fit3 = fit(FP.filtered(start:end,5),FP.filtered(start:end,4),'exp2'); 
plot(temp_fit3,FP.filtered(start:end,5),FP.filtered(start:end,4))
title('Uncorrected PAG 470','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,4);
temp_fit4 = fit(FP.filteredisos(start:end,5),FP.filteredisos(start:end,2),'exp2'); 
plot(temp_fit4,FP.filteredisos(start:end,5),FP.filteredisos(start:end,2))
title('Uncorrected CeA L isos','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,5);
temp_fit5 = fit(FP.filteredisos(start:end,5),FP.filteredisos(start:end,3),'exp2'); 
plot(temp_fit5,FP.filteredisos(start:end,5),FP.filteredisos(start:end,3))
title('Uncorrected CeA R isos','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,6);
temp_fit6 = fit(FP.filteredisos(start:end,5),FP.filteredisos(start:end,4),'exp2'); 
plot(temp_fit6,FP.filteredisos(start:end,5),FP.filteredisos(start:end,4))
title('Uncorrected PAG isos','fontsize',16)
ylabel('F (au)','fontsize',14)

subplot(4,3,7);
plot(FP.filtered(start:end,5),100*(FP.filtered(start:end,2)-temp_fit(FP.filtered(start:end,5)))./temp_fit(FP.filtered(start:end,5)))
title('Corrected CeA L 470','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)

subplot(4,3,8);
plot(FP.filtered(start:end,5),100*(FP.filtered(start:end,3)-temp_fit2(FP.filtered(start:end,5)))./temp_fit2(FP.filtered(start:end,5)))
title('Corrected CeA R 470','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)

subplot(4,3,9);
plot(FP.filtered(start:end,5),100*(FP.filtered(start:end,4)-temp_fit3(FP.filtered(start:end,5)))./temp_fit3(FP.filtered(start:end,5)))
title('Corrected PAG 470','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)

subplot(4,3,10);
plot(FP.filteredisos(start:end,5),100*(FP.filteredisos(start:end,2)-temp_fit4(FP.filteredisos(start:end,5)))./temp_fit4(FP.filteredisos(start:end,5)))
title('Corrected CeA L 415','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)

subplot(4,3,11);
plot(FP.filteredisos(start:end,5),100*(FP.filteredisos(start:end,3)-temp_fit5(FP.filteredisos(start:end,5)))./temp_fit5(FP.filteredisos(start:end,5)))
title('Corrected CeA R 415','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)

subplot(4,3,12);
plot(FP.filteredisos(start:end,5),100*(FP.filteredisos(start:end,4)-temp_fit6(FP.filteredisos(start:end,5)))./temp_fit6(FP.filteredisos(start:end,5)))
title('Corrected PAG 415','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('dF/F %','fontsize',14)
savefig(piccorrection)

%gives background subtracted data- the fit of background subtracted
%data/fit. this is df/f
FP.fitsubtract(:,1) = (FP.filtered(start:end,1));
FP.fitsubtract(:,2) = (FP.filtered(start:end,2)-temp_fit(FP.filtered(start:end,5)))./temp_fit(FP.filtered(start:end,5));
FP.fitsubtract(:,3) =(FP.filtered(start:end,3)-temp_fit2(FP.filtered(start:end,5)))./temp_fit2(FP.filtered(start:end,5));
FP.fitsubtract(:,4) = (FP.filtered(start:end,4)-temp_fit3(FP.filtered(start:end,5)))./temp_fit3(FP.filtered(start:end,5));
FP.fitsubtract(:,5) = FP.calcium_dependent(start:end,5);

FP.fitsubtractisos(:,1) = (FP.filteredisos(start:end,1));
FP.fitsubtractisos(:,2) = (FP.filteredisos(start:end,2)-temp_fit4(FP.filteredisos(start:end,5)))./temp_fit4(FP.filteredisos(start:end,5));
FP.fitsubtractisos(:,3) =(FP.filteredisos(start:end,3)-temp_fit5(FP.filteredisos(start:end,5)))./temp_fit5(FP.filteredisos(start:end,5));
FP.fitsubtractisos(:,4) = (FP.filteredisos(start:end,4)-temp_fit6(FP.filteredisos(start:end,5)))./temp_fit6(FP.filteredisos(start:end,5));
FP.fitsubtractisos(:,5) = FP.calcium_dependent(start:end,5);

%zscore
FP.zscore(:,1)= FP.fitsubtract(:,1);
FP.zscore(:,2)=zscore(FP.fitsubtract(:,2));
FP.zscore(:,3)=zscore(FP.fitsubtract(:,3));
FP.zscore(:,4)=zscore(FP.fitsubtract(:,4));
FP.zscore(:,5)=FP.fitsubtract(:,5);

FP.zscorei(:,1)= FP.fitsubtractisos(:,1);
FP.zscorei(:,2)=zscore(FP.fitsubtractisos(:,2));
FP.zscorei(:,3)=zscore(FP.fitsubtractisos(:,3));
FP.zscorei(:,4)=zscore(FP.fitsubtractisos(:,4));
FP.zscorei(:,5)=FP.fitsubtractisos(:,5);

f7=figure %plots z scored data against time 
subplot(3,1,1)
plot(FP.zscore(:,5),FP.zscore(:,2),'b')
hold on
plot(FP.zscorei(:,5),FP.zscorei(:,2),'m')
title('Z-scored CeA L','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('Z-score','fontsize',14)

subplot(3,1,2)
plot(FP.zscore(:,5),FP.zscore(:,3),'b')
hold on
plot(FP.zscorei(:,5),FP.zscorei(:,3),'m')
title('Z-scored CeA R','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('Z-score','fontsize',14)

subplot(3,1,3)
plot(FP.zscore(:,5),FP.zscore(:,4),'b')
hold on
plot(FP.zscorei(:,5),FP.zscorei(:,4),'m')
title('Z-scored PAG','fontsize',16)
xlabel('Time (m)','fontsize',14)
ylabel('Z-score','fontsize',14)
savefig(picz)

% using non negative robust linear regression
fitdata1 = fit(FP.zscorei(:,2),FP.zscore(:,2),fittype('poly1'),'Robust','on');
fitdata2 = fit(FP.zscorei(:,3),FP.zscore(:,3),fittype('poly1'),'Robust','on');
fitdata3 = fit(FP.zscorei(:,4),FP.zscore(:,4),fittype('poly1'),'Robust','on');

% Plot fit
figure
subplot(3,1,1)
hold on
plot(FP.zscorei(:,2),FP.zscore(:,2),'k.')
plot(fitdata1,'b')
title('Lin reg CEA L')
hold off

subplot(3,1,2)
hold on
plot(FP.zscorei(:,3),FP.zscore(:,3),'k.')
plot(fitdata2,'b')
title('lin reg CEA R')
hold off

subplot(3,1,3)
hold on
plot(FP.zscorei(:,4),FP.zscore(:,4),'k.')
plot(fitdata3,'b')
title('lin reg PAG')
hold off
savefig(regression)

%fit isos to ca signal
isosfitted(:,1)=fitdata1(FP.zscorei(:,2));
isosfitted(:,2)=fitdata2(FP.zscorei(:,3));
isosfitted(:,3)=fitdata3(FP.zscorei(:,4));

%subtract fitted signal
FP.corrected(:,1)=FP.zscore(:,1);
FP.corrected(:,2)=(FP.zscore(:,2)-isosfitted(:,1));
FP.corrected(:,3)=(FP.zscore(:,3)-isosfitted(:,2));
FP.corrected(:,4)=(FP.zscore(:,4)-isosfitted(:,3));
FP.corrected(:,5)=FP.zscore(:,5);

figure
subplot(2,3,1)
plot(isosfitted(:,1),'m')
hold on
plot(FP.zscore(:,2),'b')
title('fitted isos and signal CeA L')
legend('fitted isos', '470')
hold off

subplot(2,3,2)
plot(isosfitted(:,2),'m')
hold on
plot(FP.zscore(:,3),'b')
title('fitted isos and signal CeA R')
legend('fitted isos', '470')
hold off

subplot(2,3,3)
plot(isosfitted(:,3),'m')
hold on
plot(FP.zscore(:,4),'b')
title('fitted isos and signal PAG')
legend('fitted isos', '470')
hold off

subplot(2,3,4)
plot(FP.corrected(:,5),FP.corrected(:,2),'k')
title('Motion Corrected signal CeA L')
xlabel('Time')
ylabel('z-score')

subplot(2,3,5)
plot(FP.corrected(:,5),FP.corrected(:,3),'k')
title('Motion Corrected signal CeA R')
xlabel('Time')
ylabel('z-score')

subplot(2,3,6)
plot(FP.corrected(:,5),FP.corrected(:,4),'k')
title('Motion Corrected signal PAG')
xlabel('Time')
ylabel('z-score')
savefig(motioncorrected)

%create peri-event plots based on arduino timestamps

%2nd will be be tone onset, 3rd is shock. 4th tone onset and so on.
FP.tonetime(:,1)= downsample(FP.timestamps(:,1),2);
FP.bl = 400; %number of datapoints to look at before timestamp, 400= 20s
FP.seg_dur = 1800;
%number of datapoints to look at after timestamp. 1800= 90s

    FP.ERF = [];
for shock = 1:length(FP.tonetime); 
        temp = find(FP.corrected(:,1) >= FP.tonetime(shock)); %find all FP data indices that are greater than or equal to timestamp
        temp = temp(1); %take the first one
        
    for fiber = 1:3%does this for each fiber
        FP.ERF(fiber).data(:,shock) = FP.corrected(temp-FP.bl:temp+FP.seg_dur,fiber+1); %save in nice neat datastructure
      
        %and then normalize each trace so they all start at 0
      FP.ERF(fiber).data(:,shock) = FP.ERF(fiber).data(:,shock) - mean(FP.ERF(fiber).data(1:FP.bl,shock));
    end
    clear puff fiber temp
end
%and then we calculate time for these chopped traces
FP.fps = 1000/mean(diff(FP.corrected(:,1))); %note, this way of calculating FPS fails if you started and stopped your recording (from the driver box)
FP.ERF_time= linspace(-FP.bl/FP.fps,FP.seg_dur/FP.fps,length(FP.ERF(1).data));
    
f7=figure
subplot(3,1,1)
plot(FP.ERF_time,FP.ERF(1).data) 
 xlabel('Time Relative to Tone Onset (s)','fontsize',14)
 ylabel('Z-Score','fontsize',14)
 title('CeA L','fontsize',16)
 xline(0) %plots a vertical line at tone onset
 xline(28) %plots vertical line at shock start
hold on
    plot(FP.ERF_time,mean(FP.ERF(1).data,2),'k','LineWidth',2)

subplot(3,1,2)
plot(FP.ERF_time,FP.ERF(2).data) 
hold on
 xlabel('Time Relative to Tone Onset(s)','fontsize',14)
 ylabel('Z-Score','fontsize',14)
 title('CeA R','fontsize',16)   
  xline(0)
 xline(28)
hold on
    plot(FP.ERF_time,mean(FP.ERF(2).data,2),'k','LineWidth',2)
    
subplot(3,1,3)
plot(FP.ERF_time,FP.ERF(3).data) 
hold on
 xlabel('Time Relative to ToneShock (s)','fontsize',14)
 ylabel('Z-Score','fontsize',14)
 title('PAG','fontsize',16)   
  xline(0)
 xline(28)
hold on
    plot(FP.ERF_time,mean(FP.ERF(3).data,2),'k','LineWidth',2)
savefig(picchopped)


% compile data across cohort
if ismember(m,shockmice)
FP.ERF(1).data(1,1:end)=m;
FP.ERF(2).data(1,1:end)=m;
FP.ERF(3).data(1,1:end)=m;
   CEALsh=[CEALsh,FP.ERF(1).data];
   CEARsh=[CEARsh,FP.ERF(2).data];
   PAGsh=[PAGsh,FP.ERF(3).data];
   
elseif ismember(m,ctrlmice)
FP.ERF(1).data(1,1:end)=m;
FP.ERF(2).data(1,1:end)=m;
FP.ERF(3).data(1,1:end)=m;
   CEALc=[CEALc,FP.ERF(1).data];
   CEARc=[CEARc,FP.ERF(2).data];
   PAGc=[PAGc,FP.ERF(3).data];
end

time=FP.ERF_time;

trace=downsample((FP.corrected),4);
tracefilename=sprintf('%s_%dtrace.xls',exp,m)
writematrix(trace,tracefilename)
close all
clearvars -except mice ledstart time shockmice ctrlmice CEALsh CEARsh PAGsh CEALc CEARc PAGc 

end

%downsample compiled peri-event data for export
dCEALsh=downsample((CEALsh),2);
dCEARsh=downsample((CEARsh),2);
dPAGsh=downsample((PAGsh),2);
dCEALc=downsample((CEALc),2);
dCEARc=downsample((CEARc),2);
dPAGc=downsample((PAGc),2);
dtime=(downsample((time),2))';

%write processed data out in excel
filename=('Vgat2cAntagLearn.xls')
writematrix(dCEALsh,filename,'Sheet','Cea L sh');
writematrix(dCEARsh,filename,'Sheet','Cea R sh');
writematrix(dPAGsh,filename,'Sheet','PAG sh');
writematrix(dCEALc,filename,'Sheet','Cea L ctrl');
writematrix(dCEARc,filename,'Sheet','Cea R ctrl');
writematrix(dPAGc,filename,'Sheet','PAG ctrl');
writematrix(dtime,filename,'Sheet','time');

