%to detect spikes after initial processing using
%FearLearn3ROIregression_GCaMP

newmice=[1512 1514 1522 1531 1535 1537 1513 1521 1532 1533 1534 1536]

%import fear conditioning tone times calculated earlier
tonetimes=FreezingAnalysis
exp='FearLearn'

%new
for j=1:length(newmice)
m=newmice(j) 

    %load in processed photometry data
    datafile=sprintf('VgatFear2cLearn_%dtrace.xls',m)
    data = readmatrix(datafile);
    thzscore=sprintf('%s_zscoreTHMAD_%d',exp,m);
    filename=sprintf('%s_peakdataMAD_%d.xls',exp,m);
   
    %calculate threshold for spike detection based on MAD of baseline
    th=2*mad(data(1:300,4))
    
%find peaks that have a prominence above threshold
[P.pks,P.locs,P.widths,P.proms]=findpeaks(data(:,4),'MinPeakProminence',th,'Annotate','extents');

%find time of each spike
P.loc_photomtime=data(P.locs,1);
P.loc_min=data(P.locs,5);
P.loc_sec=P.loc_min*60
P.peakdata=[P.pks P.loc_photomtime P.loc_min P.loc_sec P.widths P.proms];

%plot trace with labelled peaks and save
figure
findpeaks(data(:,4),'MinPeakProminence',th,'Annotate','extents')
title('PAG')
savefig(thzscore)

%determine how many spikes/epoch

%find column with tone times for correct mouse
col=find(tonetimes(1,:)==m)

%determine which peaks are during which epoch
bin.bl=find(P.peakdata(:,4)<tonetimes(2,col));
bin.t1=find((P.peakdata(:,4)>tonetimes(2,col)) & P.peakdata(:,4)<tonetimes(3,col));
bin.iti1=find((P.peakdata(:,4)>tonetimes(3,col)) & P.peakdata(:,4)<tonetimes(4,col));
bin.t2=find((P.peakdata(:,4)>tonetimes(4,col)) & P.peakdata(:,4)<tonetimes(5,col));
bin.iti2=find((P.peakdata(:,4)>tonetimes(5,col)) & P.peakdata(:,4)<tonetimes(6,col));
bin.t3=find((P.peakdata(:,4)>tonetimes(6,col)) & P.peakdata(:,4)<tonetimes(7,col));
bin.iti3=find((P.peakdata(:,4)>tonetimes(7,col)) & P.peakdata(:,4)<tonetimes(8,col));
bin.t4=find((P.peakdata(:,4)>tonetimes(8,col)) & P.peakdata(:,4)<tonetimes(9,col));
bin.iti4=find((P.peakdata(:,4)>tonetimes(9,col)) & P.peakdata(:,4)<tonetimes(10,col));
bin.t5=find((P.peakdata(:,4)>tonetimes(10,col)) & P.peakdata(:,4)<tonetimes(11,col));
bin.consol=find(P.peakdata(:,4)>tonetimes(11,col)); 

%count the number of spikes during each epoch and put in first col
epochdata(1,:)=length(bin.bl);
epochdata(2,:)=length(bin.t1);
epochdata(3,:)=length(bin.iti1);
epochdata(4,:)=length(bin.t2);
epochdata(5,:)=length(bin.iti2);
epochdata(6,:)=length(bin.t3);
epochdata(7,:)=length(bin.iti3);
epochdata(8,:)=length(bin.t4);
epochdata(9,:)=length(bin.iti4);
epochdata(10,:)=length(bin.t5);
epochdata(11,:)=length(bin.consol);

%calculate the length of each epoch in seconds and put in second col
epochdata(1,2)=60;
for i=2:10
    epochdata(i,2)=(tonetimes(i+1,col)-tonetimes(i,col));
end
epochdata(11,2)=120;

%calculate frequency in third col
for i=1:11
epochdata(i,3)=(epochdata(i,1)/epochdata(i,2));
end

%normalize to baseline
epochdata(:,4)=(epochdata(:,3)/epochdata(1,3))

writematrix(P.peakdata,filename,'Sheet','spike details');
writematrix(epochdata,filename,'Sheet','epoch details');

clearvars -except oldmice newmice exp tonetimes
end
