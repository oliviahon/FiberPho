
%to extract behavior bout start and stop times from simba output

%will have an error with start_time if you have a trial where no behavior
%is detected at all

trials=[163 164 165 166 168 170 222 223 224 225 226 228 229 230 231 232 233 234]; 

for j=1:length(trials)
    
t=trials(j) 

filename=sprintf('Trial%dbouts.xls',t) %defines output filename
datafile=sprintf('Trial%d.csv',t)%reads in datafile


freezdata=readmatrix(datafile,'Range','A:A'); %col 1 is frame number
freezdata(:,2)=readmatrix(datafile,'Range','IN:IN'); %col 2 is true/false for behavior of interest
freezdata(1,:)=[];

freezdata(:,3)= freezdata(:,1)/30; %col 3 converts frames to time in S
freezdata(:,4)=freezdata(:,1)/1800; %col 4 converts frames to time in min

freezdata(1:1800,:)=[]

down=downsample(freezdata(:,4),3); 
down(:,2)=downsample(freezdata(:,2),3); 
down(:,3)=down(:,2)*5


down(down == 0) = NaN
down(1:300)=[]


freezdata_log=logical(freezdata(:,2))'; %makes logical array of classifier output

starts=strfind(freezdata_log,[0 1])'; %gives index number of when bout starts
ends= strfind(freezdata_log,[1 0])'; %gives index of bout end

%give time in min of transition start and end based on bout start and end
%indicies determined above
for i=1:length(starts)
    TS=(starts(i)); %time start
    starts_time(i,1)=freezdata(TS,4);%finds start time in min
end

for i=1:length(ends)
    TE=ends(i);%time end
    ends_time(i,1)=freezdata(TE,4); %finds bout end in min
end

%if you have an uneven number of starts and stops, this makes start_time
%and ends_time the same length
if length(starts)>length(ends)
    n=length(starts)
    ends_time(n,1)=NaN
end

bouts=[starts_time,ends_time]; %compiles into single matrix with two columns. start time in col 1, end in col 2
bouts(:,3)=bouts(:,2)-bouts(:,1); %calculates bout length (in minutes)
bouts(:,4)=bouts(:,3)*60; %converts bout length to seconds

writematrix(bouts,filename)
clearvars -except trials 
end