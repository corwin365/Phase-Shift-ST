clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% balena operation settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name of the programme to be used
ToCall.ProgName = 'lambdaz_analysis_fake_systematic';

%directory it will run in
ToCall.WorkingDirectory = '/home/f/cw785/Matlab/run_systematic_analysis/';

Queue = 'batch-devel'  ;%options are batch-short, batch-sky, batch-all, batch-devel
RunTime = 360; %this will be overriden with 15 min if 'batch-devel' queue is used, and only one job fired


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual function calls
%
%this part will be different for every job. Aim is to create a
%cell array with each entry representing an individual job
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%variable combinations to loop over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%cases to plot: {Date,Granules,XT row, AT row}
Cases.Case1 = {datenum(2008,1,127),1,67,45}; %we're using this one to get background, but the wave will be fake

%define the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%first tuning pass - many parameters fixed or artificially narrow to constrain variable space:
PassName  = 'tuning1';
c1        = [0.125,0.25,0.5];
c2        = [0.125,0.25,0.5];
NPeaks    = 1;%1:1:3;
Steps     = 1:1:3;
Weights   = 0;%[0,1];
RealNoise = 0;%[0,1];
Amplitude = 1;%0:1:3;
Lz        = 20:20:60;
Lx        = 250:250:1000;
Ly        = 250:250:1000;
RotX      = 0; 
RotY      = 0;
RotZ      = 0;


%make all the combinations so we can do a simple 1d loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllCombos = NaN(1,14); %the 1 will expand out
k = 0;
for iA = 1:1:numel(c1);
  for iB= 1:1:numel(c2);
    if c2(iB) > c1(iA); continue; end %fine must be finer than coarse!
    for iC = 1:1:numel(NPeaks)
      for iD = 1:1:numel(Steps)
        for iE = 1:1:numel(Weights);
          for iF = 1:1:numel(RealNoise);
            for iG=1:1:numel(Amplitude)
              for iH=1:1:numel(Lz)
                for iI=1:1:numel(Lx)
                  for iJ=1:1:numel(Ly)
                    for iK=1:1:numel(RotX)
                      for iL=1:1:numel(RotX)
                        for iM=1:1:numel(RotX)
                          for iCase = 1:1:numel(Cases)
                            k = k+1;
                            AllCombos(k,:) = [c1(iA),c2(iB),NPeaks(iC),Steps(iD),Weights(iE),iCase,RealNoise(iF),Amplitude(iG),Lz(iH),Lx(iI),Ly(iJ),RotX(iK),RotY(iL),RotZ(iM)];
                          end
                        end
                      end
                    end
                  end
                end
              end
            end;
          end
        end
      end
    end
  end
end
clear iA iB iC iD iE iF iCase iG iH iI iJ iK iL iM
clear c1 c2 NPeaks Steps Weights RealNoise


%we now have a list of all the individual combos. let's split them into groups
%to submit to the task queue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CombosPerJob = 50;

for iJob = 1%:1:ceil(k./CombosPerJob);
  
  %prepare call variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  InRange = 1+ (((iJob-1)*CombosPerJob) : 1 : ((iJob*CombosPerJob)-1));
  if max(InRange) > k; InRange(InRange >k) = []; end
  
  InfoFile = ['input_',sprintf('%06d',iJob),'.mat'];
  
  ComboData = AllCombos(InRange,:);
  SaveFile  = ['/beegfs/scratch/user/f/cw785/Data/corwin/2dp1/',PassName,'/variations_',sprintf('%06d',iJob),'.mat'];
  save(InfoFile,'ComboData','SaveFile','Cases');

  ThisCall = ['InfoFile =''',InfoFile,''';'];
       
  clear InRange InfoFile ComboData SaveFile

  %generate command
  %%%%%%%%%%%%%%%%%%
  
  JobNames{iJob} = ['2dp1_',sprintf('%03d',iJob)];
  Calls{iJob} = ThisCall;
end; clear iJob k PassName


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual job scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%special settings for certain queues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%number of cores
switch Queue
  case 'batch-sky';   Cores = 24;
  otherwise;    Cores = 16;
end

%quality of service
switch Queue
  case 'batch-devel'; QOS = '#SBATCH --qos=devel\n';
  otherwise;    QOS = ''; 
end

%max runtime
switch Queue
  case 'batch-devel'; RunTime = 15;
  otherwise; %already set, don't override
end

%max jobs
switch Queue
  case 'batch-devel'; MaxJobs = 1;
  otherwise;          MaxJobs = numel(Calls);
end

%ok, make the scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Count = 1;
for iScript=1:1:MaxJobs
  
  %general setup for Matlab balena jobs the random pause at the end 
  %is to spread out parpool calls if used in the scripts, as multiple calls
  %in rapid succession can cause scripts to crash when creating the pool 
 Start = ['#!/bin/bash\n', ...
         '#SBATCH --job-name=',JobNames{iScript},'\n', ...
         '#SBATCH --account=free\n', ...
         '#SBATCH --time=',num2str(RunTime),':00\n', ...
         '#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=',num2str(Cores),'\n', ...
         QOS, ...
         '#SBATCH --output=log.%%j.out\n', ...
         '#SBATCH --error=log.%%j.err\n', ...
         '#SBATCH --partition=',Queue,'\n', ...
         'module load matlab\n',...
         'matlab -nodisplay -r "cd ',ToCall.WorkingDirectory,';', ...
         'pause(',num2str(randi(30)),');'];
       End = ';exit;"';

 
  Commands = [Calls{iScript},';',ToCall.ProgName];
  
  %create scripts
  fid = fopen(['job',sprintf('%04d',Count),'_wrapper.txt'],'wt');
  fprintf(fid, Start);fprintf(fid, Commands);fprintf(fid, End);
  fclose(fid);
  Count = Count+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate master script to call the above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%generate file to fire the whole lot off
fid = fopen(['fire_wrappers.sh'],'wt');
for i=1:1:Count-1;
  fprintf(fid,['sbatch job',sprintf('%04d',i),'_wrapper.txt\n']);
  fprintf(fid,['rm job',sprintf('%04d',i),'_wrapper.txt\n']);
end
fprintf(fid,['rm fire_wrappers.sh\n']);
fclose(fid);

disp(['Written ',num2str(Count),' files (probably)'])

