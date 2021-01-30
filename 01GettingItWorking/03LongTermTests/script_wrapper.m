clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters to vary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define all days, then remove those in the wrong part of the year
% Days = datenum(2002,1,242):1:datenum(2019,1,297);
DaysPerRun = 5;
Days = [datenum(2002,9,1):DaysPerRun:datenum(2019,12,31)];
[~,m,~] = datevec(Days);
Days = Days(m == 4 | m == 12);
%  %  Days(m >=  5 | m <= 11 ) = [];
%  %  [~,m,~] = datevec(Days);
%  %  Days(m < 4) = [];
clear m


%clobber?
NoClobber  = 1;

Count = 1;

for iDay=1:1:numel(Days);
  JobName = num2str(Days(iDay));
  NThreads = 16;
  
  Time = num2str(round(120/5*4));
  
  
  Text1 = ['#!/bin/bash\n## Name of the job\n#SBATCH --job-name=',JobName,'\n## Account to charge to\n#SBATCH --account=free\n\n'];
  Text2 = ['\n#SBATCH --time=',Time,':00\n## Number of node required and tasks per node\n'];
  Text3 = ['#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=',num2str(NThreads),'\n\n#SBATCH --output=log.%%j.out\n#SBATCH --error=log.%%j.err\n#SBATCH --partition=batch-all'];
  
  
  %Load Matlab environment
  Text4 = ['\n\n\nmodule load matlab\n\n'];
  
  
  Text5 = ['\nmatlab -nodisplay -r "cd /home/f/cw785/Matlab/20200303SH_SSWs/01GenerateGWData;maxNumCompThreads(',num2str(NThreads),');'];
  Commands = ['Settings.DayScale=',num2str(Days(iDay)),'+(0:1:',num2str(DaysPerRun-1),');process_data'];
  
  Text6 = [';exit"'];
  
  fid = fopen(['job',sprintf('%04d',Count),'_wrapper.txt'],'wt');
  fprintf(fid, Text1);
  fprintf(fid, Text2);
  fprintf(fid, Text3);
  fprintf(fid, Text4);
  fprintf(fid, Text5);
  fprintf(fid, Commands);
  fprintf(fid, Text6);
  fclose(fid);
  
  Count = Count+1;

end


%generate file to fire the whole lot off
fid = fopen(['fire_wrappers.sh'],'wt');
%  for i=1:1:Count-1;fprintf(fid,['sbatch --begin=now+7hour job',sprintf('%04d',i),'.txt\n']);end
for i=1:1:Count-1;
  %    fprintf(fid,['sbatch job',sprintf('%04d',i),'_wrapper.txt\n']);
  %    fprintf(fid,['sbatch --begin=now+6hour job',sprintf('%04d',i),'_wrapper.txt\n']);
  fprintf(fid,['sbatch --dependency singleton job',sprintf('%04d',i),'_wrapper.txt\n']);
  %    fprintf(fid,['sbatch --dependency singleton job',sprintf('%04d',i),'_wrapper.txt\n']);
  fprintf(fid,['rm job',sprintf('%04d',i),'_wrapper.txt\n']);
end
fprintf(fid,['rm fire_wrappers.sh\n']);
fclose(fid);

disp(['Written ',num2str(Count),' files (probably)'])

