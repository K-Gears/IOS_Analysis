%% motion_cont.m
% by Bingxing Huo, 2013
% This script finds segments of continuous motion/rest
% 
% Input:
%   - imp_bin: output from velocity_proc_auto.m
%   - Fs: sampling frequency
%
% Output:
%   - T_run: 2-by-nr matrix. First row contains indices of start of
%   movement; second row contains indices of end of movements
%   - T_stand: 2-by-ns matrix. First row contains indices of start of
%   rest; second row contains indices of end of rest.
%%
function [T_run,T_stand,new_T_run,run_frac]=motion_cont_SleepScore(imp_bin,Fs,T_seg,T_beg)
runcont(1:length(imp_bin))=imp_bin;
L=length(imp_bin);
rind=find(imp_bin==1);
T_run=[];
T_stand=[];
new_T_run=[];
for k=2:length(rind)
    if le((rind(k)-rind(k-1)),(Fs))
        runcont(rind(k-1):rind(k))=1; %fuse together any running events that are less than 1sec apart
    end
end
nr=1;
ns=1;
T_run=zeros(2,1);
T_stand=zeros(2,1);

for k=2:L
    if runcont(k-1)==0
        if runcont(k)==1
            T_run(1,nr)=k;
            T_stand(2,ns)=k-1;
            ns=ns+1;
        end
    elseif runcont(k-1)==1
        if runcont(k)==0
            T_run(2,nr)=k-1;
            nr=nr+1;
            T_stand(1,ns)=k;
        end
    end
end

if T_run(1,1)==0
    T_run(1,1)=1;
end
if T_stand(1,1)==0
    T_stand(1,1)=1;
end
if T_run(2,end)==0
    T_run(2,end)=L;
end
if T_stand(2,end)==0 
    T_stand(2,end)=L;
end

run_length=T_run(2,:)-T_run(1,:)+1;
run_frac=sum(run_length)/length(imp_bin);

ShortRun=run_length<T_seg*Fs; %find running events less than a desired length
temp_T_run=T_run;
temp_T_run(:,ShortRun)=[]; %hold on to events of desired duration

if ~isempty(temp_T_run)
if (temp_T_run(1,1)==1)
    temp_T_run=temp_T_run(:,(2:end));
end

run_lengths=temp_T_run(2,:)-temp_T_run(1,:);

if (max(run_lengths)>(Fs*T_seg)>0)
    j=0;
    for k=1:size(temp_T_run,2)
        stand_ind=find(T_stand(2,:)==(temp_T_run(1,k)-1));
        if gt((T_stand(2,stand_ind)-T_stand(1,stand_ind)),(T_beg*Fs)) %determine if there is sufficient rest before run onset
            if gt((temp_T_run(2,k)-temp_T_run(1,k)),(T_seg*Fs)) 
                j=j+1;
                new_T_run(2,j)=temp_T_run(2,k);
                new_T_run(1,j)=temp_T_run(1,k);
            end
        else
%             if gt(T_stand(2,k)-T_stand(1,k),(2*Fs))
%                 if k>1
%                     if lt(temp_T_run(2,(k-1))-temp_T_run(1,(k-1)),Fs)
%                         j=j+1;
%                         new_T_run(2,j)=temp_T_run(2,k);
%                         new_T_run(1,j)=temp_T_run(1,k);
%                     end
%                 else
%                     j=j+1;
%                     new_T_run(2,j)=temp_T_run(2,k);
%                     new_T_run(1,j)=temp_T_run(1,k);
%                 end
%             end
        end
    end
end
end



% Yurong Gao
% 2014/7/16
% get the period of standing time with T_before ahead and T_end behind
% without running

% T_run: matrix 2Xn, start and end of the running time
% T_stand: maxtrix 2Xn, start and end of the standing time
% new_T_stand: the new matrix of the standing priod of T_before no running
% and T_end no running
% if (isempty(T_stand))      % if the T_stand is empty new_T_stand is empty
%     new_T_stand=[];
% else
%     if (T_stand(1,1)==1)              % delete the first T_stand period if it start from the first frame
%         T_stand=T_stand(:,2:end);
%     end
%     
%     if (isempty(T_stand))            % after the delete of the first T_stand if the T_stand is empty, the new_T_stand is
%         new_T_stand=[];
%     else  
%         
%         if (T_run(1,end)>T_stand(2,end))    % if the last period is running,delete it, so there will be run, stand, run, stand
%             T_run=T_run(:,1:end-1);
%         end
%         
%         run_lengths=T_run(2,:)-T_run(1,:)+1;
%         big_T_run=T_run;
%         big_T_run(:,find(run_lengths<2))=[]; % delete the 1frame running
%         con_T_stand=[];
%         
%         if (isempty(big_T_run))          % if there is no running, the new_T_stand is the same as T_stand
%             new_T_stand=T_stand;
%             
%         else           
%             for i=1:size(big_T_run,2)-1                    
%                 con_T_stand(2,i)=big_T_run(1,i+1)-1;
%                 con_T_stand(1,i)=big_T_run(2,i)+1;
%             end
%             con_T_stand(1,i)=big_T_run(2,end)+1;
%             con_T_stand(2,i)=T_stand(2,end);
%             new_T_stand=[];
%             j=0;
%             stand_lengths=con_T_stand(2,:)-con_T_stand(1,:);
%             if (sum(stand_lengths>Fs*(T_before+T_end))>0)
%                 for i=1:size(con_T_stand,2)
%                     if (con_T_stand(2,i)-con_T_stand(1,i)>=(T_before+T_end)*Fs)
%                         j=j+1;
%                         new_T_stand(2,j)=con_T_stand(2,i)-T_end*Fs;
%                         new_T_stand(1,j)=con_T_stand(1,i)+Fs*T_before;
%                     end
%                 end
%             end
%         end
%     end
% end
end


