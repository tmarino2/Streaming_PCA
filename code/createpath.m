%% Function creates hierarchical directory structure for RESULTS
function s=createpath(fullpath)
s=1;
list=textscan(fullpath,'%s','Delimiter','/');
currpath=pwd;
list2=textscan(currpath,'%s','Delimiter','/');
partialpath='';
minlen=min(length(list{1}),length(list2{1}));
for i=1:minlen
    partialpath=[partialpath,char(list{1}(i)),'/']; %#ok<AGROW>
    if(~strcmp(list{1}(i),list2{1}(i)))
        if(~exist(partialpath,'dir'))
            % unix(sprintf('mkdir ./%s',partialpath));
            [s,~,~]=mkdir(partialpath);
        end
    end
end
for i=minlen+1:length(list{1})
    partialpath=[partialpath,char(list{1}(i)),'/']; %#ok<AGROW>
    if(~exist(partialpath,'dir'))
        % unix(sprintf('mkdir ./%s',partialpath));
        [s,~,~]=mkdir(partialpath);
    end
end

end
