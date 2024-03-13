dbLoc='/Users/dlevy/Dropbox (UCSF Department of Neurological Surgery)/ChangLab/General Patient Info';
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

maybeNoFile=[];
maybeWeirdTaskInfo=[]; 

notesFNames=readtable('~/Desktop/notesNames.txt','ReadVariableNames',false); % result of ls EC*/*rial* > notesNames.txt in Dropbox

notesSfxs=cellfun(@(x) split(x,{'_','.',' '}),notesFNames.Var2,'UniformOutput',false);

ecTypes=[];
notesTypes=[];
fTypes=[];
for s=1:length(notesSfxs)   


    thisName=notesSfxs{s};
    ecLbl=thisName(contains(thisName,'ec','IgnoreCase',true));
    notesLbl=thisName(contains(thisName,'ote','IgnoreCase',true));
    fType=thisName(end);

    ecTypes=[ecTypes;ecLbl];
    notesTypes=[notesTypes;notesLbl];
    fTypes=[fTypes;fType];

end

notesTypes=unique(notesTypes);
notesTypes(1)={'Trial_Notes'};
notesTypes(4)={'trial_notes'};
[~, a] = unique(lower(notesTypes));
notesTypes = notesTypes(a);

fTypes=unique(fTypes);

timitTypes=cell(1,1);
mochaTypes=cell(1,1);
burscTypes=cell(1,1);
interviewTypes=cell(1,1);

ecIDs=[cellstr(strcat('EC',string(1:309)))];

taskLbls={'TIMIT','MOCHA','BURSC','interview','conversation'};

taskInfo=cell(size(ecIDs,2),5,length(taskLbls));

allPossTaskNames=[];
allHtkPaths=cell(length(ecIDs),length(taskLbls));
allAninInfo=[];
for e=1:length(ecIDs)

    for n=1:size(notesTypes,1)
        for f=1:(size(fTypes,1))
            fName_1=sprintf('%s/%s/%s_%s.%s',dbLoc,ecIDs{e},ecIDs{e},notesTypes{n},fTypes{f});
            fName_2=sprintf('%s/%s/%s_%s.%s',dbLoc,ecIDs{e},notesTypes{n},ecIDs{e},fTypes{f});

            if ~exist(fName_1,'file') && ~exist(fName_2,'file')
                continue;
            else
                try
                    tbl=readtable(fName_1);
                    successName=fName_1;
                catch
                    try
                        tbl=readtable(fName_2);
                        successName=fName_2;
                    catch
                        tbl=table;
                        tbl.Task='';
                        successName='none';
                        maybeNoFile=[maybeNoFile; ecIDs(e)];
                    end
                end
            end

            try
                if ~isnan(tbl.Task{1})
                    allPossTaskNames=[allPossTaskNames; tbl.Task];
                    allAninInfo=[allAninInfo; tbl.Inputs];
                end
            catch
            end


            for t=1:length(taskLbls)
                switch t
                    case 1
                        alsoInclude={};
                        doNotInclude={'read','song','chinese','vocoded','reversed'};
                    case 2
                        alsoInclude={};
                        doNotInclude={'mime','covert','italian'};
                    case 3
                        alsoInclude={};
                        doNotInclude={};
                    case 4
                        alsoInclude={'frog','princess','pea','alice','prince','hare','tortoise','sleeping','stories'};
                        doNotInclude={'arabic','spanish'};
                    case 5
                        alsoInclude={};
                        doNotInclude={};
                end
                try
                    if sum(contains(tbl.Task,taskLbls{t},'IgnoreCase',true))>0 
                        thisEC=ecIDs{e};

                        allTaskCandidates=contains(tbl.Task,taskLbls{t},'IgnoreCase',true);
                        
                        for d=1:length(alsoInclude)
                            alsoUse=contains(tbl.Task,alsoInclude{d},'IgnoreCase',true);
                            allTaskCandidates=allTaskCandidates+alsoUse;
                        end

                        allTaskCandidates=allTaskCandidates>0;

                        exclude=zeros(size(allTaskCandidates));

                        for dn=1:length(doNotInclude)
                            isBad=contains(tbl.Task,doNotInclude{dn},'IgnoreCase',true);
                            exclude=exclude+isBad;
                        end

                        exclude=exclude~=0;

                        allTasks=logical(allTaskCandidates.*~exclude);

                        theseBlocks=tbl.Block(allTasks);
                        taskInfo{e,1,t}=thisEC;
                        taskInfo{e,2,t}=theseBlocks;
                        try
                            taskInfo{e,4,t}=tbl.Valid(allTasks);
                        catch
                            taskInfo{e,4,t}=repmat('unknown',size(theseBlocks));
                        end
                        try
                            taskInfo{e,5,t}=tbl.Inputs(allTasks);
                        catch
                            taskInfo{e,5,t}=repmat('unknown',size(theseBlocks));
                        end
                        htkPaths=cellstr(strcat(strcat(sprintf('/data_store1/prcsd_data/human/%s/%s_B',thisEC,thisEC),string(theseBlocks)),'/RawHTK'));
                        allHtkPaths(e,t)={htkPaths};
                        try
                            theseTypes=tbl.TaskType(contains(tbl.Task,taskLbls{t},'IgnoreCase',true));
                        catch
                            try
                                splitTypes=split(strrep(tbl.Task(contains(tbl.Task,taskLbls{t},'IgnoreCase',true)),' ',''),taskLbls{t});
                                theseTypes=str2double(splitTypes(:,2));
                            catch
                                theseTypes='checkManually';
                            end
                        end
                        taskInfo{e,3,t}=theseTypes;
                    end

                catch
                    maybeWeirdTaskInfo=[maybeWeirdTaskInfo; ecIDs(e)];
                end
            end

        end
    end
end

unique(allPossTaskNames)

unique(maybeNoFile) % people for whom it seems like we can't find a trial notes file
unique(maybeWeirdTaskInfo) % people for whom it seems like we can't get task info for some reason

%a(cellfun(@isempty,a))={'hi'};
%a(ismember(a,unique(maybeNoFile)'))

multTasks=~cellfun(@isempty,squeeze(taskInfo(:,1,:)));

sumAnyTasks=sum(any(multTasks,2));
sumAllTasks=sum(all(multTasks,2));
allCounts=nan(length(taskLbls),length(taskLbls));
for t1=1:length(taskLbls)
    for t2=1:length(taskLbls)
        thisCount=sum(all(multTasks(:,[t1 t2]),2));
        allCounts(t1,t2)=thisCount;
    end
end

imagesc(allCounts); 
xlabel(''); ylabel('');
xticks(1:length(taskLbls)); yticks(1:length(taskLbls));
xticklabels(taskLbls); yticklabels(taskLbls);
colorbar;

for t=1:length(taskLbls)
    numSubs=sum(~cellfun(@isempty,taskInfo(:,1,t)));
    numBlocks=sum(cellfun(@length,taskInfo(:,2,t)));
    fprintf('%s subs: %d\n',taskLbls{t},numSubs);
    fprintf('%s blocks: %d\n',taskLbls{t},numBlocks);
end

totalBlocks=0;
for t=1:length(taskLbls)
    totalBlocks=totalBlocks+sum(cellfun(@length,taskInfo(:,2,t)));
end

%% write out csv
elecLoc='/Users/dlevy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/imaging';

allPossLocs=[];
for s=1:length(ecIDs)
    thisSub=ecIDs(s);
    elecPath=sprintf('%s/subjects/%s/elecs/TDT_elecs_all.mat',elecLoc,thisSub{1});
    try
        load(elecPath);
    catch
        fprintf('no elec file for %s\n',thisSub{1});
    end
    allPossLocs=[allPossLocs; anatomy(:,4)];
end

elecCats=categories(categorical(allPossLocs));
elecCatCounts=countcats(categorical(allPossLocs));

[a,b]=sort(elecCatCounts,'descend');

elecCats=elecCats(b(1:10));

csvLine=table;
for t=1:length(taskLbls)

    thisTask=taskLbls(t);
    taskSubs=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),1,t);
    taskBlocks=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),2,t);
    taskTypes=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),3,t); 
    taskValidity=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),4,t); 
    taskInputs=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),5,t);

    for s=1:length(taskSubs)
        thisSub=taskSubs(s);
        elecPath=sprintf('%s/subjects/%s/elecs/TDT_elecs_all.mat',elecLoc,thisSub{1});
        try
            load(elecPath);
            thisSubAnat=anatomy(:,4);
            locCounts=zeros(size(elecCats));
            for c=1:length(elecCats)
                locCounts(c)=sum(strcmp(thisSubAnat,elecCats{c}));
            end
            xCoords=elecmatrix(:,1);
            if sum(xCoords(~isnan(xCoords))<0)/length(xCoords(~isnan(xCoords)))>0.9
                thisSubHem='left';
            elseif sum(xCoords(~isnan(xCoords))>0)/length(xCoords(~isnan(xCoords)))>0.9
                thisSubHem='right';
            else
                thisSubHem='bilateralOrUnknown';
            end
            elecTbl=array2table(locCounts','VariableNames',elecCats);
            elecTbl=addvars(elecTbl,{thisSubHem},'NewVariableNames','hem','Before',elecCats{1});
        catch
            fprintf('no elec file for %s\n',thisSub{1});
            elecTbl=array2table(nan(size(elecCats))','VariableNames',elecCats);
            thisSubHem='bilateralOrUnknown';
            elecTbl=addvars(elecTbl,{thisSubHem},'NewVariableNames','hem','Before',elecCats{1});
        end
        theseBlocks=taskBlocks{s};
        theseTypes=taskTypes{s};
        try
            thisValidity=cellstr(taskValidity{s});
        catch
            thisValidity=cellstr(repmat('unknown',size(taskValidity{s})));
        end
        try
            thisInput=cellstr(taskInputs{s});
        catch
            thisInput=cellstr(repmat('unknown',size(taskInputs{s})));
        end
        if isa(theseTypes,'char')
            theseTypes=nan(size(theseBlocks));
        end
        if isa(theseTypes,'cell') || all(isnan(theseTypes))
            theseTypes=str2double(theseTypes);
            if length(theseTypes)~=length(theseBlocks)
                theseTypes=nan(size(theseBlocks));
            end
        end
        if isa(theseBlocks,'double') && ~any(isnan(theseBlocks))
            for b=1:length(theseBlocks)
                csvLine=[csvLine; [table(thisSub,theseBlocks(b),thisTask,theseTypes(b),thisValidity(b),thisInput(b)) elecTbl]];
            end
        end
    end
end

csvLine.Properties.VariableNames(1:6)={'SubjectID','Block','Task','Set','Valid','Inputs'};

csvLine=csvLine(~any(isnan(csvLine.superiortemporal),2),:);

writetable(csvLine,'dataCSV_031224.csv');

%% visualization of electrodes

elecLoc='/Users/dlevy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/imaging';
allElecs_timit=[];
allElecs_mocha=[];
allElecs_bursc=[];
allElecs_int=[];
allElecs_conv=[];
noWarpedElecsFile=[];
for e=1:size(ecIDs,2)
    theseSets=taskInfo{e,3,t};
    if ~isempty(theseSets)
        try
            clear elecmatrix
            elecPath=sprintf('%s/subjects/%s/elecs/TDT_elecs_all_warped.mat',elecLoc,ecIDs{e});
            load(elecPath);
            switch t
                case 1
                    allElecs_timit=[allElecs_timit;elecmatrix];
                case 2
                    allElecs_mocha=[allElecs_mocha;elecmatrix];
                case 3
                    allElecs_bursc=[allElecs_bursc;elecmatrix];
                case 4
                    allElecs_int=[allElecs_int;elecmatrix];
                case 5
                    allElecs_conv=[allElecs_conv;elecmatrix];
            end
        catch
            noWarpedElecsFile=[noWarpedElecsFile; ecIDs(e)];
        end
    end
end

addpath(genpath('~/Documents/MATLAB/utils'));
addpath(genpath('spinningGIF'));

% timit

figure; title('timit');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs_timit=allElecs_timit(allElecs_timit(:,2,:)~=-128,:); % get rid of weird outliers (190 of these?)
timit_left=allElecs_timit(allElecs_timit(:,1,:)<0,:,:);
scatter3(timit_left(:,1)-2,timit_left(:,2),timit_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
timit_right=allElecs_timit(allElecs_timit(:,1,:)>0,:,:);
scatter3(timit_right(:,1)+2,timit_right(:,2),timit_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('timit_coverage.gif');
close all;

% mocha

figure; title('mocha');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs_mocha=allElecs_mocha(allElecs_mocha(:,2,:)~=-128,:); % get rid of weird outliers
mocha_left=allElecs_mocha(allElecs_mocha(:,1,:)<0,:,:);
scatter3(mocha_left(:,1)-2,mocha_left(:,2),mocha_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
mocha_right=allElecs_mocha(allElecs_mocha(:,1,:)>0,:,:);
scatter3(mocha_right(:,1)+2,mocha_right(:,2),mocha_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('mocha_coverage.gif');
close all;

% bursc

figure; title('bursc');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs_bursc=allElecs_bursc(allElecs_bursc(:,2,:)~=-128,:); % get rid of weird outliers
bursc_left=allElecs_bursc(allElecs_bursc(:,1,:)<0,:,:);
scatter3(bursc_left(:,1)-2,bursc_left(:,2),bursc_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
bursc_right=allElecs_bursc(allElecs_bursc(:,1,:)>0,:,:);
scatter3(bursc_right(:,1)+2,bursc_right(:,2),bursc_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('bursc_coverage.gif');
close all;

% interview

figure; title('interview');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs_int=allElecs_int(allElecs_int(:,2,:)~=-128,:); % get rid of weird outliers
int_left=allElecs_int(allElecs_int(:,1,:)<0,:,:);
scatter3(int_left(:,1)-2,int_left(:,2),int_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
int_right=allElecs_int(allElecs_int(:,1,:)>0,:,:);
scatter3(int_right(:,1)+2,int_right(:,2),int_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('interview_coverage.gif');
close all;

% conv

figure; title('conversation');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs_conv=allElecs_conv(allElecs_conv(:,2,:)~=-128,:); % get rid of weird outliers
conv_left=allElecs_conv(allElecs_conv(:,1,:)<0,:,:);
scatter3(conv_left(:,1)-2,conv_left(:,2),conv_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
conv_right=allElecs_conv(allElecs_conv(:,1,:)>0,:,:);
scatter3(conv_right(:,1)+2,conv_right(:,2),conv_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('conv_coverage.gif');
close all;

% all

figure; title('all');
mni.lh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_lh_pial.mat');
mni.rh = load('~/Documents/MATLAB/utils/cvs_avg35_Meshes/cvs_avg35_inMNI152_rh_pial.mat');
ctmr_gauss_plot(mni.lh.cortex,[0 0 0],0,'lh'); hold on; ctmr_gauss_plot(mni.rh.cortex,[0 0 0],0,'rh');
alpha(1);
view(270,0);
lightangle(270,0);
allElecs=[allElecs_timit; allElecs_mocha; allElecs_bursc; allElecs_int; allElecs_conv];
all_left=allElecs(allElecs(:,1,:)<0,:,:);
scatter3(all_left(:,1)-2,all_left(:,2),all_left(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
all_right=allElecs(allElecs(:,1,:)>0,:,:);
scatter3(all_right(:,1)+2,all_right(:,2),all_right(:,3)-2,'filled','MarkerFaceAlpha',0.7,'MarkerFaceColor','b');
spinningGIF('all_coverage.gif');
close all;

% A patients

aPts=readtable('~/Downloads/ec_pts.csv');
for t=1:length(taskLbls)
    thisTaskSum=sum(ismember(aPts.ec_id,taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),1,t)));
    fprintf('%s - %d subs with A and EC\n',taskLbls{t},thisTaskSum);
end


% demo

demo=readtable('ec_demo.csv');
demo.ID=strcat('EC',cellstr(string(demo.ID)));

for t=1:length(taskLbls)
    taskDemo=demo(ismember(demo.ID,taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),1,t)),:);
    countSex=countcats(categorical(taskDemo.Sex));
    fprintf('%s - %d F, %d M, %d unknown, meanAge %0.3f, stdAge %0.3f\n',taskLbls{t},countSex(1),countSex(2),sum(cellfun(@isempty,taskDemo.Sex)),nanmean(taskDemo.Age),nanstd(taskDemo.Age)); 
end
