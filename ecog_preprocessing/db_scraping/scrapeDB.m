%% user input: type in filepaths

% local location of General Patient Info folder on Dropbox
dbLoc="/Users/dlevy/Dropbox (UCSF Department of Neurological Surgery)/ChangLab/General Patient Info";

% local location of Neuropixels Management spreadsheet
npLoc='Neuropixels management - Subjects.csv';

% local location of /data_store2/imaging/ directory (for me I mount this via Mountain Duck)
elecLoc='/Users/dlevy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/imaging';

%% scraping code

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

maybeNoFile=[];
maybeWeirdTaskInfo=[]; 
maybeWeirdQuestionInfo=[]; 

cmd=sprintf("cd '%s'; ls EC*/*rial* > ~/Desktop/notesNames.txt",dbLoc);
system(cmd);
notesFNames=readtable('~/Desktop/notesNames.txt','ReadVariableNames',false); 
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

npTbl=readtable(npLoc); % downloaded from Google sheets; cannot read directly due to restricted sheet

taskLbls={'TIMIT','MOCHA','BURSC','interview','conversation'};

taskInfo=cell(size(ecIDs,2),12,length(taskLbls));

allPossTaskNames=[];
allHtkPaths=cell(length(ecIDs),length(taskLbls));
allAninInfo=[];
for e=1:length(ecIDs)

    npID=npTbl.Var1(strcmp(npTbl.EC,ecIDs{e}));
    aID=npTbl.A(strcmp(npTbl.EC,ecIDs{e}));

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
                qTbl=readtable(successName,'Sheet','Questionnaire',Range="A1:M60");
                hand=qTbl.Var2(strcmp(qTbl.Var1,'Handedness'));
                imped=qTbl.Var2(strcmp(qTbl.Var1,'Impediments'));
                music={[qTbl.Var2(strcmp(qTbl.Var1,'Instrument / Singing 1')) qTbl.Var3(strcmp(qTbl.Var1,'Instrument / Singing 1'))]};
                native={[qTbl.Var2(strcmp(qTbl.Var1,'Native')) qTbl.Var3(strcmp(qTbl.Var1,'Native')) qTbl.Var4(strcmp(qTbl.Var1,'Native')) qTbl.Var5(strcmp(qTbl.Var1,'Native'))]};
                secondLang={[qTbl.Var2(strcmp(qTbl.Var1,'Second')) qTbl.Var3(strcmp(qTbl.Var1,'Second')) qTbl.Var4(strcmp(qTbl.Var1,'Second')) qTbl.Var5(strcmp(qTbl.Var1,'Second'))]};
            catch
                maybeWeirdQuestionInfo=[maybeWeirdQuestionInfo; ecIDs(e)];
                hand={'nan'};
                imped={'nan'};
                music={'nan'};
                native={'nan'};
                secondLang={'nan'};
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
                        thisNP=npID;
                        thisA=aID;

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
                        taskInfo{e,6,t}=thisNP;
                        taskInfo{e,7,t}=thisA;
                        taskInfo{e,8,t}=hand;
                        taskInfo{e,9,t}=imped;
                        taskInfo{e,10,t}=music;
                        taskInfo{e,11,t}=native;
                        taskInfo{e,12,t}=secondLang;
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
                        %htkPaths=cellstr(strcat(strcat(sprintf('/data_store1/prcsd_data/human/%s/%s_B',thisEC,thisEC),string(theseBlocks)),'/RawHTK'));
                        %allHtkPaths(e,t)={htkPaths};
                        try
                            theseTypes=tbl.TaskType(allTasks);
                            if iscell(theseTypes)
                                theseTypes=str2double(theseTypes);
                            elseif all(isnan(theseTypes))
                                theseTypes=tbl.TaskSubType(allTasks);
                            end
                        catch
                            try
                                splitTypes=split(strrep(tbl.Task(allTasks),' ',''),taskLbls{t});
                                theseTypes=str2double(splitTypes(:,2));
                            catch
                                try
                                    theseTypes=str2double(strrep((strrep(lower(tbl.Task(allTasks)),lower(taskLbls{t}),'')),' ',''));
                                catch
                                    theseTypes='checkManually';
                                end
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
unique(maybeWeirdQuestionInfo) % people for whom it seems like we can't get questionairre info for some reason

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
    taskNPSubs=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),6,t);
    taskASubs=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),7,t);
    taskHand=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),8,t);
    taskImped=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),9,t);
    taskMusic=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),10,t);
    taskNative=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),11,t);
    taskSecond=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),12,t);
    taskBlocks=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),2,t);
    taskTypes=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),3,t); 
    taskValidity=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),4,t); 
    taskInputs=taskInfo(~cellfun(@isempty,taskInfo(:,1,t)),5,t);

    for s=1:length(taskSubs)
        dataStruct=struct;
        dataStruct.thisSub=taskSubs(s);
        dataStruct.thisNP=taskNPSubs{s};
        dataStruct.thisA=taskASubs{s};
        dataStruct.thisHand=taskHand{s};
        dataStruct.thisImped=taskImped{s};
        dataStruct.thisNative=taskNative{s};
        dataStruct.thisSecond=taskSecond{s};
        dataStruct.thisMusic=taskMusic{s};
        fn=fieldnames(dataStruct);
        for i=1:numel(fn)
            switch fn{i}
                case 'thisNative'
                    try
                        dataStruct.thisNative_lang=dataStruct.thisNative{1}(1);
                        dataStruct.thisNative_prof=dataStruct.thisNative{1}(2);
                        dataStruct.thisNative_freq=dataStruct.thisNative{1}(3);
                        dataStruct.thisNative_age=dataStruct.thisNative{1}(4);
                    catch
                        dataStruct.thisNative_lang={'na'};
                        dataStruct.thisNative_prof={'na'};
                        dataStruct.thisNative_freq={'na'};
                        dataStruct.thisNative_age={'na'};
                    end
                case 'thisSecond'
                    try
                        dataStruct.thisSecond_lang=dataStruct.thisSecond{1}(1);
                        dataStruct.thisSecond_prof=dataStruct.thisSecond{1}(2);
                        dataStruct.thisSecond_freq=dataStruct.thisSecond{1}(3);
                        dataStruct.thisSecond_age=dataStruct.thisSecond{1}(4);
                    catch
                        dataStruct.thisSecond_lang={'na'};
                        dataStruct.thisSecond_prof={'na'};
                        dataStruct.thisSecond_freq={'na'};
                        dataStruct.thisSecond_age={'na'};
                    end
                case 'thisMusic'
                    try
                        dataStruct.thisMusic_inst=dataStruct.thisMusic{1}(1);
                        dataStruct.thisMusic_dur=dataStruct.thisMusic{1}(2);
                        if strcmp(dataStruct.thisMusic_inst,'n')
                            dataStruct.thisMusic_inst={'na'};
                        end
                        if strcmp(dataStruct.thisMusic_dur,'a')
                            dataStruct.thisMusic_dur={'na'};
                        end
                    catch
                        dataStruct.thisMusic_inst={'na'};
                        dataStruct.thisMusic_dur={'na'};
                    end
                otherwise
                    if isempty(dataStruct.(fn{i}))
                        dataStruct.(fn{i})={'na'};
                    end
            end
        end
        elecPath=sprintf('%s/subjects/%s/elecs/TDT_elecs_all.mat',elecLoc,dataStruct.thisSub{1});
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
            fprintf('no elec file for %s\n',dataStruct.thisSub{1});
            elecTbl=array2table(nan(size(elecCats))','VariableNames',elecCats);
            thisSubHem='bilateralOrUnknown';
            elecTbl=addvars(elecTbl,{thisSubHem},'NewVariableNames','hem','Before',elecCats{1});
        end
        theseBlocks=taskBlocks{s};
        theseTypes=taskTypes{s};
        try
            thisValidity=cellstr(taskValidity{s});
        catch
            thisValidity=cellstr(repmat('unknown',size(theseBlocks)));
        end
        try
            thisInput=cellstr(taskInputs{s});
        catch
            thisInput=cellstr(repmat('unknown',size(theseBlocks)));
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
                csvLine=[csvLine; [table(dataStruct.thisSub,theseBlocks(b),thisTask,theseTypes(b),thisValidity(b),thisInput(b)),dataStruct.thisNP,dataStruct.thisA,dataStruct.thisHand,dataStruct.thisImped,dataStruct.thisNative_lang,dataStruct.thisNative_prof,dataStruct.thisSecond_lang,dataStruct.thisSecond_prof,dataStruct.thisMusic_inst,dataStruct.thisMusic_dur,elecTbl]];
            end
        end
    end
end

csvLine.Properties.VariableNames(1:16)={'SubjectID','Block','Task','Set','Valid','Inputs','NPID','AID','Handedness','SpeechImped','NativeLang','NativeLangProf','SecondLang','SecondLangProf','MusicInstrument','MusicDur'};

csvLine=csvLine(~any(isnan(csvLine.superiortemporal),2),:);

writetable(csvLine,sprintf('dataCSV_%s.csv',datetime('now','Format','uuuu-MM-dd')));