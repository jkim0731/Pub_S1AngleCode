
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')

%% Fig 1A - Task
%% Fig S1 - Trial and session structure

%%
%% Fig 1B Learning curves
%%

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[1:18],[1:7,15],[1:20],[1:16], [1:24], [1:31], [1:21], [1:30], [1:20], [1:21], [1:30], [1:13]};

% initialization
pcbysessions = cell(length(mice),1); % percent correct by sessions
databytrials = cell(length(mice),1); % data by trials
pcbytrials = cell(length(mice),1); % percent correct by trials
sessionDurations = cell(length(mice),1);

for mi = 1 : length(mice)
    mouse = sprintf('JK%03d',mice(mi));
    fn = sprintf('behavior_%s.mat',mouse);
    load([baseDir, 'Behavior\', fn])
    
    pcbysessions{mi} = zeros(length(sessions{mi}),1);
    databytrials{mi} = [];
    sessionDurations{mi} = zeros(length(sessions{mi}),1);
    for si = 1 : length(sessions{mi})
        bi = find(cellfun(@(x) strcmp(x.sessionName, sprintf('S%02d',sessions{mi}(si))), b));
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{bi}.trials)); % 'oo' trials are catch trials (where the object was too far away so that mice could not touch it)
        pcbysessions{mi}(si) = sum(b{bi}.hitTrialInds(notooind))/(sum(b{bi}.hitTrialInds(notooind)) + sum(b{bi}.faTrialInds(notooind))) * 100;
        tempCorrects = b{bi}.trialCorrects;
        tempCorrects(tempCorrects==-1) = []; % remove misses
        databytrials{mi} = [databytrials{mi}, tempCorrects];
        sessionDurations{mi}(si) = length(tempCorrects);        
    end
    pcbytrials{mi} = zeros(length(databytrials{mi})-99,1);
    for ti = 1 : length(databytrials{mi})-99
        pcbytrials{mi}(ti) = sum(databytrials{mi}(ti:ti+99));
    end
end

%% Calculate mean +- sem of # of sessions to expert
learnedInd = [1:4,7,9];
threshold = 75;
expertSessions = zeros(length(learnedInd),1);
for i = 1 : length(learnedInd)
    temp = pcbysessions{learnedInd(i)};
    si = find(temp > 75);
    diffsi = diff(si);
    sumdiffsi = zeros(length(diffsi)-1,1);
    for j = 1 : length(sumdiffsi)
        sumdiffsi(j) = sum(diffsi(j:j+1));
    end
    finalInd = find(sumdiffsi==2,1,'first');
    expertSessions(i) = si(finalInd);
end

mean(expertSessions)
sem(expertSessions)

%%
figure, hold all,
for mi = 1 : length(mice)
    if length(find(pcbysessions{mi} > 75 )) > 2
        plot(pcbysessions{mi}, 'k-', 'linewidth', 1)
    else
        plot(pcbysessions{mi}, 'color',[0.7 0.7 0.7], 'linewidth', 1)
    end
end
plot(0:35, ones(36,1)*75, 'r--', 'linewidth', 1)
xlabel('Sessions'), ylabel('% Correct')
set(gca,'fontsize',12, 'fontname', 'Arial')

errorbar(mean(expertSessions), 95, sem(expertSessions), 'horizontal', 'ko')


%%
%% Fig 1C. Left lick probability in 7-angle test sessions, before (left) and after (right) learning
%%
angles = 45:15:135;
learned = [25,27,30,36,39,52];
sessionsNaive = [4,3,3,1,1,3];
sessionsExpert = [19,10,21,17,23,21];
plickLBefore = zeros(length(learned),length(angles));
plickLAfter = zeros(length(learned),length(angles));

for mi = 1 : length(learned)
    mouse = sprintf('JK%03d',learned(mi));
    fn = sprintf('behavior_%s.mat',mouse);
    load([baseDir, 'Behavior\', fn])
    
    bindNaive = find(cellfun(@(x) strcmpi(x.sessionName, sprintf('S%02d',sessionsNaive(mi))), b));
    bNaive = b{bindNaive};
    
    for ai = 1 : length(angles)
        angleTrials = find(cellfun(@(x) x.servoAngle == angles(ai), bNaive.trials));
        leftLicks = find(cellfun(@(x) strcmp(x.choice, 'l'), bNaive.trials(angleTrials)));
        rightLicks = find(cellfun(@(x) strcmp(x.choice, 'r'), bNaive.trials(angleTrials)));
        plickLBefore(mi,ai) = length(leftLicks) / (length(leftLicks) + length(rightLicks));
    end
    
    bindExpert = find(cellfun(@(x) strcmpi(x.sessionName, sprintf('S%02d',sessionsExpert(mi))), b));
    bExpert = b{bindExpert};
    
    for ai = 1 : length(angles)
        angleTrials = find(cellfun(@(x) x.servoAngle == angles(ai), bExpert.trials));
        leftLicks = find(cellfun(@(x) strcmp(x.choice, 'l'), bExpert.trials(angleTrials)));
        rightLicks = find(cellfun(@(x) strcmp(x.choice, 'r'), bExpert.trials(angleTrials)));
        plickLAfter(mi,ai) = length(leftLicks) / (length(leftLicks) + length(rightLicks));
    end
end
%%
figure('units','inch','pos',[2 2 7 3]),
subplot(121), hold on
for mi = 1 : length(learned)    
    plot(angles, plickLBefore(mi,:), 'k-', 'linewidth', 1)
end
errorbar(angles, mean(plickLBefore),sem(plickLBefore), 'r-', 'linewidth', 2)
ylabel('Left lick probability'), xlabel('Angle (\circ)'), xticks(angles), title('Before learning'), ylim([0 1])
xlim([40 140]), xticks(45:30:135)


subplot(122), hold on
for mi = 1 : length(learned)    
    plot(angles, plickLAfter(mi,:), 'k-', 'linewidth', 1)
end
errorbar(angles, mean(plickLAfter),sem(plickLAfter), 'r-', 'linewidth', 2)
ylabel('Left lick probability'), xlabel('Angle (\circ)'), xticks(angles), title('After learning'), ylim([0 1])
xlim([40 140]), xticks(45:30:135)


%%
%% Fig 1D. Choice dependence on angle difference, in Expert 7-angle test sessions. 
%%

%% Lick rate difference vs Angle difference 
lickRateDiff = zeros(length(learned), length(angles)-1);
c = nchoosek(1:length(angles), 2);
for mi = 1 : length(learned)
    for ai = 1 : length(angles)-1
        tempInd = find(c(:,2) - c(:,1) == ai);
        lickRateDiff(mi,ai) = mean(plickLAfter(mi,c(tempInd,2)) - plickLAfter(mi,c(tempInd,1)),2);
    end
end

mean(lickRateDiff)
std(lickRateDiff)

%%
figure,
errorbar(15:15:90,mean(lickRateDiff), std(lickRateDiff), 'ko-')
ylim([0 1])
xticks(15:15:90)
xlim([10 95])
xlabel('Angle difference (\circ)')
ylabel('Left lick probability difference')
set(gca, 'fontsize', 12, 'fontname', 'Arial', 'box', 'off')


%% statistical test
[~,p] = ttest(lickRateDiff(:,1))



%%
%% Fig 1E. Performance of expert mice at different experimental stages.
%%
learned = [25,27,30,36,39,52];
discreteCorrect = zeros(length(learned),5);

for mi = 1 : length(learned)
    mouse = sprintf('JK%03d',learned(mi));
    fn = sprintf('behavior_%s.mat',mouse);
    load([baseDir, 'Behavior\', fn])
    
    binds = find(cellfun(@(x) strcmpi(x.taskTarget, 'Angle-Discrete'), b));
    notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(1)}.trials));
    discreteCorrect(mi,1) = sum(b{binds(1)}.hitTrialInds(notooind))/(sum(b{binds(1)}.hitTrialInds(notooind)) + sum(b{binds(1)}.faTrialInds(notooind))) * 100;
    notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(2)}.trials));
    discreteCorrect(mi,3) = sum(b{binds(2)}.hitTrialInds(notooind))/(sum(b{binds(2)}.hitTrialInds(notooind)) + sum(b{binds(2)}.faTrialInds(notooind))) * 100;
    
    binds = find(cellfun(@(x) strcmpi(x.distractor, 'Discrete'), b)); 
    if learned(mi) == 52 % for mice 52-56, I had radial distractor sessions (distance from the face) before learning as well, so have to select the second session for Expert sessions
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(2)}.trials));
        discreteCorrect(mi,4) = sum(b{binds(2)}.hitTrialInds(notooind))/(sum(b{binds(2)}.hitTrialInds(notooind)) + sum(b{binds(2)}.faTrialInds(notooind))) * 100;
    else
        notooind = find(cellfun(@(x) ~strcmp(x.trialType, 'oo'), b{binds(1)}.trials));
        discreteCorrect(mi,4) = sum(b{binds(1)}.hitTrialInds(notooind))/(sum(b{binds(1)}.hitTrialInds(notooind)) + sum(b{binds(1)}.faTrialInds(notooind))) * 100;
    end
    
    if mi == 3
        discreteCorrect(mi,5) = sum(b{end-2}.hitTrialInds)/(sum(b{end-2}.hitTrialInds) + sum(b{end-2}.faTrialInds)) * 100;
    elseif mi == 6
        discreteCorrect(mi,5) = NaN; % One mouse was lost before 'no whisker' session
    else
        discreteCorrect(mi,5) = sum(b{end-1}.hitTrialInds)/(sum(b{end-1}.hitTrialInds) + sum(b{end-1}.faTrialInds)) * 100;
    end
    mind = find(mice == learned(mi));
    tempPC = sort(pcbysessions{mind}, 'descend');
    discreteCorrect(mi,2) = mean(tempPC(1:3)); % For Expert 2-angle, take the average of 3 sessions before Expert 7-angle sessions
end

nanmean(discreteCorrect)
sem(discreteCorrect)

%%
figure, hold all
for mi = 1 : length(learned)
    plot(discreteCorrect(mi,1:5), 'k.', 'markersize',20)
end

xticks([1:5]), xlim([0.5 5.5]), ylabel('% Correct')
labels = ({'Naive 7-angle', 'Expert 2-angle', 'Expert 7-angle','Radial distractor', 'No whisker'});
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
XTickLabel = labels;
set(gca,'fontsize',12, 'xticklabel', labels)





%%
%% Fig S1B. Training structure.
%%

% Manual curation of # of sessions in each training phase
learner = {[1,3,1,10,4,1,1,3],    [1,2,1,1,3,3,1,3],    [2,2,1,14,3,1,1,3],   [1,0,1,12,3,1,1,3], [1,0,1,14,6,2,2,3], [1,2,1,13,3,4,1,0]};
nonlearner = {[2,6,1,13,0,0,0,0],    [1,1,1,29,0,0,0,0],   [1,2,1,27,0,0,0,0],   [1,2,1,18,0,0,0,0], [1,2,1,24,0,0,0,0], [1,2,2,8,0,0,0,0]};

figure, hold on
for ni = 1 : length(nonlearner)
    tempVec = [];
    for i = 1 : length(nonlearner{ni})
        tempVec = [tempVec, repmat(i, [1,nonlearner{ni}(i)])];
    end
    plot(1:sum(nonlearner{ni}), tempVec+0.05*ni, 'color', [0.7 0.7 0.7], 'linewidth', 1)
end
for ni = 1 : length(learner)
    tempVec = [];
    for i = 1 : length(learner{ni})
        tempVec = [tempVec, repmat(i, [1,learner{ni}(i)])];
    end
    plot(1:sum(learner{ni}), tempVec+0.05*ni, 'k', 'linewidth', 1)
end
xlabel('Sessions')
ylim([0 9])
yticks(1.15:1:8.15)
set(gca, 'fontsize', 12, 'fontname', 'Arial')




%%
%% Fig S1C
%%

%% Whisking rate
whiskerBaseDir = [baseDir, 'Whisker\'];
mice = [25,27,30,36,39,52];
sessions = {[4,19], [3,10], [3,21], [1,17], [1,23], [3,21]};

maxTimepoint = 1900; % hardcoded. 4 s will be about 1240 points. 1900 ~ 6.12 s
thetaAll = cell(length(mice),2);
amplitudeAll = cell(length(mice),2);
midpointAll = cell(length(mice),2);
whiskRateAll = cell(length(mice),2);
for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        d = sprintf('%sJK%03dS%02d',whiskerBaseDir,mouse,session);
        wla = Whisker.WhiskerTrialLite_2padArray(d);
        catchTrialInd = find(cellfun(@(x) strcmp(x.trialType,'oo'), wla.trials));
        poleTrialInd = setdiff(1:length(wla.trials),catchTrialInd);
        theta = cellfun(@(x) x.thetaAtBase{1}',wla.trials(poleTrialInd)','un',0);
        
        maxLength = max(cellfun(@length, theta));
        whiskRate = zeros(length(theta),1);
        thetaMat = nan(length(theta),maxLength);
        amplitudeMat = nan(length(theta),maxLength);
        midpointMat = nan(length(theta),maxLength);
        for ti = 1 : length(theta)
            trialInd = poleTrialInd(ti);
            [onsetFrames, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(theta{ti});
            tempLength = length(theta{ti});
            thetaMat(ti,1:tempLength) = theta{ti};
            amplitudeMat(ti,1:tempLength) = amplitude;
            midpointMat(ti,1:tempLength) = midpoint;
            frameRate = 1 / wla.trials{trialInd}.framePeriodInSec;
            poleUpFrame = wla.trials{trialInd}.poleUpFrames(1);
            samplingPeriodFrames = floor(poleUpFrame-0.1*frameRate):floor(poleUpFrame+frameRate);
            numOnset = sum(ismember(onsetFrames,samplingPeriodFrames));
            whiskRate(ti) = numOnset / 1.1;
        end        
        amplitudeAll{mi,si} = nanmean(amplitudeMat(:,1:maxTimepoint));
        midpointAll{mi,si} = nanmean(midpointMat(:,1:maxTimepoint));
        whiskRateAll{mi,si} = whiskRate;
        thetaAll{mi,si} = nanmean(thetaMat(:,1:maxTimepoint));
    end
end

whiskRateMat = zeros([size(whiskRateAll)]);
for i = 1 : 6
    for j = 1 : 2
        whiskRateMat(i,j) = mean(whiskRateAll{i,j});
    end
end
figure('units','normalized','position',[0.2 0.2 0.2 0.4]), hold on
for i = 1 : size(whiskRateMat,1)
    plot([1,2],whiskRateMat(i,:), 'ko-')
end
for i = 1 : 2
    errorbar(i,mean(whiskRateMat(:,i)), sem(whiskRateMat(:,i)), 'ro')
end
xlim([0.5 2.5]), ylim([0 12])
xticks([1,2])
xticklabels({'Naive', 'Expert'})
ylabel('Whisk rate (Hz)')
% title('Mean whisk rate during sampling period')
set(gca,'fontsize',12, 'fontname', 'Arial')

[~,p] = ttest(whiskRateMat(:,1)-whiskRateMat(:,2))
% p = 0.4124

%% Touch duration and counts
% regardless of object angle
calciumBaseDir = [baseDir, 'Calcium\'];
mice = [25,27,30,36,39,52];
sessions = {[4,19], [3,10], [3,21], [1,17], [1,23], [3,21]};

% features: (1) maxDtheta, (2) maxDphi, (3) maxDkappaH, (4) maxDkappaV, (5) slide distance,
% (6) touch duration, (7) theta, (8) phi, (9) kappaH, (10) kappaV, (11) arc length, (12) touch counts

tDuration = zeros(length(mice),2); % (:,1) naive, (:,2) expert
tCounts = zeros(length(mice),2);
for mi = 1 : length(mice)
    mouse = mice(mi);
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
        load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumBaseDir,mouse,mouse,session),'u')
        % select trials having touch before first lick
        noPoleTrialInds = find(cellfun(@(x) strcmp(x.trialType,'oo'), u.trials));
        poleTrialInds = setdiff(1:length(u.trials), noPoleTrialInds);
        answerTrialInds = find(cellfun(@(x) ~isempty(x.answerLickTime), u.trials));
        testInds = intersect(poleTrialInds, answerTrialInds);
        
        poleUpTime = cellfun(@(x) x.poleUpTime(1), u.trials(testInds), 'un', 0);
        allLickTimes = cellfun(@(x) union(x.leftLickTime, x.rightLickTime), u.trials(testInds), 'un', 0);
        firstLickTimes = cellfun(@(x,y) x(find(x > y, 1, 'first')), allLickTimes, poleUpTime, 'un', 0);
        allTouchTimes = cellfun(@(x) x.whiskerTime(cellfun(@(y) y(1), x.protractionTouchChunksByWhisking)), u.trials(testInds), 'un', 0);
        touchBeforeFirstLickInds = cellfun(@(x,y) find(x<y), allTouchTimes, firstLickTimes, 'un', 0); % there can be empty cells (no touch before first lick)
        
        tDuration(mi,si) = nanmean(cellfun(@(x,y) nanmean(x.protractionTouchDurationByWhisking(y)), u.trials(testInds), touchBeforeFirstLickInds));
        
        % Touch count can have 0
        tCounts(mi,si) = mean(cellfun(@length, touchBeforeFirstLickInds));
    end
end

%%
figure('units','normalized','position',[0.2 0.2 0.2 0.4]), hold on
for i = 1 : size(tDuration,1)
    plot([1,2],tDuration(i,:)*1000, 'ko-')
end
for i = 1 : 2
    errorbar(i,mean(tDuration(:,i))*1000, sem(tDuration(:,i))*1000, 'ro')
end
xlim([0.5 2.5])
xticks([1,2])
xticklabels({'Naive', 'Expert'})
ylim([0 27])
ylabel('Mean touch duration (ms)')
% title('Mean whisk rate during sampling period')
set(gca,'fontsize',12, 'fontname', 'Arial')

[~,p] = ttest(tDuration(:,1)-tDuration(:,2))
% 0.2532

%%
figure('units','normalized','position',[0.2 0.2 0.2 0.4]), hold on
for i = 1 : size(tCounts,1)
    plot([1,2],tCounts(i,:), 'ko-')
end
for i = 1 : 2
    errorbar(i,mean(tCounts(:,i)), sem(tCounts(:,i)), 'ro')
end
xlim([0.5 2.5])
xticks([1,2])
xticklabels({'Naive', 'Expert'})
ylim([0 10])
ylabel('Mean touch counts')
% title('Mean whisk rate during sampling period')
set(gca,'fontsize',12, 'fontname', 'Arial')

[~,p] = ttest(tCounts(:,1)-tCounts(:,2))
% 0.6489

%%
%% Fig S1D-F. Operant aspect of learning - lick pattern stereotypy
%%

trialNums = cell(length(mice),1);
correctLicks = cell(length(mice),1); % lick timing on the correct lick port (regardless of the trial type)
wrongLicks = cell(length(mice),1); % lick timing on the wrong lick port

for mi = 1 : length(mice)
    mouseName = sprintf('JK%03d',mice(mi));
    load([baseDir, '\Behavior\', 'behavior_', mouseName])
    
    trialNums{mi} = cell(length(sessions{mi}),1);
    correctLicks{mi} = cell(length(sessions{mi}),1);
    wrongLicks{mi} = cell(length(sessions{mi}),1);
    
    for si = 1 : length(sessions{mi})
        sessionName = sprintf('S%02d',sessions{mi}(si));
        bind = find(cellfun(@(x) strcmp(x.sessionName, sessionName), b));
        currB = b{bind};
        catchInd = find(cellfun(@(x) strcmp(x.trialType,'oo'), currB.trials));
        missInd = find(currB.trialCorrects == -1);
        behaviorTns = currB.trialNums(setdiff(1:length(currB.trials),union(missInd, catchInd)));
        trialNums{mi}{si} = behaviorTns;

        btninds = find(ismember(currB.trialNums, trialNums{mi}{si}));
       
        correctLicks{mi}{si} = cell(length(trialNums{mi}{si}),1);
        wrongLicks{mi}{si} = cell(length(trialNums{mi}{si}),1);
        
        for ti = 1 : length(trialNums{mi}{si})
            if currB.trials{btninds(ti)}.trialType(1) == 'r'
                wrongLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesLeft;
                correctLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesRight;
            elseif currB.trials{btninds(ti)}.trialType(1) == 'l'
                wrongLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesRight;
                correctLicks{mi}{si}{ti} = currB.trials{btninds(ti)}.beamBreakTimesLeft;
            else
                error('trialType is wrong')
            end
        end
    end
end

%%
%% Fig S1D. Averaged correct lick rates
%%

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
learner = [1,2,3,4,7,9];
nonlearner = [5,6,8,10,11,12];
sessions = {[1:17],[1:7],[1:20],[1:16],[1:24],[1:31],[1:18],[1:30],[1:20],[1:21],[1:24],[1:13]};
step = 0.1; % lick histogram step, in s
maxTime = 8; % in s
bins = 0:step:maxTime;
naiveWrongLick = zeros(length(bins)-1,length(learner));
naiveCorrectLick = zeros(length(bins)-1,length(learner));
expertWrongLick = zeros(length(bins)-1,length(learner));
expertCorrectLick = zeros(length(bins)-1,length(learner));
firstWrongLick = zeros(length(bins)-1,length(nonlearner));
firstCorrectLick = zeros(length(bins)-1,length(nonlearner));
lastWrongLick = zeros(length(bins)-1,length(nonlearner));
lastCorrectLick = zeros(length(bins)-1,length(nonlearner));
for mi = 1 : length(learner)
    for si = 1 : 3
        temp = wrongLicks{learner(mi)}{si};
        tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        naiveWrongLick(:,mi) = naiveWrongLick(:,mi) + tempwl/3;
        temp = correctLicks{learner(mi)}{si};
        tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        naiveCorrectLick(:,mi) = naiveCorrectLick(:,mi) + tempcl/3;
    end
    
    for si = length(sessions{learner(mi)}) - 4 : length(sessions{learner(mi)}) - 2
        temp = wrongLicks{learner(mi)}{si};
        tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        expertWrongLick(:,mi) = expertWrongLick(:,mi) + tempwl/3;
        temp = correctLicks{learner(mi)}{si};
        tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        expertCorrectLick(:,mi) = expertCorrectLick(:,mi) + tempcl/3;
    end
end

for mi = 1 : length(nonlearner)
    for si = 1 : 3
        temp = wrongLicks{nonlearner(mi)}{si};
        tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        firstWrongLick(:,mi) = firstWrongLick(:,mi) + tempwl/3;
        temp = correctLicks{nonlearner(mi)}{si};
        tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        firstCorrectLick(:,mi) = firstCorrectLick(:,mi) + tempcl/3;
    end
    
    for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
        temp = wrongLicks{nonlearner(mi)}{si};
        tempwl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        lastWrongLick(:,mi) = lastWrongLick(:,mi) + tempwl/3;
        temp = correctLicks{nonlearner(mi)}{si};
        tempcl = smooth(histcounts(cell2mat(temp), bins)) / step / length(temp);
        lastCorrectLick(:,mi) = lastCorrectLick(:,mi) + tempcl/3;
    end
end

c = parula(10);
c1 = c(4,:);
c2 = c(2,:);
w1 = c(10,:);
w2 = c(8,:);
step = 0.1; % lick histogram step, in s
maxTime = 8; % in s
bins = 0:step:maxTime;
figure, subplot(211), 
plot(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
hold on,
plot(bins(1:end-1)+step/2, mean(expertCorrectLick,2), 'color', c1, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(naiveWrongLick,2), '--', 'color', w2, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(expertWrongLick,2), 'color', w2, 'linewidth', 3);

boundedline(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), std(naiveCorrectLick, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
boundedline(bins(1:end-1)+step/2, mean(expertCorrectLick,2), std(expertCorrectLick, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
boundedline(bins(1:end-1)+step/2, mean(naiveWrongLick,2), std(naiveWrongLick, [], 2), '--', 'cmap', w2, 'alpha', 'transparency', 0.2);
boundedline(bins(1:end-1)+step/2, mean(expertWrongLick,2), std(expertWrongLick, [], 2), 'cmap', w2, 'alpha', 'transparency', 0.2);

plot(bins(1:end-1)+step/2, mean(naiveCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(expertCorrectLick,2), 'color', c1, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(naiveWrongLick,2), '--', 'color', w2, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(expertWrongLick,2), 'color', w2, 'linewidth', 3);
plot([1 1], [0 6], 'k--')
plot([2 2], [0 6], 'k--')
set(gca, 'fontsize', 12, 'fontname', 'Arial', 'box', 'off')
title('Learners'), ylabel('Lick rate (Hz)'), ylim([0 6]), xlim([0 6]), xticklabels({})
legend({'Correct (first 3 sessions)'; 'Correct (last 3 sessions)'; 'Wrong (first 3 sessions)'; 'Wrong (last 3 sessions)'}, 'box', 'off')

subplot(212), boundedline(bins(1:end-1)+step/2, mean(firstCorrectLick,2), std(firstCorrectLick, [], 2), '--', 'cmap', c1, 'alpha', 'transparency', 0.2);
hold on,
boundedline(bins(1:end-1)+step/2, mean(lastCorrectLick,2), std(lastCorrectLick, [], 2), 'cmap', c1, 'alpha', 'transparency', 0.2);
boundedline(bins(1:end-1)+step/2, mean(firstWrongLick,2), std(firstWrongLick, [], 2), '--', 'cmap', w2, 'alpha', 'transparency', 0.2);
boundedline(bins(1:end-1)+step/2, mean(lastWrongLick,2), std(lastWrongLick, [], 2), 'cmap', w2, 'alpha', 'transparency', 0.2);

plot(bins(1:end-1)+step/2, mean(firstCorrectLick,2), '--', 'color', c1, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(lastCorrectLick,2), 'color', c1, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(firstWrongLick,2), '--','color', w2, 'linewidth', 3);
plot(bins(1:end-1)+step/2, mean(lastWrongLick,2), 'color', w2, 'linewidth', 3);
plot([1 1], [0 6], 'k--')
plot([2 2], [0 6], 'k--')
set(gca, 'fontsize', 12, 'fontname', 'Arial', 'box', 'off')
title('Non-learners'), ylabel('Lick rate (Hz)'), ylim([0 6]), xlim([0 6]), xlabel('Time (s)')

%%
%% Fig S1E. Correlation of averaged lick patterns from each session to the template
%%

step = 0.1; % lick histogram step, in s
maxTime = 8; % in s
bins = 0:step:maxTime;
learnerWlicktemplate = zeros(length(learner),length(bins)-1); % wrong lick templates for learners (lick rate on the wrong side of the lick port)
nonlearnerWlicktemplate = zeros(length(nonlearner),length(bins)-1);% wrong lick templates for non-learners
learnerClicktemplate = zeros(length(learner),length(bins)-1); % correct lick templates for learners (lick rate on the correct side of the lick port)
nonlearnerClicktemplate = zeros(length(nonlearner),length(bins)-1); % correct lick templates for non-learners
for mi = 1 : length(learner)
    for si = length(sessions{learner(mi)}) - 2 : length(sessions{learner(mi)})
        tempwl = wrongLicks{learner(mi)}{si};
        tempcl = correctLicks{learner(mi)}{si};

        temphc = histcounts(cell2mat(tempwl), bins) / step / length(tempwl);
        learnerWlicktemplate(mi,:) = learnerWlicktemplate(mi,:) + temphc/3;
        temphc = histcounts(cell2mat(tempcl), bins) / step / length(tempcl);
        learnerClicktemplate(mi,:) = learnerClicktemplate(mi,:) + temphc/3;
    end
end

for mi = 1 : length(nonlearner)
    for si = length(sessions{nonlearner(mi)}) - 2 : length(sessions{nonlearner(mi)})
        tempwl = wrongLicks{nonlearner(mi)}{si};
        tempcl = correctLicks{nonlearner(mi)}{si};

        temphc = histcounts(cell2mat(tempwl), bins) / step / length(tempwl);
        nonlearnerWlicktemplate(mi,:) = nonlearnerWlicktemplate(mi,:) + temphc/3;
        temphc = histcounts(cell2mat(tempcl), bins) / step / length(tempcl);
        nonlearnerClicktemplate(mi,:) = nonlearnerClicktemplate(mi,:) + temphc/3;
    end
end

learnerWLickCorr = cell(length(learner),1);
for mi = 1 : length(learner)
    learnerWLickCorr{mi} = zeros(length(sessions{learner(mi)})-2,1);
    for si = 1 : length(sessions{learner(mi)})-2
        tempwl = wrongLicks{learner(mi)}{si};
        tempcl = correctLicks{learner(mi)}{si};
        temphc = histcounts(cell2mat(tempwl), bins) / step / length(tempwl);
        learnerWLickCorr{mi}(si) = corr(temphc', learnerWlicktemplate(mi,:)');
    end
end

nonlearnerWLickCorr = cell(length(nonlearner),1);
for mi = 1 : length(nonlearner)
    nonlearnerWLickCorr{mi} = zeros(length(sessions{nonlearner(mi)})-2,1);
    for si = 1 : length(sessions{nonlearner(mi)})
        tempwl = wrongLicks{nonlearner(mi)}{si};
        tempcl = correctLicks{nonlearner(mi)}{si};
        temphc = histcounts(cell2mat(tempwl), bins) / step / length(tempwl);
        nonlearnerWLickCorr{mi}(si) = corr(temphc', nonlearnerWlicktemplate(mi,:)');
    end
end

learnerCLickCorr = cell(length(learner),1);
for mi = 1 : length(learner)
    learnerCLickCorr{mi} = zeros(length(sessions{learner(mi)})-2,1);
    for si = 1 : length(sessions{learner(mi)})-2
        tempwl = wrongLicks{learner(mi)}{si};
        tempcl = correctLicks{learner(mi)}{si};
        temphc = histcounts(cell2mat(tempcl), bins) / step / length(tempcl);
        learnerCLickCorr{mi}(si) = corr(temphc(1:20)', learnerClicktemplate(mi,1:20)');
    end
end

nonlearnerCLickCorr = cell(length(nonlearner),1);
for mi = 1 : length(nonlearner)
    nonlearnerCLickCorr{mi} = zeros(length(sessions{nonlearner(mi)})-2,1);
    for si = 1 : length(sessions{nonlearner(mi)})
        tempwl = wrongLicks{nonlearner(mi)}{si};
        tempcl = correctLicks{nonlearner(mi)}{si};
        temphc = histcounts(cell2mat(tempcl), bins) / step / length(tempcl);
        nonlearnerCLickCorr{mi}(si) = corr(temphc(1:20)', nonlearnerClicktemplate(mi,1:20)');
    end
end

% wrong lick
figure,
hold on

for i = 1 : length(learner)

    plot(1:length(learnerWLickCorr{i}), learnerWLickCorr{i}(1:end), 'k-', 'linewidth', 3)
end
for i = 1 : length(nonlearner)
    plot(1:length(nonlearnerWLickCorr{i}), nonlearnerWLickCorr{i}(1:end), 'color', [0.7 0.7 0.7], 'linewidth', 3)
end

title('Wrong lick'), ylabel('Correlation'), xlabel('Session')
set(gca, 'fontsize', 14, 'box', 'off')

% correct lick frequency, before drinking

figure, hold on

for i = 1 : length(learner)

    plot(1:length(learnerCLickCorr{i}), learnerCLickCorr{i}(1:end), 'k-', 'linewidth', 3)
end
for i = 1 : length(nonlearner)
    plot(1:length(nonlearnerCLickCorr{i}), nonlearnerCLickCorr{i}(1:end), 'color', [0.7 0.7 0.7], 'linewidth', 3)
end

title('Correct lick before drinking'), ylabel('Correlation'), xlabel('Session')
set(gca, 'fontsize', 14, 'box', 'off')

%% 
%% Fig S1F. Lick pattern stereotypy, only from learners
%%

figure, hold on
totalWlickCorr = nan(6,max(cellfun(@length, learnerWLickCorr)));
totalClickCorr = nan(6,max(cellfun(@length, learnerWLickCorr)));
for i = 1 : 6
    temp = learnerWLickCorr{i};
    totalWlickCorr(i,end-length(temp)+1:end) = temp;
    temp = learnerCLickCorr{i};
    totalClickCorr(i,end-length(temp)+1:end) = temp;
end

plot(1:size(totalClickCorr,2), nanmean(totalClickCorr), 'color', c1, 'linewidth', 3)
plot(1:size(totalWlickCorr,2), nanmean(totalWlickCorr), 'color', w2, 'linewidth', 3)
boundedline(1:size(totalClickCorr,2), nanmean(totalClickCorr), nanstd(totalClickCorr)/sqrt(6),'cmap', c1, 'alpha', 'transparency', 0.2)
boundedline(1:size(totalWlickCorr,2), nanmean(totalWlickCorr), nanstd(totalWlickCorr)/sqrt(6),'cmap', w2, 'alpha', 'transparency', 0.2)
plot(1:size(totalClickCorr,2), nanmean(totalClickCorr), 'color', c1, 'linewidth', 3)
plot(1:size(totalWlickCorr,2), nanmean(totalWlickCorr), 'color', w2, 'linewidth', 3)

xticks(1:size(totalClickCorr,2))
xticklabels(-size(totalClickCorr,2)+2+(1:size(totalClickCorr,2)))
title('Lick pattern stereotypy'), ylabel('Correlation'), xlabel('Session')
set(gca, 'fontsize', 14, 'box', 'off'), legend({'Correct lick', 'Wrong lick'}, 'box', 'off', 'location', 'southeast')




