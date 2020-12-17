
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% basic settings
calciumDir = [baseDir, 'Calcium\'];
mouse = 39;
session = 1;
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

naiveInd = find(mice == mouse);

repeat = 10;
% %% dependent settings
ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_lasso_NC_R', mouse, session);
angletuningFn = sprintf('JK%03dS%02dangle_tuning_lasso_predecision_NC_permTestCorrected', mouse, session);
% %% load files
load([calciumDir, 'glmResults_devExp_touch_NC'], 'naive', 'expert')
% cd(sprintf('%s%03d',baseDir, mouse))
load(sprintf('%s%03d\\%s',calciumDir, mouse, ufn), 'u') % loading u
load(sprintf('%s%03d\\%s01',calciumDir, mouse, glmfnBase), 'allPredictors', 'posShift');
load(sprintf('%s%03d\\%s',calciumDir, mouse, angletuningFn))

touch = load([calciumDir, 'glmResults_devExp_touch_NC']);
wkv = load([calciumDir, 'glmResults_devExp_WKV_touchCell_NC']);

%%
%% Fig 5A
%%
%% Figure of Whisker GLM
cID = 5087;
wkvglmfnBase = sprintf('%s%03d\\glmWhisker_lasso_touchCell_NC_JK%03dS%02d_R', calciumDir, mouse, mouse, session);
load([wkvglmfnBase, '01'], 'allPredictors', 'posShift', 'cIDAll', 'fitCoeffs', 'indPartial')
ci = find(cIDAll == cID);

wkv.naive(7).allDE(ci)

testCoeff = zeros(size(fitCoeffs{ci},1),10);
for i = 1 : 10
    load([wkvglmfnBase,sprintf('%02d',i)], 'fitCoeffs')
    testCoeff(:,i) = fitCoeffs{ci};
end

traces = get_traces_per_cell(u, cID, allPredictors, testCoeff, posShift);
figure('units', 'normalized', 'outerposition', [0.3 0.3 0.4 0.4]), hold on
plot(min_max_normalization(traces.calcium)*4 + 8 , 'color', [0.1 0.8 0.1], 'linewidth', 1)
plot(min_max_normalization(traces.spikes)*2.5 + 6.2, 'color', [0.8 0.1 0.1], 'linewidth', 1)

plot(min_max_normalization(traces.model)*2.5 + 6.2, 'k-', 'linewidth', 1)
plot(min_max_normalization(traces.predictors(:,10)) + 4.4, 'color', 'g', 'linewidth', 1) % dKv
plot(min_max_normalization(traces.predictors(:,4)) + 3.4, 'color', 'r', 'linewidth', 1) % dPhi
plot(min_max_normalization(traces.predictors(:,13)) + 3.6, 'color', 'b', 'linewidth', 1) % slide distance

plot(min_max_normalization(traces.predictors(:,51)) + 2.5, 'color', [0.7 0.7 0.7], 'linewidth', 1)
plot(min_max_normalization(traces.predictors(:,58)) + 1.8, 'color', [0.7 0.7 0.7], 'linewidth', 1)
plot(min_max_normalization(traces.predictors(:,65)) + 1.1, 'color', [0.7 0.7 0.7], 'linewidth', 1)
lickL = min_max_normalization(traces.predictors(:,72));
lickR = min_max_normalization(traces.predictors(:,76));
plot(lickL, 'm', 'linewidth', 1)
plot(lickR, 'c', 'linewidth', 1)
rewardL = min_max_normalization(traces.predictors(:,41));
rewardL(rewardL<max(rewardL)) = 0;
rewardL(isnan(rewardL)) = 0;
rewardR = min_max_normalization(traces.predictors(:,46));
rewardR(rewardR<max(rewardR)) = 0;
rewardR(isnan(rewardR)) = 0;
reward = rewardL + rewardR;
rewarded = find(reward);
rdiffOne = find(diff(find(reward))==1);
reward(rewarded(rdiffOne+1)) = nan;
reward(reward==0) = nan;
sound = min_max_normalization(traces.predictors(:,37));
sound(sound<max(sound)) = nan;
plot(reward-1.5, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'b')
plot(sound-1.5, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'k')
onedF = 4 / max(traces.calcium);

plot([930 930], [10 10+onedF], 'color', [0.1 0.8 0.1], 'linewidth', 1) % 1 dF/F0 scale bar
plot([900 900 + u.frameRate*5], [10 10], 'k-', 'linewidth', 1) % 5 s scale bar
xlim([240 940]), ylim([-0.5 12])

axis off
set(gcf, 'InvertHardCopy', 'off', 'color', 'w');



%%
%%
%%




%%
%% S4A
%%
%% object model GOF - sensorimotor model GOF

touch = load([calciumDir, 'glmResults_devExp_touch_NC']);
wkv = load([calciumDir, 'glmResults_devExp_WKV_touchCell_NC']);

expertInd = [1:4,7,9];

%% all naive

diffDE = cell(length(wkv.naive),1);
for i = 1 : length(wkv.naive)
    tempCID = wkv.naive(i).cellID;
    pfTemp = zeros(length(tempCID),1); % parfor temp
    for ci = 1 : length(tempCID)    
        wkvInd = find(wkv.naive(i).cellID == tempCID(ci));
        touchInd = find(touch.naive(i).cellID == tempCID(ci)); % this way is safer than assuming sorted
        pfTemp(ci) = touch.naive(i).allDE(touchInd) - wkv.naive(i).allDE(wkvInd);
    end
    diffDE{i} = pfTemp;
end

histRange = [-0.15:0.01:0.25];
pdfDiff = zeros(length(diffDE), length(histRange)-1);
for i = 1 : length(diffDE)
    pdfDiff(i,:) = histcounts(diffDE{i}, histRange, 'normalization', 'probability');
end
figure, hold on
boundedline(histRange(1:end-1), mean(pdfDiff), std(pdfDiff), 'cmap', [0 0 0])
graphMaxVal = max(mean(pdfDiff) + std(pdfDiff));

graphMean = mean(cell2mat(diffDE))
graphStd = std(cell2mat(diffDE))

ylim([0 0.2])
yval = ylim();
errorbar(graphMean, graphMaxVal + (yval(2)-graphMaxVal)/2, graphStd, 'ro-','horizontal')
plot([0.1 0.1], yval, '--', 'color', [0.7 0.7 0.7])
plot([-0.1 -0.1], yval, '--', 'color', [0.7 0.7 0.7])

xlabel('Difference between models')
ylabel('Proportion')
set(gca,'fontsize',12,'fontname','Arial')

%%
%% S4B
%%
%% compare angle-tuning correlation in angle-tuned neurons
% between full whisker model and whisker-only model
%% naive all

% 
% featureNames = {'maxDq', 'maxDf', 'maxDkH', 'maxDkV', 'max(slide distance)', 'max(duration)', ...    
%                             'q', 'f', 'kH', 'kV', 'arc length', 'touch count'};

loadFn = 'wkv_angle_tuning_model_v9';
load(sprintf('%s%s',calciumDir, loadFn)); 

angles = 45:15:135;
learnerInd = [1:4,7,9];

% correct for JK027 when comparing between naive and expert (n = 6)
naiveLearner = naive(learnerInd);
jk027inds = find(naiveLearner(2).cIDAll > 5000);
jk027AllLength = length(naiveLearner(2).cIDAll);

fn = fieldnames(naiveLearner);
for i = 1 : length(fn)
    if length(naiveLearner(2).(fn{i})) == jk027AllLength
        naiveLearner(2).(fn{i}) = naiveLearner(2).(fn{i})(jk027inds);
    end
end

atCorr = cell(length(naive),2); %(:,1) full model, (:,2) whisker-only model

histRange = 0:0.02:1;
distCorrFull = zeros(length(naive),length(histRange)-1);
distCorrWhisker = zeros(length(naive),length(histRange)-1);

for mi = 1 : length(naive)
    ind = find(naive(mi).tuned);
    fullCorrMat = cell2mat(naive(mi).atCorrFull(ind));
    atCorr{mi,1} = fullCorrMat(:,1);
    whiskerCorrMat = cell2mat(naive(mi).atCorrWhisker(ind));
    atCorr{mi,2} = whiskerCorrMat(:,1);
    
    distCorrFull(mi,:) = histcounts(atCorr{mi,1}, histRange, 'norm', 'prob');
    distCorrWhisker(mi,:) = histcounts(atCorr{mi,2}, histRange, 'norm', 'prob');
end
%
figure, hold on
plot(histRange(2:end), mean(distCorrFull), 'k')
plot(histRange(2:end), mean(distCorrWhisker), 'c')
legend({'Full model', 'Whisker-only model'}, 'autoupdate', false, 'location', 'northwest')
boundedline(histRange(2:end), mean(distCorrFull), sem(distCorrFull), 'k')
boundedline(histRange(2:end), mean(distCorrWhisker), sem(distCorrFull), 'c')
xlabel('Correlation')
ylabel('Proportion')
title({'Naive all (n = 12)'; sprintf('Full model mean %.2f \\pm %.2f', ...
    mean(cellfun(@mean, atCorr(:,1))), sem(cellfun(@mean, atCorr(:,1)))) ;
    sprintf('Whisker-only model mean %.2f \\pm %.2f', ...
    mean(cellfun(@mean, atCorr(:,2))), sem(cellfun(@mean, atCorr(:,2))))})







%%
%% 4B - example plots for impact on angle tuning
%%

load([calciumDir, 'wkv_angle_tuning_model_v9'], 'naive');
mouse = naive(7);

%%
angles = 45:15:135;
tempInd = find(mouse.cIDAll == 5104); % 5104
% cID = cellIDs(tempInd);

spikesTemp = (cellfun(@mean, mouse.spikeAngleAll{tempInd,1}));
whiskerTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(1,:)));
woDkVTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(5,:)));
woSDTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(6,:)));

spikes = (spikesTemp-mean(spikesTemp)) / std(spikesTemp);
whiskerModel = (whiskerTemp-mean(whiskerTemp)) / std(whiskerTemp);
woDkV = (woDkVTemp-mean(woDkVTemp)) / std(woDkVTemp);
woSD = (woSDTemp-mean(woSDTemp)) / std(woSDTemp);

figure, hold on
plot(angles, spikes, 'ro-', 'MarkerFaceColor','r', 'linewidth', 1);
ylabel('Response (standardized)')
plot(angles, whiskerModel, 'o-', 'color',[0.6 0.6 0.6], 'linewidth', 3);
plot(angles, woDkV, 'go--', 'linewidth', 1)
plot(angles, woSD, 'bo--', 'linewidth', 1)
xticks(angles)
xlabel('Object angle (\circ)')
legend({'Inferred spikes', 'Full model', '-maxDkV', '-max(slide distance)'}, 'location', 'northwest', 'box', 'off')


%%
%% 4B (top right)
%%
%%
corrVal = zeros(1,3);
corrVal(1) = mouse.atCorrWhisker{tempInd}(1);
corrVal(2) = mouse.atCorrWhisker{tempInd}(5);
corrVal(3) = mouse.atCorrWhisker{tempInd}(6);
figure, hold on
bar(1, corrVal(1), 'k', 'baseValue', -1)
bar(2, corrVal(2), 'g', 'baseValue', -1)
bar(3, corrVal(3), 'b', 'baseValue', -1)
yticks([-1:0.5:1])
ylabel('Correlation')
xticks([1:3])
xticklabels({'Full model', '-vertical distance', '-slide distance'})
xtickangle(45)


%% confirmation
corrVal = zeros(1,3);
corrVal(1) = corr(spikesTemp', whiskerTemp');
corrVal(2) = corr(spikesTemp', woDkVTemp');
corrVal(3) = corr(spikesTemp', woSDTemp');
figure, hold on
bar(1, corrVal(1), 'k', 'baseValue', -1)
bar(2, corrVal(2), 'g', 'baseValue', -1)
bar(3, corrVal(3), 'b', 'baseValue', -1)
ylabel('Correlation')
xticks([1:3])
xticklabels({'Full model', '-maxDkV', '-max(slide distance)'})
xtickangle(45)






%%
%% 4B (bottome left)
%%
%% listening to slide distance, but not to dKv

tempInd = find(mouse.cIDAll == 5326);
spikesTemp = (cellfun(@mean, mouse.spikeAngleAll{tempInd,1}));
whiskerTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(1,:)));
woDkVTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(5,:)));
woSDTemp = (cellfun(@mean, mouse.whiskerOnlyAngleAll{tempInd}(6,:)));

spikes = (spikesTemp-mean(spikesTemp)) / std(spikesTemp);
whiskerModel = (whiskerTemp-mean(whiskerTemp)) / std(whiskerTemp);
woDkV = (woDkVTemp-mean(woDkVTemp)) / std(woDkVTemp);
woSD = (woSDTemp-mean(woSDTemp)) / std(woSDTemp);

figure, hold on
plot(angles, spikes, 'ro-', 'MarkerFaceColor','r', 'linewidth', 1);
ylabel('Response (standardized)')
plot(angles, whiskerModel, 'o-', 'color',[0.6 0.6 0.6], 'linewidth', 3);
plot(angles, woDkV, 'go--', 'linewidth', 1)
plot(angles, woSD, 'bo--', 'linewidth', 1)
xticks(angles)
xlabel('Object angle (\circ)')
legend({'Inferred spikes', 'Full model', '-maxDkV', '-max(slide distance)'}, 'location', 'northwest', 'box', 'off')

%%
%% 4B (bottom right)
%%
%%
corrVal = zeros(1,3);
corrVal(1) = corr(spikesTemp', whiskerTemp');
corrVal(2) = corr(spikesTemp', woDkVTemp');
corrVal(3) = corr(spikesTemp', woSDTemp');
figure, hold on
bar(1, corrVal(1), 'k', 'baseValue', -1)
bar(2, corrVal(2), 'g', 'baseValue', -1)
bar(3, corrVal(3), 'b', 'baseValue', -1)
yticks([-1:0.5:1])
ylabel('Correlation')
xticks([1:3])
xticklabels({'Full model', '-vertical bending', '-slide distance'})
xtickangle(45)








%%
%% 4C - Impact on angle tuning by each features from whisker-only models
%%

impactWhiskerNaiveAll = cell(length(naive),1);

for mi = 1 : length(naive)
    tunedInd = find(naive(mi).tuned);
    tempMat = cell2mat(naive(mi).atCorrWhisker(tunedInd));
    impactWhiskerNaiveAll{mi} = tempMat(:,1) - tempMat(:,2:end);
end

featureXpos = [8,9,10,11,13,12,1,2,3,4,5,6];
figure, hold on
tempMatWhisker = cell2mat(cellfun(@nanmean, impactWhiskerNaiveAll, 'un', 0));
bar(featureXpos, mean(tempMatWhisker), 'facecolor', 'k')
errorbar(featureXpos, mean(tempMatWhisker), sem(tempMatWhisker), 'k', 'lines', 'no')
ylabel('Impact on tuning')
ylim([-0.02 0.2]), yticks([0:0.05:0.2])
title('Naive all (n = 12)')
xticks([1:6, 8:13])
xticklabels({'Azimuthal angle', 'Vertical angle', 'Horizontal curvature', 'Vertical curvature', 'Arc length', 'Touch count', ...
    'Azimuthal push angle', 'Vertical displacement', 'Horizontal bending', 'Vertical bending', 'Touch duration', 'Slide distance'})
xtickangle(45)



%% (Impact on angle tuning by each features from whisker-only models)
%% (expert mice)
impactWhiskerNaiveAll = cell(length(expert),1);

for mi = 1 : length(expert)
    tunedInd = find(expert(mi).tuned);
    tempMat = cell2mat(expert(mi).atCorrWhisker(tunedInd));
    impactWhiskerNaiveAll{mi} = tempMat(:,1) - tempMat(:,2:end);
end

featureXpos = [8,9,10,11,13,12,1,2,3,4,5,6];
figure, hold on
tempMatWhisker = cell2mat(cellfun(@nanmean, impactWhiskerNaiveAll, 'un', 0));
bar(featureXpos, mean(tempMatWhisker), 'facecolor', 'k')
errorbar(featureXpos, mean(tempMatWhisker), sem(tempMatWhisker), 'k', 'lines', 'no')
ylabel('Impact on tuning')
ylim([-0.02 0.2]), yticks([0:0.05:0.2])
title('Expert (n = 6)')
xticks([1:6, 8:13])
xticklabels({'Azimuthal angle', 'Vertical angle', 'Horizontal curvature', 'Vertical curvature', 'Arc length', 'Touch count', ...
    'Azimuthal push angle', 'Vertical displacement', 'Horizontal bending', 'Vertical bending', 'Touch duration', 'Slide distance'})
xtickangle(45)




%% Fig 4D
%% Impact of dynamics whisker features on angle tuning
%% across tuned angles (defined from inferred spikes)
atCorrThreshold = 0.1;

sortedFeatureImportance = cell(length(naive), 2); %(:,1) from full model, (:,2) from whisker-only model
sortedFeatureInd = cell(length(naive),2);
cellDepths = cell(length(naive),1);
cellTunedAngle = cell(length(naive),2); %(:,1) from spikes, (:,2) from whisker-only model

for mi = 1 : length(naive)
    tempInd = find(naive(mi).tuned);
    cellDepths{mi} = naive(mi).depth(tempInd);
    tempSpikeAll = cell2mat(naive(mi).atSpikeAll(tempInd));
    [~, tempAngleInd] = max(tempSpikeAll, [], 2);
    cellTunedAngle{mi,1} = 30 + tempAngleInd * 15;
    tempWhiskerAll = cell2mat(cellfun(@(x) x(1,:), naive(mi).atWhiskerAll(tempInd), 'un', 0));
    [~, tempAngleInd] = max(tempWhiskerAll, [], 2);
    cellTunedAngle{mi,2} = 30 + tempAngleInd * 15;
    
    tempMat = cell2mat(naive(mi).atCorrFull(tempInd));
    [sortedFeatureImportance{mi,1}, sortedFeatureInd{mi,1}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
    tempMat = cell2mat(naive(mi).atCorrWhisker(tempInd));
    [sortedFeatureImportance{mi,2}, sortedFeatureInd{mi,2}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
end


impactFeature = zeros(length(naive), 12, length(angles));
propFeature = zeros(length(naive), 12, length(angles));

for mi = 1 : length(naive)
    tempInd = find(naive(mi).tuned);
    tunedAngleInd = cellfun(@(x) find(x == max(x)), naive(mi).atSpikeAll(tempInd));
    for ai = 1 : length(angles)
        angleInd = find(tunedAngleInd == ai);
        
        % proportion  (Whisker-only)
        tempMat = sortedFeatureImportance{mi,2}(angleInd,:);
        tempMatInd = find(tempMat(:) > atCorrThreshold);
        featureIndexMat = sortedFeatureInd{mi,2}(angleInd,:);
        allIndices = featureIndexMat(tempMatInd);
        for fi = 1 : 12
            propFeature(mi,fi,ai) = length(find(allIndices == fi)) / length(angleInd);
        end
        
        % mean impact of each feature from multi feature neurons (Whisker-only)
        for fi = 1 : 12
            tempImpact = zeros(length(angleInd),1);
            for i = 1 : length(angleInd)
                tempImpactVec = sortedFeatureImportance{mi,2}(angleInd(i),:);
                tempIndVec = sortedFeatureInd{mi,2}(angleInd(i),:);
                tempImpact(i) = tempImpactVec(find(tempIndVec == fi));
            end
            impactFeature(mi,fi,ai) = mean(tempImpact);
        end
        
    end
end

figure('units','inch','pos',[2 1 3 6]),
for ai = 1 : length(angles)
    subplot(7,1,ai), hold on
    tempMat = squeeze(impactFeature(:,[1,2,3,4,6,5],ai));
    bar(nanmean(tempMat), 'facecolor', 'k')
    errorbar(nanmean(tempMat), sem(tempMat), 'k', 'lines', 'no')

    if ai == 4
        ylabel('Mean impact on tuning')
    end
    
    if ai < length(angles)
        xticks([])
    else
        xticks(1:6)
        xticklabels({'Push angle', 'Vertical displacement', 'Horizontal bending', ...
            'Vertical bending', 'Touch duration', 'Slide distance'})
        xtickangle(45)
    end
    ylim([0 0.35])
    xlim([0.5 6.5])
end
sgtitle('Naive all (n = 12)')

























































%%
%% Fig 5C (top left)
%%
%% Show which features affect the correlation the most.

loadFn1 = 'modelAngleTuning_NC';
loadFn2 = 'modelAngleTuning_NC_combinations';
data1 = load(sprintf('%s%s',calciumDir, loadFn1), 'naive', 'expert'); % 
data2 = load(sprintf('%s%s',calciumDir, loadFn2), 'naive', 'expert'); % 
numMiceNaive = length(data1.naive);
numMiceExpert = length(data1.expert);
% whiskerGLM = load('Y:\Whiskernas\JK\suite2p\glmResults_devExp_WKV_touchCell_NC.mat', 'naive');

%% listening to dKv, but not to slide distance





%%
%% Fig 5D
%%
%% map position
% load u first
%% Whisker feature-angle tuning impact map
baseDir = 'D:\TPM\JK\suite2p\';
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
mouse = 39;
session = 1;
plane = 5;
load(sprintf('%s%03d\\angle_tuning_model_lasso_NC_preAnswer_perTouch_JK%03dS%02d.mat', baseDir, mouse, mouse, session), 'spkValAllCell', 'tuneAngleAllCell')
ufn = sprintf('%s%03d\\UberJK%03dS%02d_NC', baseDir, mouse, mouse, session);
load(ufn, 'u') % loading u
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC', baseDir, mouse, mouse,session), 'spk', 'info')
%%
numFeature = 12;
corrFeaturesOutNaive = cell(1,length(spkValAllCell));

for ci = 1 : length(spkValAllCell)    
    corrFeaturesOutNaive{ci} = zeros(1,numFeature+1);
    if isfinite(tuneAngleAllCell(ci,1)) % only in tuned neurons, because not-tuned neurons have NaN angle
        tempVal = corr(cellfun(@mean, spkValAllCell{(ci),1}), cellfun(@mean, spkValAllCell{(ci),3}));
        if isnan(tempVal) || tempVal < 0
            tempVal = 0;
        end
        corrFeaturesOutNaive{ci}(1) = tempVal;
        for fi = 1 : numFeature            
            tempVal = corr(cellfun(@mean, spkValAllCell{(ci),1}), cellfun(@mean, spkValAllCell{(ci),3+fi}));
            if isnan(tempVal) || tempVal < 0
                tempVal = 0;
            end
            corrFeaturesOutNaive{ci}(fi+1) = corrFeaturesOutNaive{ci}(1) - tempVal;
        end
    end
end
maxCorrVal = cellfun(@(x) max(x(2:end)), corrFeaturesOutNaive, 'un', 0);
maxCorrInd = cellfun(@(x,y) find(ismember(x(2:end),y),1), corrFeaturesOutNaive, maxCorrVal);

colors = jet(7); % because we know there is only up to 7 different features (differ between examples)
mimg = (mat2gray(u.mimg{plane}));
h = figure;
imshow((mimg))
hold on,
plot([u.c2xpoints, u.c2xpoints(1)], [u.c2ypoints, u.c2ypoints(1)], 'w--', 'linewidth', 2)

cIDlist = spk.touchID(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
ciList = find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000);
tunedAngles = spk.tunedAngle(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
indList = find(ismember(u.cellNums, cIDlist));

cellmap = u.cellmap{plane};
maxCorrIndsMap = zeros(length(cIDlist),1);
for ci = 1 : length(cIDlist)
    cID = cIDlist(ci);
    [ypix,xpix] = ind2sub(size(cellmap),find(cellmap == cID));    
    k = boundary(ypix, xpix);
    if tunedAngles(ci)
        if maxCorrVal{ciList(ci)} < 0.1
            patch(xpix(k), ypix(k), 'k', 'edgecolor', 'none')
        else
            tempFeatureInd = maxCorrInd(ciList(ci)); % from 1 to 12
            patch(xpix(k), ypix(k), colors(tempFeatureInd,:), 'edgecolor', 'none')
            maxCorrIndsMap(ci) = tempFeatureInd;
        end
    end
end

% scale bar edge. scale bar length 100 um
scaleBarPix = 100 / u.pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(u.mimg{plane});
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
expectedTextLength = 120;
text(sizes(2)-expectedTextLength, margins, [num2str(-round(mean(mean(u.fovdepth{plane})))), ' \mum'], 'fontsize', 18, 'fontweight', 'bold', 'color', 'w')
set(gcf, 'InvertHardCopy', 'off', 'color', 'w');
set(gca, 'fontsize', 12, 'fontname', 'Arial')

imcontrast(h)

%%
% saveDir = 'C:\Users\jinho\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig5 Sensory input for angle tuning\';
% fn = sprintf('example_FOV(JK%03dS%02dp%d)_whiskerImpactMap.eps', mouse, session, plane);
% export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
% fix_eps_fonts([saveDir, fn])

% and then...

cmtemp = zeros(size(u.cellmap{5}));
cmtemp(find(u.cellmap{5} == 5104)) = 1;
cmtemp(find(u.cellmap{5} == 5326)) = 2;

figure, imagesc(cmtemp)


%%
%% Fig 5E
%%
%% maxDkV and max(slide distance) affects the angle tuning the most.

%% but first off, how many of whisker models are angle-tuned?
tunedWhiskerModel = zeros(numMiceNaive,1);
for i = 1 : numMiceNaive
    indTunedSpikes = find(data1.naive(i).tunedAllCell(:,1));
    tunedWhiskerModel(i) = sum(data1.naive(i).tunedAllCell(indTunedSpikes,3)) / length(indTunedSpikes);
end

mean(tunedWhiskerModel)
sem(tunedWhiskerModel)

%%
tunedWhiskerModel = zeros(numMiceExpert,1);
for i = 1 : numMiceExpert
    indTunedSpikes = find(data1.expert(i).tunedAllCell(:,1));
    tunedWhiskerModel(i) = sum(data1.expert(i).tunedAllCell(indTunedSpikes,3)) / length(indTunedSpikes);
end

mean(tunedWhiskerModel)
sem(tunedWhiskerModel)

%%
%% touch, full, 12 drop-outs, others only (16 total)
corrTouchNaive = cell(numMiceNaive,1);
corrWhiskerNaive = cell(numMiceNaive,1);
corrIndsNaive = cell(numMiceNaive,1);
for i = 1 : numMiceNaive
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    corrIndsNaive{i} = indWhisker;
    
    corrTouchNaive{i} = zeros(length(indTouch),1);
    corrWhiskerNaive{i} = zeros(length(indWhisker),1);
    for j = 1 : length(indTouch)
        tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indTouch(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indTouch(j),2}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrTouchNaive{i}(j) = tempVal;
    end
    for j = 1 : length(indWhisker)
        tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),3}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrWhiskerNaive{i}(j) = tempVal;
    end
end

corrFeaturesOutNaive = cell(numMiceNaive,12);
for i = 1 : numMiceNaive
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesOutNaive{i,fi} = zeros(length(indWhisker),1);        
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),3+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesOutNaive{i,fi}(j) = tempVal;
        end
    end
end

corrOtherNaive = cell(numMiceNaive,1);
for i = 1 : numMiceNaive
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for j = 1 : length(indWhisker)
        tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j), 1}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrOtherNaive{i}(j) = tempVal;
    end    
end


corrFeaturesCombOutNaive = cell(numMiceNaive,4);
corrFeaturesCombInNaive = cell(numMiceNaive,4);
for i = 1 : numMiceNaive
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 4
        corrFeaturesCombOutNaive{i,fi} = zeros(length(indWhisker),1);
        corrFeaturesCombInNaive{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),13+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesCombOutNaive{i,fi}(j) = tempVal;
            
            tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),21+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesCombInNaive{i,fi}(j) = tempVal;
            
        end
    end
end


corrFeaturesInNaive = cell(numMiceNaive,12);
for i = 1 : numMiceNaive
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesInNaive{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.naive(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.naive(i).spkValAllCell{indWhisker(j),1+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesInNaive{i,fi}(j) = tempVal;
        end
    end
end





%% touch, full, 12 drop-outs, others only (16 total) - Expert
corrTouchExpert = cell(numMiceExpert,1);
corrWhiskerExpert = cell(numMiceExpert,1);
corrIndsExpert = cell(numMiceExpert,1);
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,2));
    indTouch = indTuned(indTemp);
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    corrIndsExpert{i} = indWhisker;
    
    corrTouchExpert{i} = zeros(length(indTouch),1);
    corrWhiskerExpert{i} = zeros(length(indWhisker),1);
    for j = 1 : length(indTouch)
        tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indTouch(j),1}), cellfun(@mean, data1.expert(i).spkValAllCell{indTouch(j),2}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrTouchExpert{i}(j) = tempVal;
    end
    for j = 1 : length(indWhisker)
        tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),3}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrWhiskerExpert{i}(j) = tempVal;
    end
end

corrFeaturesOutExpert = cell(numMiceExpert,12);
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesOutExpert{i,fi} = zeros(length(indWhisker),1);        
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),3+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesOutExpert{i,fi}(j) = tempVal;
        end
    end
end

corrOtherExpert = cell(numMiceExpert,1);
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for j = 1 : length(indWhisker)
        tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.expert(i).spkValAllCell{indWhisker(j), 1}));
        if isnan(tempVal)
            tempVal = 0;
        end
        corrOtherExpert{i}(j) = tempVal;
    end    
end


corrFeaturesCombOutExpert = cell(numMiceExpert,4);
corrFeaturesCombInExpert = cell(numMiceExpert,4);
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 4
        corrFeaturesCombOutExpert{i,fi} = zeros(length(indWhisker),1);
        corrFeaturesCombInExpert{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.expert(i).spkValAllCell{indWhisker(j),13+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesCombOutExpert{i,fi}(j) = tempVal;
            
            tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.expert(i).spkValAllCell{indWhisker(j),21+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesCombInExpert{i,fi}(j) = tempVal;
            
        end
    end
end


corrFeaturesInExpert = cell(numMiceExpert,12);
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    for fi = 1 : 12
        corrFeaturesInExpert{i,fi} = zeros(length(indWhisker),1);
        for j = 1 : length(indWhisker)
            tempVal = corr(cellfun(@mean, data1.expert(i).spkValAllCell{indWhisker(j),1}), cellfun(@mean, data2.expert(i).spkValAllCell{indWhisker(j),1+fi}));
            if isnan(tempVal)
                tempVal = 0;
            end
            corrFeaturesInExpert{i,fi}(j) = tempVal;
        end
    end
end






%%
%% Fig S8G
%%
%% Distribution of whisker correlation
histRange = 0:0.1:1;
corrWhiskerDist = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'prob'), corrWhiskerNaive, 'un', 0));
figure, hold on
% boundedline(histRange(2:end)-0.05, mean(corrWhiskerDist), sem(corrWhiskerDist), 'k')
bar(histRange(2:end)-0.05, mean(corrWhiskerDist), 'k'), hold on
errorbar(histRange(2:end)-0.05, mean(corrWhiskerDist), sem(corrWhiskerDist), 'k', 'lines', 'no')
xlabel('Correlation with inferred spike')
ylabel('Proportion')
title('Naive')
set(gca, 'fontsize', 12, 'fontname','Arial', 'box', 'off');



%% Expert
histRange = 0:0.1:1;
corrWhiskerDist = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'prob'), corrWhiskerExpert, 'un', 0));
figure, hold on
% boundedline(histRange(2:end)-0.05, mean(corrWhiskerDist), sem(corrWhiskerDist), 'k')
bar(histRange(2:end)-0.05, mean(corrWhiskerDist), 'k'), hold on
errorbar(histRange(2:end)-0.05, mean(corrWhiskerDist), sem(corrWhiskerDist), 'k', 'lines', 'no')
xlabel('Correlation with inferred spike')
ylabel('Proportion')
title('Expert')
set(gca, 'fontsize', 12, 'fontname','Arial', 'box', 'off');


%%
%% Fig S8H
%%
% Impacts of each features
eachFeatureOut = cellfun(@mean, corrWhiskerNaive) - cellfun(@mean, corrFeaturesOutNaive);

figure, hold on
bar(1:12, mean(eachFeatureOut), 'k')
errorbar(1:12, mean(eachFeatureOut), sem(eachFeatureOut), 'k', 'lines', 'no')
ylim([-0.001 0.2])
yticks(0:.05:0.2)
xticks(1:12)
xticklabels({'maxDq', 'maxDf', 'maxDkH', 'maxDkV', 'Slide distance', 'Touch duration', ...
    'q', 'f', 'kH', 'kV', 'Arc length', 'Touch count'})
xtickangle(45)
ylabel('Impact on tuning curve')
title('Naive')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%% Expert
eachFeatureOut = cellfun(@mean, corrWhiskerExpert) - cellfun(@mean, corrFeaturesOutExpert);

figure, hold on
bar(1:12, mean(eachFeatureOut), 'k')
errorbar(1:12, mean(eachFeatureOut), sem(eachFeatureOut), 'k', 'lines', 'no')
ylim([-0.001 0.2])
yticks(0:.05:0.2)
xticks(1:12)
xticklabels({'maxDq', 'maxDf', 'maxDkH', 'maxDkV', 'Slide distance', 'Touch duration', ...
    'q', 'f', 'kH', 'kV', 'Arc length', 'Touch count'})
xtickangle(45)
ylabel('Impact on tuning curve')
title('Expert')
set(gca, 'fontsize', 12, 'fontname', 'Arial')
%%
% saveDir = 'C:\Users\jinho\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig5-S2 Variable importance and impact on angle tuning\';
% fn = 'each_feature_impacts.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
% fix_eps_fonts([saveDir, fn])

%%
%% Fig 5E
%%
% Distribution of maximum single feature impact
maxImpact = cell(numMiceNaive,1);
for mi = 1 : numMiceNaive
    maxImpact{mi} = max(corrWhiskerNaive{mi} - cell2mat(corrFeaturesOutNaive(mi,:)),[],2);
end
maxImpactDist = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'prob'), maxImpact, 'un', 0)); 
% figure,
%boundedline(histRange(2:end)-0.05, mean(maxImpactDist), sem(maxImpactDist), 'b')

%% Normalized by whisker correlation
maxReductionProportion = cellfun(@(x,y) y./x, corrWhiskerNaive, maxImpact, 'un', 0);
maxRPdist = cell2mat(cellfun(@(x) histcounts(x,histRange, 'norm', 'prob'), maxReductionProportion, 'un', 0));
%figure,
%boundedline(histRange(2:end)-0.05, mean(maxRPdist), sem(maxRPdist))

%% Distribution of all feature impact
allImpact = cell(numMiceNaive,1);
for mi = 1 : numMiceNaive
    allImpact{mi} = max(corrWhiskerNaive{mi} - corrOtherNaive{mi}',[],2);
end
allImpactDist = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'prob'), allImpact, 'un', 0)); 
% figure,
% boundedline(histRange(2:end)-0.05, mean(allImpactDist), sem(allImpactDist), 'r')


%% Compare mean impacts of 
% mean whisker feature, each candidate features, max single feature, all whisker features
tempFeatureOut = cellfun(@mean, corrFeaturesOutNaive(:, [10,2,5,4]));

allImpactMean = cellfun(@mean, allImpact);
maxImpactMean = cellfun(@mean, maxImpact);
temptempFeatureOutMean = cellfun(@mean, corrFeaturesOutNaive);
tempFeatureOutMean = mean(temptempFeatureOutMean,2);
meanFeatureImpact = cellfun(@mean, corrWhiskerNaive) - tempFeatureOutMean;
eachFeatureImpact = cellfun(@mean, corrWhiskerNaive) - tempFeatureOut;

bothFeatureImpact = cellfun(@mean, corrWhiskerNaive) - cellfun(@mean, corrFeaturesCombOutNaive(:,3));
allBut2FeatureImpact = cellfun(@mean, corrWhiskerNaive) - cellfun(@mean, corrFeaturesCombInNaive(:,3));

figure, hold on
bar(1, mean(allImpactMean), 'facecolor', ones(1,3)*0.7)
bar(2, mean(maxImpactMean), 'facecolor', ones(1,3)*0.7)
bar(3, mean(meanFeatureImpact), 'facecolor', ones(1,3)*0.7)
bar(4:7, mean(eachFeatureImpact), 'k')
bar(8, mean(bothFeatureImpact), 'k')
bar(9, mean(allBut2FeatureImpact), 'k')

errorbar(1, mean(allImpactMean), sem(allImpactMean), 'k')
errorbar(2, mean(maxImpactMean), sem(maxImpactMean), 'k')
errorbar(3, mean(meanFeatureImpact), sem(meanFeatureImpact), 'k')
errorbar(4:7, mean(eachFeatureImpact), sem(eachFeatureImpact),'k', 'lines', 'no')
errorbar(8, mean(bothFeatureImpact), sem(bothFeatureImpact), 'k')
errorbar(9, mean(allBut2FeatureImpact), sem(allBut2FeatureImpact), 'k')

xticks([1:9])

xticklabels({'All features', 'Max feature', 'Mean features', '-kV', '-maxDphi','-max(slide distance)', '-maxDkV', '-(slide distance & maxDkV)', '-All features except (slide distance & maxDkV)'})
xtickangle(45)
ylabel('Impact on tuning curve')
set(gca, 'fontsize', 12, 'fontname','Arial', 'box', 'off');
%%

[~, p] = ttest(cellfun(@mean, corrWhiskerNaive), cellfun(@mean, allImpact))
mean(cellfun(@mean, allImpact))
sem(cellfun(@mean, allImpact))
%%
mean(cellfun(@mean, maxImpact))
sem(cellfun(@mean, maxImpact))

%%
t = array2table([maxImpactMean,meanFeatureImpact,eachFeatureImpact, bothFeatureImpact, allBut2FeatureImpact], 'VariableNames', {'v1','v2','v3','v4','v5','v6','v7','v8'});
rm = fitrm(t,'v1-v8~1');
tb1 = mauchly(rm);
tb2 = ranova(rm);
if tb1.pValue(1) < 0.05
    pval = tb2.pValueGG(1);
else
    pval = tb2.pValue(1);
end

multcompare(rm, 'Time')

%%
% saveDir = 'C:\Users\jinho\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig5 Sensory input for angle tuning\';
% fn = 'impacts_on_tuning_curve.eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
% fix_eps_fonts([saveDir, fn])


%%
%% Fig 5F
%%
%%
distMethod = 'cdf';
histBin = 0.05;
histRange = [-1,0:histBin:1];
corrWhiskerDist = zeros(numMiceNaive, length(histRange)-1);
corrFeaturesOutDist = zeros(numMiceNaive, length(histRange)-1, 12);
corrFeaturesOutDropDist = zeros(numMiceNaive, length(histRange)-1, 12);
corrOtherDropDist = zeros(numMiceNaive, length(histRange)-1);
corrFeaturesCombInDropDist = zeros(numMiceNaive, length(histRange)-1); % only for 4 and 5
corrFeaturesCombOutDropDist = zeros(numMiceNaive, length(histRange)-1); % only for 4 and 5

for i = 1 : numMiceNaive
    corrWhiskerDist(i,:) = histcounts(corrWhiskerNaive{i}, histRange, 'norm', distMethod);
    for j = 1 : 12
        corrFeaturesOutDist(i,:,j) = histcounts(corrFeaturesOutNaive{i,j}, histRange, 'norm', distMethod);
        corrFeaturesOutDropDist(i,:,j) = histcounts(corrWhiskerNaive{i} - corrFeaturesOutNaive{i,j}, histRange, 'norm', distMethod);
    end
    corrOtherDropDist(i,:) = histcounts(corrWhiskerNaive{i} - corrOtherNaive{i}, histRange, 'norm', distMethod);
    corrFeaturesCombInDropDist(i,:) = histcounts(corrWhiskerNaive{i} - corrFeaturesCombInNaive{i,3}, histRange, 'norm', distMethod);
    corrFeaturesCombOutDropDist(i,:) = histcounts(corrWhiskerNaive{i} - corrFeaturesCombOutNaive{i,3}, histRange, 'norm', distMethod);
end
avgCorrDropDist = mean(corrFeaturesOutDropDist,3);

%%
figure, hold on
boundedline(histRange(2:end), mean(squeeze(corrFeaturesOutDropDist(:,:,4))), sem(squeeze(corrFeaturesOutDropDist(:,:,4))), '-g')
boundedline(histRange(2:end), mean(squeeze(corrFeaturesOutDropDist(:,:,5))), sem(squeeze(corrFeaturesOutDropDist(:,:,5))), '-b')
boundedline(histRange(2:end), mean(avgCorrDropDist), sem(avgCorrDropDist), '-k')
boundedline(histRange(2:end), mean(corrFeaturesCombOutDropDist), sem(corrFeaturesCombOutDropDist), '-m')
boundedline(histRange(2:end), mean(corrFeaturesCombInDropDist), sem(corrFeaturesCombInDropDist), '-y')
xlim([0 1])
xlabel('Correlation value drop')
ylabel('Cumulative distribution')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%%
%% Fig S8I
%%
tunedAngleDiffOut = cell(numMiceNaive, 15); % 1 for touch model, 2 for full whisker model, 3:14 for each individual whisker features, 15 for no whisker features
tunedAngleDiffIn = cell(numMiceNaive, 12); % 1:12 for each individual whisker features, + others
tunedAngleDiffCombOut = cell(numMiceNaive, 4); %
tunedAngleDiffCombIn = cell(numMiceNaive, 4); %
for i = 1 : 12
    indTuned = find(data1.naive(i).tunedAllCell(:,1));
    indTemp = find(data1.naive(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    tunedAngleDiffOut{i,1} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data1.naive(i).tuneAngleAllCell(indWhisker,2); 
    tunedAngleDiffOut{i,2} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data1.naive(i).tuneAngleAllCell(indWhisker,3); 
    
    for fi = 1 : 12
        tunedAngleDiffOut{i,2+fi} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data1.naive(i).tuneAngleAllCell(indWhisker, 3+fi);
        tunedAngleDiffIn{i,fi} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data2.naive(i).tuneAngleAllCell(indWhisker, 1+fi); 
    end
    tunedAngleDiffOut{i,15} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data2.naive(i).tuneAngleAllCell(indWhisker,1); 
    
    for fi = 1 : 4
        tunedAngleDiffCombOut{i,fi} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data2.naive(i).tuneAngleAllCell(indWhisker, 13+fi); 
        tunedAngleDiffCombIn{i,fi} = data1.naive(i).tuneAngleAllCell(indWhisker,1) - data2.naive(i).tuneAngleAllCell(indWhisker, 21+fi); 
    end
end

fullMatch = cellfun(@(x) length(find(abs(x) <=15))/length(x), tunedAngleDiffOut(:,2));
eachMatch = cellfun(@(x) length(find(abs(x) <=15))/length(x), tunedAngleDiffOut(:,3:14));

figure, hold on
bar(1, 1-mean(fullMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, 1-mean(fullMatch), sem(fullMatch), 'color', ones(1,3)*0.7)
bar(2:13, 1-mean(eachMatch), 'k')
errorbar(2:13, 1-mean(eachMatch), sem(eachMatch), 'k', 'lines','no')
xticks([1:13]), xticklabels({'Full model', 'MaxDq', 'MaxDf', 'MaxDkH', 'MaxDkV', 'Slide distance', 'Touch duration', ...
    'q', 'f', 'kH', 'kV', 'Arc length', 'Touch count'})

xtickangle(45)
ylabel('Prop. tuning changed')
title('Naive: Tuned & \Deltaangle <= 15\circ')

%% Expert
tunedAngleDiffOutExpert = cell(numMiceExpert, 15); % 1 for touch model, 2 for full whisker model, 3:14 for each individual whisker features, 15 for no whisker features
tunedAngleDiffInExpert = cell(numMiceExpert, 12); % 1:12 for each individual whisker features, + others
tunedAngleDiffCombOutExpert = cell(numMiceExpert, 4); %
tunedAngleDiffCombInExpert = cell(numMiceExpert, 4); %
for i = 1 : numMiceExpert
    indTuned = find(data1.expert(i).tunedAllCell(:,1));
    indTemp = find(data1.expert(i).tunedAllCell(indTuned,3));
    indWhisker = indTuned(indTemp);
    
    tunedAngleDiffOutExpert{i,1} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data1.expert(i).tuneAngleAllCell(indWhisker,2); 
    tunedAngleDiffOutExpert{i,2} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data1.expert(i).tuneAngleAllCell(indWhisker,3); 
    
    for fi = 1 : 12
        tunedAngleDiffOutExpert{i,2+fi} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data1.expert(i).tuneAngleAllCell(indWhisker, 3+fi);
        tunedAngleDiffInExpert{i,fi} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data2.expert(i).tuneAngleAllCell(indWhisker, 1+fi); 
    end
    tunedAngleDiffOutExpert{i,15} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data2.expert(i).tuneAngleAllCell(indWhisker,1); 
    
    for fi = 1 : 4
        tunedAngleDiffCombOutExpert{i,fi} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data2.expert(i).tuneAngleAllCell(indWhisker, 13+fi); 
        tunedAngleDiffCombInExpert{i,fi} = data1.expert(i).tuneAngleAllCell(indWhisker,1) - data2.expert(i).tuneAngleAllCell(indWhisker, 21+fi); 
    end
end

fullMatch = cellfun(@(x) length(find(abs(x) <=15))/length(x), tunedAngleDiffOutExpert(:,2));
eachMatch = cellfun(@(x) length(find(abs(x) <=15))/length(x), tunedAngleDiffOutExpert(:,3:14));

figure, hold on
bar(1, 1-mean(fullMatch), 'facecolor', ones(1,3)*0.7)
errorbar(1, 1-mean(fullMatch), sem(fullMatch), 'color', ones(1,3)*0.7)
bar(2:13, 1-mean(eachMatch), 'k')
errorbar(2:13, 1-mean(eachMatch), sem(eachMatch), 'k', 'lines','no')
xticks([1:13]), xticklabels({'Full model', 'MaxDq', 'MaxDf', 'MaxDkH', 'MaxDkV', 'Slide distance', 'Touch duration', ...
    'q', 'f', 'kH', 'kV', 'Arc length', 'Touch count'})

xtickangle(45)
ylabel('Prop. tuning changed')
title('Expert: Tuned & \Deltaangle <= 15\circ')

%%
% saveDir = 'C:\Users\jinho\Dropbox\Works\Manuscripts\Object Angle Coding in vS1\Figures\Fig5-S2 Variable importance and impact on angle tuning\';
% fn = 'prop_tuned_angle_changed(15d).eps';
% export_fig([saveDir, fn], '-depsc', '-painters', '-r600', '-transparent')
% fix_eps_fonts([saveDir, fn])





%% Distribution of 1st, 2nd, and 3rd features in impact on tuning
impactNaive = cell(numMiceNaive, 12);
for mi = 1 : numMiceNaive
    for fi = 1 : 12
        impactNaive{mi,fi} = corrWhiskerNaive{mi} - corrFeaturesOutNaive{mi,fi};
    end
end

%%
impactSortedNaive = cell(numMiceNaive,1);
impactSortiNaive = cell(numMiceExpert,1);
for mi = 1 : numMiceNaive
    tempMat = cell2mat(impactNaive(mi,:));
    [impactSortedNaive{mi}, impactSortiNaive{mi}] = sort(tempMat,2,'descend');
end
%%
histRange = 0:0.02:1;
distImpactNaive = cell(numMiceNaive,1);
for mi = 1 : numMiceNaive
    distImpactNaive{mi} = zeros(length(histRange)-1, 12);
    for fi = 1 : 12
        distImpactNaive{mi}(:,fi) = histcounts(impactSortedNaive{mi}(:,fi), histRange,'norm','cdf');
    end
end
%%
threshold = 0.1;
figure
for i = 1 : 4
    tempMat = zeros(numMiceNaive, length(histRange)-1);
    for mi = 1 : numMiceNaive
        tempMat(mi,:) = distImpactNaive{mi}(:,i)';
    end
    subplot(2,2,i), hold on
    boundedline(histRange(1:end-1), mean(tempMat), sem(tempMat), 'k')
    plot([threshold threshold], [0 1], 'r--')
    ind = find(histRange>=threshold, 1, 'first');
    val = mean(tempMat(:,ind));
    ylim([0 1])
    xticks([0, threshold, 0.5, 1])
    title({sprintf('Order = %d', i); sprintf('> threshold = %.2f', 1-val)})
end
    
sgtitle('Naive')

%%
impactExpert = cell(numMiceExpert, 12);
for mi = 1 : numMiceExpert
    for fi = 1 : 12
        impactExpert{mi,fi} = corrWhiskerExpert{mi} - corrFeaturesOutExpert{mi,fi};
    end
end

impactSortedExpert = cell(numMiceExpert,1);
impactSortiExpert = cell(numMiceExpert,1);
for mi = 1 : numMiceExpert
    tempMat = cell2mat(impactExpert(mi,:));
    [impactSortedExpert{mi}, impactSortiExpert{mi}] = sort(tempMat,2,'descend');
end
histRange = 0:0.02:1;
distImpactExpert = cell(numMiceExpert,1);
for mi = 1 : numMiceExpert
    distImpactExpert{mi} = zeros(length(histRange)-1, 12);
    for fi = 1 : 12
        distImpactExpert{mi}(:,fi) = histcounts(impactSortedExpert{mi}(:,fi), histRange,'norm','cdf');
    end
end

%%
threshold = 0.1;
figure
for i = 1 : 4
    tempMat = zeros(numMiceExpert, length(histRange)-1);
    for mi = 1 : numMiceExpert
        tempMat(mi,:) = distImpactExpert{mi}(:,i)';
    end
    subplot(2,2,i), hold on
    boundedline(histRange(1:end-1), mean(tempMat), sem(tempMat), 'k')
    plot([threshold threshold], [0 1], 'r--')
    ind = find(histRange>=threshold, 1, 'first');
    val = mean(tempMat(:,ind));
    ylim([0 1])
    xticks([0, threshold, 0.5, 1])
    title({sprintf('Order = %d', i); sprintf('> threshold = %.2f', 1-val)})
end
    
sgtitle('Expert')



%% Layer specificity
% the ones with second impact feature > threshold, where are they?
% Record their ID, to compare between naive and expert

tune = load([baseDir, 'angle_tuning_summary_predecision_NC'], 'naive', 'expert');

%% naive
threshold = 0.1;
multiFeatureIdNaive = cell(numMiceNaive,1);
multiFeatureDepthNaive = cell(numMiceNaive,1);
for mi = 1 : numMiceNaive
    tempInd = impactSortedNaive{mi}(:,2) > threshold;
    multiFeatureIdNaive{mi} = tune.naive(mi).touchID(corrIndsNaive{mi}(tempInd));
    multiFeatureDepthNaive{mi} = tune.naive(mi).depth(corrIndsNaive{mi}(tempInd));
end


%%
multiL23propNaive = cellfun(@(x) length(find(x<350))/length(x), multiFeatureDepthNaive);
mean(multiL23propNaive)
sem(multiL23propNaive)

%% Expert

multiFeatureIdExpert = cell(numMiceExpert,1);
multiFeatureDepthExpert = cell(numMiceExpert,1);
for mi = 1 : numMiceExpert
    tempInd = impactSortedExpert{mi}(:,2) > threshold;
    multiFeatureIdExpert{mi} = tune.expert(mi).touchID(corrIndsExpert{mi}(tempInd));
    multiFeatureDepthExpert{mi} = tune.expert(mi).depth(corrIndsExpert{mi}(tempInd));
end


%%
multiL23propExpert = cellfun(@(x) length(find(x<350))/length(x), multiFeatureDepthExpert);
mean(multiL23propExpert)
sem(multiL23propExpert)



%% How these multi feature neurons change during learning?
% Follow-up multi-feature angle-tuned neurons from naive sessions.
% Calculate their angle selectivity

load([baseDir, 'cellMatching_beforeNafter'], 'match');
learnerInd = [1:4,7,9];
multiSelectivity = cell(length(learnerInd),1);
multiPersistency = zeros(length(learnerInd),1);
for mi = 1 : length(learnerInd)
    multiId = multiFeatureIdNaive{learnerInd(mi)};
    multiInd = find(ismember(match{mi,1}, multiId));
    multiExpertId = match{mi,2}(multiInd);
    numMultiTransient = length(find(multiExpertId == 0));
    indMultiPersistent = find(match{mi,2}(multiInd));
    multiPersistency(mi) = 1 - numMultiTransient / length(multiId);
    multiSelectivity{mi} = zeros(length(indMultiPersistent),2);
    
    for ci = 1 : length(indMultiPersistent)
        tempInd = multiInd(indMultiPersistent(ci));
        idNaive = match{mi,1}(tempInd);
        indNaive = find(tune.naive(learnerInd(mi)).touchID == idNaive);
        selectivityNaive = cellfun(@(x) max(cellfun(@mean, x))*7/6 - sum(cellfun(@mean,x))/6, tune.naive(learnerInd(mi)).val(indNaive));
        idExpert = match{mi,2}(tempInd);
        indExpert = find(tune.expert(mi).touchID == idExpert);
        selectivityExpert = cellfun(@(x) max(cellfun(@mean, x))*7/6 - sum(cellfun(@mean,x))/6, tune.expert(mi).val(indExpert));
        multiSelectivity{mi}(ci,:) = [selectivityNaive, selectivityExpert];
    end
end


%%
mean(multiPersistency)
sem(multiPersistency)

%%
selectivityChange = cell2mat(cellfun(@(x) x(:,2) - x(:,1), multiSelectivity, 'un', 0));
mean(selectivityChange)
sem(selectivityChange)
[~,p,m] = paired_test(selectivityChange)


%% Examples of multi-feature angle-tuning
% find max 2nd feature impact from all neurons
% from impactSortedNaive & corrIndsNaive
all2ndImpactNaive = cell2mat(cellfun(@(x) x(:,2), impactSortedNaive, 'un', 0));
[~,sortInd] = sort(all2ndImpactNaive, 'descend');
numEach = cellfun(@length, impactSortedNaive);
numCumsum = cumsum(numEach);
%%
i = 9;
mi = find(sortInd(i)<numCumsum,1);
tunedCi = sortInd(i)-numCumsum(mi-1);
ci = corrIndsNaive{mi}(tunedCi);
cellSorti = impactSortiNaive{mi}(tunedCi,:);
%%
angles = 45:15:135;
spikes = (cellfun(@mean, data1.naive(mi).spkValAllCell{ci,1}));
whiskerGLM = (cellfun(@mean, data1.naive(mi).spkValAllCell{ci,3}));
wo1st = (cellfun(@mean, data1.naive(mi).spkValAllCell{ci,3+cellSorti(1)}));
wo2dn = (cellfun(@mean, data1.naive(mi).spkValAllCell{ci,3+cellSorti(2)}));

figure, hold on
ax1 = plot(angles, spikes, 'ro-', 'MarkerFaceColor','r', 'linewidth', 1);
ylabel('DSpikes/touch')
yyaxis right
ax2 = plot(angles, whiskerGLM, 'ko-', 'MarkerFaceColor','k', 'linewidth', 1);
plot(angles, wo1st, 'mo--', 'linewidth', 1)
plot(angles, wo2dn, 'co--', 'linewidth', 1)
ylabel('DSpikes/touch (model)')
xticks(angles)
xlabel('Object angle (\circ)')
legend({'Inferred spikes', 'Full model', '-maxDkV', '-max(slide distance)'}, 'location', 'northeast', 'box', 'off')
set(ax1.Parent, 'ycolor', 'k')
yyaxis left
set(ax1.Parent, 'ycolor', 'r')


corrVal = zeros(1,3);
corrVal(1) = corr(spikes, whiskerGLM);
corrVal(2) = corr(spikes, wo1st);
corrVal(3) = corr(spikes, wo2dn);
figure, hold on
bar(1, corrVal(1), 'k')
bar(2, corrVal(2), 'm')
bar(3, corrVal(3), 'c')
ylabel('Correlation')
xticks([1:3])
xticklabels({'Full model', '-maxDkv', '-max(slide distance)'})
xtickangle(45)

%%
% cellfun(@mean, data1.naive(mi).spkValAllCell{ci,3+fi})
cellfun(@mean, data1.naive(mi).spkValAllCell{ci,3+fi})






















































%%
%%
%% 5E - Multiple whisker features for angle-tuning
% How many whisker features contribute to angle-tuning? (S4D)
% Where do they reside (L2/3 vs L4) (S4E)
% How do they differ across tuned angles? (4E)

%%
loadFn = 'wkv_angle_tuning_model_v9';
load([calciumDir, loadFn])

learnerInd = [1:4,7,9];
naiveLearner = naive(learnerInd);
jk027inds = find(naiveLearner(2).cIDAll > 5000);
jk027AllLength = length(naiveLearner(2).cIDAll);

fn = fieldnames(naiveLearner);
for i = 1 : length(fn)
    if length(naiveLearner(2).(fn{i})) == jk027AllLength
        naiveLearner(2).(fn{i}) = naiveLearner(2).(fn{i})(jk027inds);
    end
end

angles = 45:15:135;

atCorrThreshold = 0.1; % any feature that reduces the correlation more than this amount is considered important in angle-tuning
depthThreshold = 350; % 

%% S4D
%% Cumulative distribution of neurons fit with multiple features
sortedFeatureImportance = cell(length(naive), 2); %(:,1) from full model, (:,2) from whisker-only model
sortedFeatureInd = cell(length(naive),2);
cellDepths = cell(length(naive),1);
cellTunedAngle = cell(length(naive),2); %(:,1) from spikes, (:,2) from whisker-only model

for mi = 1 : length(naive)
    tempInd = find(naive(mi).tuned);
    cellDepths{mi} = naive(mi).depth(tempInd);
    tempSpikeAll = cell2mat(naive(mi).atSpikeAll(tempInd));
    [~, tempAngleInd] = max(tempSpikeAll, [], 2);
    cellTunedAngle{mi,1} = 30 + tempAngleInd * 15;
    tempWhiskerAll = cell2mat(cellfun(@(x) x(1,:), naive(mi).atWhiskerAll(tempInd), 'un', 0));
    [~, tempAngleInd] = max(tempWhiskerAll, [], 2);
    cellTunedAngle{mi,2} = 30 + tempAngleInd * 15;
    
    tempMat = cell2mat(naive(mi).atCorrFull(tempInd));
    [sortedFeatureImportance{mi,1}, sortedFeatureInd{mi,1}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
    tempMat = cell2mat(naive(mi).atCorrWhisker(tempInd));
    [sortedFeatureImportance{mi,2}, sortedFeatureInd{mi,2}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
end

% distribution 1st to 3rd feature importance
% and L2/3prop
histRange = [-1, 0:0.02:1.8];
distFull = zeros(length(naive),length(histRange)-1, 3);
distWhisker = zeros(length(naive),length(histRange)-1, 3);
L23prop = zeros(length(naive),2); % (:,1) full model, (:,2) whisker-only model
L4prop = zeros(length(naive),2);
for mi = 1 : length(naive)
    for i = 1 : 3
        distFull(mi,:,i) = histcounts(sortedFeatureImportance{mi,1}(:,i), histRange, 'norm', 'cdf');
        distWhisker(mi,:,i) = histcounts(sortedFeatureImportance{mi,2}(:,i), histRange, 'norm', 'cdf');
    end
    for i = 1 : 2
        multiInd = find(sortedFeatureImportance{mi,i}(:,2) > atCorrThreshold);
        L23prop(mi,i) = length(find(cellDepths{mi}(multiInd) <= depthThreshold)) / length(find(cellDepths{mi} <= depthThreshold));
        L4prop(mi,i) = length(find(cellDepths{mi}(multiInd) > depthThreshold)) / length(find(cellDepths{mi} > depthThreshold));
    end
end

% Draw
colors = {'k','r','b'};
figure, hold on
props = zeros(1,3);
sems = zeros(1,3);
for i = 1 : 3
    tempMat = squeeze(distWhisker(:,:,i));
    plot(histRange(2:end), mean(tempMat), 'color', colors{i})
    props(i) = mean(tempMat(:,find(histRange == atCorrThreshold)-1));
    sems(i) = sem(tempMat(:,find(histRange == atCorrThreshold)-1));
end
legend({'1^s^t', '2^n^d', '3^r^d'}, 'location','southeast','autoupdate', false)
for i = 1 : 3
    tempMat = squeeze(distWhisker(:,:,i));
    boundedline(histRange(2:end),mean(tempMat), sem(tempMat), colors{i})
end
plot([atCorrThreshold atCorrThreshold], [0 1], '--', 'color', [0.6 0.6 0.6])
xlabel('Impact on angle-tuning')
xticks([0 0.1 0.5 1 1.5 1.8])
ylabel('Cumulative proportion')
title({'Naive all (n = 12)'; sprintf('Props = %2.1f \\pm %2.1f, %2.1f \\pm %2.1f, %2.1f \\pm %2.1f %%', ...
    (1-props(1)) * 100, sems(1)*100, (1-props(2)) * 100, sems(2)*100, (1-props(3)) * 100, sems(3)*100 ); ...
    sprintf('Prop L2/3 = %2.1f \\pm %2.1f %%', mean(L23prop(:,2))*100, sem(L23prop(:,2)) * 100); ...
    sprintf('Prop L4 = %2.1f \\pm %2.1f %%', mean(L4prop(:,2))*100, sem(L4prop(:,2)) * 100)})





%% 4E
%% Proportion of multi-whisker angle-tuned neurons across tuned angles
propMultiAnglePerSpike = zeros(length(naive), length(angles)); 

for mi = 1 : length(naive)
    multiInd = find(sortedFeatureImportance{mi,2}(:,2) > atCorrThreshold);
    for ai = 1: length(angles)
        angleIndSpike = find(cellTunedAngle{mi,1} == angles(ai)); % perSpike
        propMultiAnglePerSpike(mi,ai) = length(intersect(multiInd, angleIndSpike)) / length(angleIndSpike);
     end
end

figure, hold on
bar(angles, nanmean(propMultiAnglePerSpike), 'facecolor', 'k')
errorbar(angles, nanmean(propMultiAnglePerSpike), sem(propMultiAnglePerSpike),'k', 'lines', 'no')
xlim([35 145]), xticks(angles), xlabel('Tuned angle (\circ)')
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel({'Proportion';'(/ each preferred angle bin)'})
title('Multi-feature neurons naive all (n = 12)')

%% ANOVA
anovaVal = propMultiAnglePerSpike(:);
anovaGroup = zeros(length(anovaVal),1);
for ai = 1 : length(angles)
    anovaGroup(length(naive)*(ai-1)+1:length(naive)*ai) = deal(ai);
end
[anovaP, anovaTbl, stats] = anova1(anovaVal, anovaGroup, 'off');

%% paired test
middleVal = mean(propMultiAnglePerSpike(:,3:5),2);
edgeVal = mean(propMultiAnglePerSpike(:,[1,2,6,7]),2);
[~,p,m] = paired_test(middleVal, edgeVal);




%% Extra
% proportion of multi-whisker feature neurons 
% across layers and learning

% naive learner (n = 6)
sortedFeatureImportanceNaive = cell(length(naiveLearner), 2); %(:,1) from full model, (:,2) from whisker-only model
sortedFeatureIndNaive = cell(length(naiveLearner),2);
cellDepthsNaive = cell(length(naiveLearner),1);
cellTunedAngleNaive = cell(length(naiveLearner),2); %(:,1) from spikes, (:,2) from whisker-only model

for mi = 1 : length(naiveLearner)
    tempInd = find(naiveLearner(mi).tuned);
    cellDepthsNaive{mi} = naiveLearner(mi).depth(tempInd);
    tempSpikeAll = cell2mat(naiveLearner(mi).atSpikeAll(tempInd));
    [~, tempAngleInd] = max(tempSpikeAll, [], 2);
    cellTunedAngleNaive{mi,1} = 30 + tempAngleInd * 15;
    tempWhiskerAll = cell2mat(cellfun(@(x) x(1,:), naiveLearner(mi).atWhiskerAll(tempInd), 'un', 0));
    [~, tempAngleInd] = max(tempWhiskerAll, [], 2);
    cellTunedAngleNaive{mi,2} = 30 + tempAngleInd * 15;
    
    tempMat = cell2mat(naiveLearner(mi).atCorrFull(tempInd));
    [sortedFeatureImportanceNaive{mi,1}, sortedFeatureIndNaive{mi,1}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
    tempMat = cell2mat(naiveLearner(mi).atCorrWhisker(tempInd));
    [sortedFeatureImportanceNaive{mi,2}, sortedFeatureIndNaive{mi,2}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
end
%
% distribution 1st to 3rd feature importance
% and L2/3prop
histRange = 0:0.02:1;
distFullNaive = zeros(length(naiveLearner),length(histRange)-1, 3);
distWhiskerNaive = zeros(length(naiveLearner),length(histRange)-1, 3);
allPropNaive = zeros(length(naiveLearner),2);
L23propNaive = zeros(length(naiveLearner),2); % (:,1) full model, (:,2) whisker-only model
L4propNaive = zeros(length(naiveLearner),2);
for mi = 1 : length(naiveLearner)
    for i = 1 : 3
        distFullNaive(mi,:,i) = histcounts(sortedFeatureImportanceNaive{mi,1}(:,i), histRange, 'norm', 'cdf');
        distWhiskerNaive(mi,:,i) = histcounts(sortedFeatureImportanceNaive{mi,2}(:,i), histRange, 'norm', 'cdf');
    end
    for i = 1 : 2
        multiInd = find(sortedFeatureImportanceNaive{mi,i}(:,2) > atCorrThreshold);
        allPropNaive(mi,i) = length(multiInd) / size(sortedFeatureImportanceNaive{mi,i},1);
        L23propNaive(mi,i) = length(find(cellDepthsNaive{mi}(multiInd) <= depthThreshold)) / length(find(cellDepthsNaive{mi} <= depthThreshold));
        L4propNaive(mi,i) = length(find(cellDepthsNaive{mi}(multiInd) > depthThreshold)) / length(find(cellDepthsNaive{mi} > depthThreshold));
    end
end

% expert (n = 6)
sortedFeatureImportanceExpert = cell(length(expert), 2); %(:,1) from full model, (:,2) from whisker-only model
sortedFeatureIndExpert = cell(length(expert),2);
cellDepthsExpert = cell(length(expert),1);
cellTunedAngleExpert = cell(length(expert),2); %(:,1) from spikes, (:,2) from whisker-only model

for mi = 1 : length(expert)
    tempInd = find(expert(mi).tuned);
    cellDepthsExpert{mi} = expert(mi).depth(tempInd);
    tempSpikeAll = cell2mat(expert(mi).atSpikeAll(tempInd));
    [~, tempAngleInd] = max(tempSpikeAll, [], 2);
    cellTunedAngleExpert{mi,1} = 30 + tempAngleInd * 15;
    tempWhiskerAll = cell2mat(cellfun(@(x) x(1,:), expert(mi).atWhiskerAll(tempInd), 'un', 0));
    [~, tempAngleInd] = max(tempWhiskerAll, [], 2);
    cellTunedAngleExpert{mi,2} = 30 + tempAngleInd * 15;
    
    tempMat = cell2mat(expert(mi).atCorrFull(tempInd));
    [sortedFeatureImportanceExpert{mi,1}, sortedFeatureIndExpert{mi,1}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
    tempMat = cell2mat(expert(mi).atCorrWhisker(tempInd));
    [sortedFeatureImportanceExpert{mi,2}, sortedFeatureIndExpert{mi,2}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
end
%
% distribution 1st to 3rd feature importance
% and L2/3prop
histRange = 0:0.02:1;
distFullExpert = zeros(length(expert),length(histRange)-1, 3);
distWhiskerExpert = zeros(length(expert),length(histRange)-1, 3);
allPropExpert = zeros(length(expert),2);
L23propExpert = zeros(length(expert),2); % (:,1) full model, (:,2) whisker-only model
L4propExpert = zeros(length(expert),2);
for mi = 1 : length(expert)
    for i = 1 : 3
        distFullExpert(mi,:,i) = histcounts(sortedFeatureImportanceExpert{mi,1}(:,i), histRange, 'norm', 'cdf');
        distWhiskerExpert(mi,:,i) = histcounts(sortedFeatureImportanceExpert{mi,2}(:,i), histRange, 'norm', 'cdf');
    end
    for i = 1 : 2
        multiInd = find(sortedFeatureImportanceExpert{mi,i}(:,2) > atCorrThreshold);
        allPropExpert(mi,i) = length(multiInd) / size(sortedFeatureImportanceExpert{mi,i},1);
        L23propExpert(mi,i) = length(find(cellDepthsExpert{mi}(multiInd) <= depthThreshold)) / length(find(cellDepthsExpert{mi} <= depthThreshold));
        L4propExpert(mi,i) = length(find(cellDepthsExpert{mi}(multiInd) > depthThreshold)) / length(find(cellDepthsExpert{mi} > depthThreshold));
    end
end


% Draw
figure('units','inch','pos',[3 3 3 5]), hold on
for mi = 1 : 6
    plot([allPropNaive(mi,2), allPropExpert(mi,2)], 'ko-')
end
errorbar(1, mean(allPropNaive(:,2)), sem(allPropNaive(:,2)), 'ro')
errorbar(2, mean(allPropExpert(:,2)), sem(allPropExpert(:,2)), 'ro')
xlim([0.5 2.5]), xticks(1:2), xticklabels({'Naive', 'Expert'})
ylim([0 0.3]), yticks(0:0.1:0.3)
ylabel('Proportion')

[~,p,m] = paired_test(allPropNaive(:,2), allPropExpert(:,2));
title(sprintf('p = %s (%s)', num2str(p,3), m))




%% bar graph
colorsTransient = [248 171 66; 40 170 225] / 255;

barOffset = 0.2;
barWidth = 0.4;
xpos = 1:2;

p = zeros(1,2);
m = cell(1,2);

figure('units','norm','pos',[0.1 0.3 0.2 0.3]), hold on
bar(1-barOffset, mean(L23propNaive(:,2)), barWidth, 'facecolor', colorsTransient(1,:), 'edgecolor', 'none');
bar(1+barOffset, mean(L4propNaive(:,2)), barWidth, 'facecolor', 'w', 'edgecolor', colorsTransient(1,:));
legend({'L2/3', 'L4'}, 'autoupdate', false)
errorbar(1-barOffset, mean(L23propNaive(:,2)), sem(L23propNaive(:,2)), 'k')
errorbar(1+barOffset, mean(L4propNaive(:,2)), sem(L4propNaive(:,2)), 'k')
[~,p(1),m{1}] = paired_test(L23propNaive(:,2), L4propNaive(:,2));

bar(2-barOffset, mean(L23propExpert(:,2)), barWidth, 'facecolor', colorsTransient(2,:), 'edgecolor', 'none');
bar(2+barOffset, mean(L4propExpert(:,2)), barWidth, 'facecolor', 'w', 'edgecolor', colorsTransient(2,:));
errorbar(2-barOffset, mean(L23propExpert(:,2)), sem(L23propExpert(:,2)), 'k')
errorbar(2+barOffset, mean(L4propExpert(:,2)), sem(L4propExpert(:,2)), 'k')
[~,p(2),m{2}] = paired_test(L23propExpert(:,2), L4propExpert(:,2));

xlim([0.2 2.8]), xticks(1:2), xticklabels({'Naive', 'Expert'})
ylabel('Proportion')
title(sprintf('p = %.3f (m = %s); p = %.3f (m = %s)',p(1), m{1}, p(2), m{2}))







%%
%% S4E - From all naive mice
%%

sortedFeatureImportanceNaiveAll = cell(length(naive), 2); %(:,1) from full model, (:,2) from whisker-only model
sortedFeatureIndNaive = cell(length(naive),2);
cellDepthsNaiveAll = cell(length(naive),1);
cellTunedAngleNaiveAll = cell(length(naive),2); %(:,1) from spikes, (:,2) from whisker-only model

for mi = 1 : length(naive)
    tempInd = find(naive(mi).tuned);
    cellDepthsNaiveAll{mi} = naive(mi).depth(tempInd);
    tempSpikeAll = cell2mat(naive(mi).atSpikeAll(tempInd));
    [~, tempAngleInd] = max(tempSpikeAll, [], 2);
    cellTunedAngleNaiveAll{mi,1} = 30 + tempAngleInd * 15;
    tempWhiskerAll = cell2mat(cellfun(@(x) x(1,:), naive(mi).atWhiskerAll(tempInd), 'un', 0));
    [~, tempAngleInd] = max(tempWhiskerAll, [], 2);
    cellTunedAngleNaiveAll{mi,2} = 30 + tempAngleInd * 15;
    
    tempMat = cell2mat(naive(mi).atCorrFull(tempInd));
    [sortedFeatureImportanceNaiveAll{mi,1}, sortedFeatureIndNaive{mi,1}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
    tempMat = cell2mat(naive(mi).atCorrWhisker(tempInd));
    [sortedFeatureImportanceNaiveAll{mi,2}, sortedFeatureIndNaive{mi,2}] = sort(tempMat(:,1) - tempMat(:,2:end),2, 'descend');
end
%
% distribution 1st to 3rd feature importance
% and L2/3prop
histRange = 0:0.02:1;
distFullNaiveAll = zeros(length(naive),length(histRange)-1, 3);
distWhiskerNaiveAll = zeros(length(naive),length(histRange)-1, 3);
L23propNaiveAll = zeros(length(naive),2); % (:,1) full model, (:,2) whisker-only model
L4propNaiveAll = zeros(length(naive),2);
for mi = 1 : length(naive)
    for i = 1 : 3
        distFullNaiveAll(mi,:,i) = histcounts(sortedFeatureImportanceNaiveAll{mi,1}(:,i), histRange, 'norm', 'cdf');
        distWhiskerNaiveAll(mi,:,i) = histcounts(sortedFeatureImportanceNaiveAll{mi,2}(:,i), histRange, 'norm', 'cdf');
    end
    for i = 1 : 2
        multiInd = find(sortedFeatureImportanceNaiveAll{mi,i}(:,2) > atCorrThreshold);
        L23propNaiveAll(mi,i) = length(find(cellDepthsNaiveAll{mi}(multiInd) <= depthThreshold)) / length(find(cellDepthsNaiveAll{mi} <= depthThreshold));
        L4propNaiveAll(mi,i) = length(find(cellDepthsNaiveAll{mi}(multiInd) > depthThreshold)) / length(find(cellDepthsNaiveAll{mi} > depthThreshold));
    end
end

% Draw
barOffset = 0.2;
barWidth = 0.4;
figure('units','inch','pos',[5 5 3 4]), hold on
bar(1-barOffset, mean(L23propNaiveAll(:,2)), barWidth, 'facecolor', 'k')
bar(1+barOffset, mean(L4propNaiveAll(:,2)), barWidth, 'facecolor', 'w')
errorbar(1-barOffset, mean(L23propNaiveAll(:,2)), sem(L23propNaiveAll(:,2)), 'k')
errorbar(1+barOffset, mean(L4propNaiveAll(:,2)), sem(L4propNaiveAll(:,2)), 'k')
[~,p,m] = paired_test(L23propNaiveAll(:,2), L4propNaiveAll(:,2));
title({'Multi-whisker feature neurons, Naive (n = 12)'; sprintf('p = %.4f m = %s', p, m)})
ylabel('Proportion')
xlim([0.5 1.5]), yticks([0:0.05:0.2])
xticks([])

