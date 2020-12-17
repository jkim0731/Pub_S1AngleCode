
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% Load calcium data

calciumDir = [baseDir, 'Calcium\'];
% mouse = 25;
% session = 4;
mouse = 39;
session = 1;
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

naiveInd = find(mice == mouse);
learnerInd = [1:4,7,9];
nonlearnerInd = [5,6,8,10:12];

repeat = 10;
% %% dependent settings
ufn = sprintf('UberJK%03dS%02d_NC',mouse, session);
glmfnBase = sprintf('glmResponseType_JK%03dS%02d_lasso_NC_R', mouse, session);
angletuningFn = sprintf('JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC_permTestCorrected', mouse, session);
% %% load files
load([calciumDir, 'glmResults_devExp_touch_NC'], 'naive', 'expert')

naiveLearner = naive(learnerInd);
jk027inds = find(naiveLearner(2).cellID > 5000);
jk027AllLength = length(naiveLearner(2).cellID);

fn = fieldnames(naiveLearner);
for i = 1 : length(fn)
    if length(naiveLearner(2).(fn{i})) == jk027AllLength
        naiveLearner(2).(fn{i}) = naiveLearner(2).(fn{i})(jk027inds);
    end
end

load(sprintf('%s%03d\\%s',calciumDir, mouse,ufn), 'u') % loading u
load(sprintf('%s%03d\\%s01',calciumDir, mouse,glmfnBase), 'allPredictors', 'posShift');
load(sprintf('%s%03d\\%s',calciumDir, mouse,angletuningFn))

%%
%% Fig 3A
%%
%% figure of touch GLM
cID = 5087;

ci = find(u.cellNums == cID);
traces = get_traces_per_cell(u, cID, allPredictors, naive(naiveInd).coeffs{ci}, posShift);
figure('units', 'normalized', 'outerposition', [0.3 0.3 0.4 0.4]), hold on
plot(min_max_normalization(traces.calcium)*4 + 8 , 'color', [0.1 0.8 0.1], 'linewidth', 2)
plot(min_max_normalization(traces.spikes)*2.5 + 6.2, 'color', [0.8 0.1 0.1], 'linewidth', 2)

plot(min_max_normalization(traces.model)*2.5 + 6.2, 'k-', 'linewidth', 2)

angleColor = jet(7);

plot(min_max_normalization(traces.predictors(:,8)) + 3.2+0.8, 'color', [0.7 0.7 0.7], 'linewidth', 2)
for i = 1:7
    plot(min_max_normalization(traces.predictors(:,i)) + 3.2+0.8-0.1*i, 'color', angleColor(i,:), 'linewidth', 2)
end
plot(min_max_normalization(traces.predictors(:,41)) + 2.2, 'color', [0.7 0.7 0.7], 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,48)) + 1.5, 'color', [0.7 0.7 0.7], 'linewidth', 2)
plot(min_max_normalization(traces.predictors(:,55)) + 0.8, 'color', [0.7 0.7 0.7], 'linewidth', 2)
lickL = min_max_normalization(traces.predictors(:,62));
lickR = min_max_normalization(traces.predictors(:,66));
plot(lickL, 'color', 'm', 'linewidth', 2)
plot(lickR, 'color', 'c', 'linewidth', 2)
rewardL = min_max_normalization(traces.predictors(:,29));
rewardL(rewardL<max(rewardL)) = 0;
rewardL(isnan(rewardL)) = 0;
rewardR = min_max_normalization(traces.predictors(:,34));
rewardR(rewardR<max(rewardR)) = 0;
rewardR(isnan(rewardR)) = 0;
reward = rewardL + rewardR;
rewarded = find(reward);
rdiffOne = find(diff(find(reward))==1);
reward(rewarded(rdiffOne+1)) = nan;
reward(reward==0) = nan;
sound = min_max_normalization(traces.predictors(:,25));
sound(sound<max(sound)) = nan;
plot(reward-1.5, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'b')
plot(sound-1.5, '^', 'markersize', 4, 'markeredgecolor', 'none', 'markerfacecolor', 'k')
onedF = 4 / max(traces.calcium);
oneSpk = 2.5/max(traces.spikes);
plot([930 930], [10 10+onedF], 'color', [0.1 0.8 0.1], 'linewidth', 2) % 1 dF/F0 scale bar
plot([900 900 + u.frameRate*5], [10 10], 'k-', 'linewidth', 2) % 5 s scale bar
plot([930 930], [10 10 - oneSpk], 'color', [0.8 0.1 0.1], 'linewidth', 2) % 1 spk scale bar
xlim([240 940]), ylim([-0.5 12])

axis off
set(gcf, 'InvertHardCopy', 'off', 'color', 'w');
DE = naive(naiveInd).allDE(ci)


%%
%% S3A - whisking parameters
%%
mouse = 25;
session = 4;
load(sprintf('%s%03d\\UberJK%03dS%02d_NC', calciumDir, mouse, mouse, session), 'u')
%%
[~,uind] = sort(cellfun(@(x) nanstd(x.theta), u.trials), 'descend');
ii = 50;
utrial = u.trials{uind(ii)};
colorList = [0 0 0; 0 0.2 0.8; 0.7 0.4 0];
figure, 
subplot(211), hold on
plot(utrial.whiskerTime, utrial.theta, '-', 'color', colorList(1,:))
[onsetFrames, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(utrial.theta);
plot(utrial.whiskerTime, amplitude, '-', 'color', colorList(2,:)), 
plot(utrial.whiskerTime, midpoint, '-', 'color', colorList(3,:))
plot(utrial.whiskerTime, 2.5 * ones(length(utrial.whiskerTime),1), '--', 'color', [0.5 0.5 0.5])
set(gca, 'box', 'off', 'fontsize', 14)
xticklabels([])
xlim([0 4.5])

tpmWhiskerFrames = cell(length(utrial.tpmTime{1}),1);
tpmTime = [0, utrial.tpmTime{1}];
for i = 1 : length(tpmWhiskerFrames)
    tpmWhiskerFrames{i} = find(utrial.whiskerTime > tpmTime(i) & utrial.whiskerTime <= tpmTime(i+1));
end
onsetWhisker = zeros(length(utrial.whiskerTime),1);
onsetWhisker(onsetFrames) = 1;
onsetTPM = cellfun(@(x) sum(onsetWhisker(x)), tpmWhiskerFrames);
ampTPM = cellfun(@(x) mean(amplitude(x)), tpmWhiskerFrames);
midTPM = cellfun(@(x) mean(midpoint(x)), tpmWhiskerFrames);

subplot(212), hold on
yyaxis right
plot(utrial.tpmTime{1}, onsetTPM, '-', 'color', colorList(1,:))
yyaxis left
plot(utrial.tpmTime{1}, ampTPM, '-', 'color', colorList(2,:))
plot(utrial.tpmTime{1}, midTPM, '-', 'color', colorList(3,:))
set(gca, 'box', 'off', 'fontsize', 14)
xlim([0 4.5])
xlabel('Time (s)')

%%
%% S3B - % deviance explained distribution. 
%%
%% From all naive mice
load([calciumDir, 'glmResults_devExp_touch_NC.mat'])
histRange = [-10,0:0.02:0.62,10];
distDE = zeros(length(naive), length(histRange)-1);
for i = 1 : length(naive)
    distDE(i,:) = histcounts(naive(i).allDE, histRange, 'normalization', 'cdf');
end

figure, hold on
boundedline(histRange(2:end-1),mean(distDE(:,2:end)), sem(distDE(:,2:end)), 'k')
xlim([0 0.6])
ylim([0.4 1])
plot([0.1 0.1], [0.4 1], '--', 'color', ones(1,3)*0.7)
xlabel('Goodness-of-fit')
ylabel('Cumulative proportion')

%%
propFit = zeros(length(naive), 1);
for i = 1 : length(naive)
    propFit(i) = length(find(naive(i).allDE > 0.1)) / length(naive(i).allDE);
end

mean(propFit)
sem(propFit)





%%
%% 3B - Proportion of neuronal response types
%%
info = load([calciumDir, 'cellFunctionLasso_NC.mat']);
cellFuncProp = zeros(length(info.naive),4); % 1: touch, 2: touch & whisking, 3: whisking, 4: other
for ni = 1 : length(info.naive)
    numFit = length(info.naive(ni).allCell.cellNums);
    tempMat = cell2mat(info.naive(ni).cellFunction);
    touchTemp = sum(tempMat(:,1));
    whiskingTemp = sum(tempMat(:,4));
    mixed = length(intersect(find(tempMat(:,1)), find(tempMat(:,4))));
    touch = touchTemp - mixed;
    whisking = whiskingTemp - mixed;
    other = sum(sum(tempMat(:,[2,3,5])));
    cellFuncProp(ni,1) = touch / numFit;
    cellFuncProp(ni,2) = mixed / numFit;
    cellFuncProp(ni,3) = whisking / numFit;
    cellFuncProp(ni,4) = other / numFit;
end

figure,
pie(flip([mean(cellFuncProp) 1-sum(mean(cellFuncProp))]))
mean(cellFuncProp)
sem(cellFuncProp)

mean(1-sum(cellFuncProp,2))
sem(1-sum(cellFuncProp,2))






%%
%% Fig 3C
%%
%% example from JK039 S01
mouse = 39;
session = 1;
mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]};

naiveInd = find(mice == mouse);

learnerInd = [1:4,7,9];
nonlearnerInd = [5,6,8,10:12];

load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_predecision_NC', calciumDir, mouse, mouse, session), 'ca')
load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC_permTestCorrected', calciumDir, mouse, mouse, session), 'spk', 'info')
load(sprintf('%s%03d\\UberJK%03dS%02d_NC', calciumDir, mouse, mouse, session))

touchGlmNaive = load([calciumDir, 'glmResults_devExp_touch_NC.mat'], 'naive');
% whiskerGlmNaive = load([baseDir, 'glmResults_devExp_WKV_touchCell_NC.mat'], 'naive');
% tuneNaive = load([calciumDir, 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat'], 'naive');

touchGLM = touchGlmNaive.naive(naiveInd);
% whiskerGLM = whiskerGlmNaive.naive(naiveInd);
% tune = tuneNaive.naive(naiveInd);

%% following cell shape
plane = 5;

touchID = spk.touchID(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);
tunedAngle = spk.tunedAngle(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);
unimodalSingle = spk.unimodalSingle(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);
unimodalBroad = spk.unimodalBroad(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);
multimodal = spk.multimodal(spk.touchID > plane*1000 & spk.touchID < (plane+1) * 1000);

mimg = (mat2gray(u.mimg{plane}));
h = figure;
imshow((mimg))
hold on,

plot([u.c2xpoints, u.c2xpoints(1)], [u.c2ypoints, u.c2ypoints(1)], 'w--', 'linewidth', 2)

colors = jet(7);
cIDlist = spk.touchID(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
tunedAngles = spk.tunedAngle(find(spk.touchID > plane * 1000 & spk.touchID < (plane+1) * 1000));
indList = find(ismember(u.cellNums, cIDlist));

cellmap = u.cellmap{plane};
for ci = 1 : length(cIDlist)
    cID = cIDlist(ci);
    [ypix,xpix] = ind2sub(size(cellmap),find(cellmap == cID));    
    k = boundary(ypix, xpix);    
    if tunedAngles(ci) == 0
        patch(xpix(k), ypix(k), 'k', 'edgecolor', 'k')
    else
        angleInd = (tunedAngles(ci)-30)/15;
        patch(xpix(k), ypix(k), colors(angleInd,:), 'edgecolor', 'none')
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
%% 3C(i)-(iv) & S3C(i)-(iv) Individual examples
%%
%% pick one from tunedAngle

cid = 5094; % specific
% cid = 5087; % broad
% cid = 5367; % complex
% cid = 5020; % non-selective
example_angle_tuning_calcium(u, ca, spk, cid)

cind = find(u.cellNums==cid);
figure(1113), 
imshow(mat2gray(u.mimg{plane}))
hold on
plot(u.cellx(cind), u.celly(cind), 'y.','markersize', 20)




%%
%% Fig 3D - tiling all angle (from all naive)
%%

tune = load([calciumDir,'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat']);

stretchFactor = 300;
sharpness = [];
tunedAngle = [];
val = {};
for i = 1 : length(tune.naive)
    indTuned = find(tune.naive(i).tuned);
    sharpness = [sharpness; tune.naive(i).sharpness(indTuned)];
    tunedAngle = [tunedAngle; tune.naive(i).tunedAngle(indTuned)];
    val = [val; tune.naive(i).val(indTuned)];
end
[~, indsortSharpness] = sort(sharpness, 'descend');
tempTunedAngle = tunedAngle(indsortSharpness);
[~, indsortAngle] = sort(tempTunedAngle);
indsort = indsortSharpness(indsortAngle);

tilingMap = zeros(length(tunedAngle), 7 * stretchFactor);
for i = 1 : length(tunedAngle)
    tempMap = cellfun(@mean, val{indsort(i)});
    normTempMap = min_max_normalization(tempMap);
    for j = 1 : 7
        tilingMap(i,(j-1)*stretchFactor+1:j*stretchFactor) = deal(normTempMap(j));
    end
end
tempMap = tilingMap(:,[1,stretchFactor*1+1, stretchFactor*2+1, stretchFactor*3+1, stretchFactor*4+1, stretchFactor*5+1, stretchFactor*6+1]);
inds = [1:size(tempMap,1)]';
numGroup = size(tempMap,2);
sortiCell = cell(numGroup,1);
for i = 1 : numGroup
    row2sort = find(tempMap(:,i)==1); % since it's min-max normalized.
    if isempty(row2sort)
        sortiCell{i} = [];
    else
        col2sort = setdiff(1:numGroup,i);
        tempsorti = nested_sorting(tempMap(row2sort,col2sort),inds(row2sort));
        if size(tempsorti,2) ~= 1
            tempsorti = tempsorti';
        end
        sortiCell{i} = tempsorti;
    end
end
sorti = cell2mat(sortiCell);
   
newMap = tilingMap(sorti,:);

figure
imshow(newMap), hold on, 
axis tight



%%
%% 3E - tuned angle proportion naive all
%%

angles = 45:15:135;
tuning = zeros(length(tune.naive), length(angles));
for i = 1 : length(tune.naive)
    temp = tune.naive(i);
    for j = 1 : length(angles)
        tuning(i,j) = length(find(temp.tunedAngle == angles(j))) / sum(temp.tuned);
    end
end

figure,
errorbar(angles, mean(tuning), sem(tuning), 'k-')
xlabel('Object angle (\circ)')
xticks([45:15:135]), ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Proportion (/all tuned cells)')
set(gca, 'fontsize', 14, 'box', 'off')






%%
%% Fig 3F - across layers (tuning type distribution, preferred angle distribution)
%%

colors = [62, 38, 168; 49, 198, 159; 252, 207, 48; 153, 153, 153] / 255;
tuning = load(sprintf('%sangle_tuning_summary_preAnswer_perTouch_NC_PTC.mat',calciumDir),'naive', 'expert');

numMice = length(tuning.naive);
angles = 45:15:135;
tuningTypesL23 = zeros(numMice, 4); % (:,1) sharp, (:,2) broad, (:,3) complex, (:,4) non-selective
tuningTypesL4 = zeros(numMice, 4);
tunedAnglesL23 = zeros(numMice, 7);
tunedAnglesL4 = zeros(numMice, 7);
for mi = 1 : numMice
    indL23 = find(tuning.naive(mi).depth < 350);
    indL4 = find(tuning.naive(mi).depth >= 350);
    
    indSharp = (find(tuning.naive(mi).unimodalSingle));
    indBroad = (find(tuning.naive(mi).unimodalBroad));
    indComplex = (setdiff(find(tuning.naive(mi).multimodal), find(tuning.naive(mi).unimodalBroad)));
    indNS = (find(tuning.naive(mi).tuned == 0));
    
    tuningTypesL23(mi,1) = length(intersect(indSharp, indL23)) / length(indL23);
    tuningTypesL23(mi,2) = length(intersect(indBroad, indL23)) / length(indL23);
    tuningTypesL23(mi,3) = length(intersect(indComplex, indL23)) / length(indL23);
    tuningTypesL23(mi,4) = length(intersect(indNS, indL23)) / length(indL23);
    
    tuningTypesL4(mi,1) = length(intersect(indSharp, indL4)) / length(indL4);
    tuningTypesL4(mi,2) = length(intersect(indBroad, indL4)) / length(indL4);
    tuningTypesL4(mi,3) = length(intersect(indComplex, indL4)) / length(indL4);
    tuningTypesL4(mi,4) = length(intersect(indNS, indL4)) / length(indL4);
    
    indTuned = find(tuning.naive(mi).tuned);
    for ai = 1 : length(angles)
        indTunedAngle = find(tuning.naive(mi).tunedAngle == angles(ai));
        tunedAnglesL23(mi,ai) = length(intersect(indTunedAngle, indL23)) / length(intersect(indTuned, indL23));
        tunedAnglesL4(mi,ai) = length(intersect(indTunedAngle, indL4)) / length(intersect(indTuned, indL4));
    end
    
end

%% show tuning type distribution across layers
pbWidth= 0.39;
pbOffset = 0.2;
figure, hold on
for i = 1 : 4
    bar(i-pbOffset, mean(tuningTypesL23(:,i)), pbWidth, 'facecolor', colors(i,:))
    bar(i+pbOffset, mean(tuningTypesL4(:,i)), pbWidth, 'facecolor', 'w', 'edgecolor', colors(i,:))
    if i == 1
        legend({'L2/3', 'L4'}, 'autoupdate', false)
    end
end
errorbar([1:4]-pbOffset, mean(tuningTypesL23), sem(tuningTypesL23), 'k', 'lines', 'no')
errorbar([1:4]+pbOffset, mean(tuningTypesL4), sem(tuningTypesL4), 'k', 'lines', 'no')

xlim([0 5])
xticks(1:4)
xticklabels({'Sharp', 'Broad', 'Complex', 'Not-significant'})
xtickangle(45)
ylim([0 0.5])
yticks(0:0.1:0.5)
ylabel('Proportion')


[hTuning, pTuning, mTuning] = paired_test(tuningTypesL23, tuningTypesL4)



%%
%% 3G - Preferred angle distribution across layers
%%
figure, hold on
errorbar(angles, mean(tunedAnglesL23), sem(tunedAnglesL23), 'k-o', 'markerfacecolor','k')
errorbar(angles, mean(tunedAnglesL4), sem(tunedAnglesL4), 'k--o')

legend({'L2/3', 'L4'})
xlim([40 140]), xticks([45:15:135])
ylim([0 0.4]), yticks([0:0.1:0.4])

[hAngle, pAngle, mAngle] = paired_test(tunedAnglesL23, tunedAnglesL4)







%%
%% S3D - Angle-tuning from models (naive all)
%%

load([calciumDir, 'touchModel_angle_tuning_model_v9'])

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

%% naive all
atDist = zeros(length(naive),length(angles),3); % (:,:,1) spike, (:,:,2) full model, (:,:,3) touch-only model
for mi = 1 : length(naive)
    tempSpike = cellfun(@(x) find(x == max(x)), naive(mi).atSpikeAll(find(naive(mi).tunedSpk)));
    tempFull = cellfun(@(x) find(x == max(x)), naive(mi).atFullAll(find(naive(mi).tunedSpk)));
    tempTouch = cellfun(@(x) find(x == max(x)), naive(mi).atTouchAll(find(naive(mi).tunedSpk)));
    for ai = 1 : length(angles)
        atDist(mi,ai,1) = length(find(tempSpike == ai)) / length(tempSpike);
        atDist(mi,ai,2) = length(find(tempFull == ai)) / length(tempFull);
        atDist(mi,ai,3) = length(find(tempTouch == ai)) / length(tempTouch);
    end
end

colors = {'r','k','c'};
figure, hold on
for i = 1 : 3
    tempMat = squeeze(atDist(:,:,i));
    errorbar(angles, mean(tempMat), sem(tempMat), 'color', colors{i});
end
legend({'Inferred spike', 'Full model', 'Touch-only model'}, 'box', 'off', 'location', 'northwest')
xlim([40 140]), xticks(angles)
xlabel('Tuned angle (\circ)')
ylabel('Proportion')
title('Naive all (n = 12)')

%%
%% S3E - Tuned angle in tuned neurons from spikes (naive all)
%%

changeTunedAngle = cell(length(naive),2); % (:,1) spk to full, (:,2) spk to touch-only

histRange = [0:15:90,100];
distTAchangeFull = zeros(length(naive), length(histRange)-1);
distTAchangeTouch = zeros(length(naive), length(histRange)-1);
for mi = 1 : length(naive)
    ind = find(naive(mi).tunedSpk);
    [~, taSpk] = max(cell2mat(naive(mi).atSpikeAll(ind)), [], 2);
    [~, taFull] = max(cell2mat(naive(mi).atFullAll(ind)), [], 2);
    [~, taTouch] = max(cell2mat(naive(mi).atTouchAll(ind)), [], 2);
    
    changeTunedAngle{mi,1} = 15 * abs(taSpk - taFull);
    changeTunedAngle{mi,2} = 15 * abs(taSpk - taTouch);
    
    distTAchangeFull(mi,:) = histcounts(changeTunedAngle{mi,1}, histRange, 'norm', 'prob');
    distTAchangeTouch(mi,:) = histcounts(changeTunedAngle{mi,2}, histRange, 'norm', 'prob');
end

figure, hold on
plot(histRange(1:end-1), mean(distTAchangeFull), 'k')
plot(histRange(1:end-1), mean(distTAchangeTouch), 'c')
legend({'Full object model', 'Touch-only model'}, 'autoupdate', false)
boundedline(histRange(1:end-1), mean(distTAchangeFull), sem(distTAchangeFull), 'k')
boundedline(histRange(1:end-1), mean(distTAchangeTouch), sem(distTAchangeTouch), 'c')

xlabel('\Delta tuned angle (\circ)')
xticks(0:15:90)
ylabel('Proportion')
title('Naive all, tuned angles only (n = 12)')







%%
%% Scnn1a-Cre mice
%% 
% proportion of tuning types, angle selectivity across types, and
% distribution of tuned angles

mice = [75,76];
sessions = [4,4];
colors = [62, 38, 168; 49, 198, 159; 252, 207, 48; 153, 153, 153] / 255;

angles = 45:15:135;
selectivity = cell(length(mice),4); % (:,1) specific, (:,2) broad, (:,3) complex, (:,4) non-selective
propSelectivity = zeros(length(mice),4);
taDistribution = zeros(length(mice), length(angles));
for mi = 1 : length(mice)
    mouse = mice(mi);
    session = sessions(mi);
    load(sprintf('%s%03d\\JK%03dS%02dangle_tuning_lasso_preAnswer_perTouch_spkOnly_NC_permTestCorrected', calciumDir, mouse, mouse, session))
    typeInd = cell(1,4);
    typeInd{1} = find(spk.unimodalSingle);
    typeInd{2} = find(spk.unimodalBroad);
    typeInd{3} = setdiff(find(spk.multimodal), find(spk.unimodalBroad));
    typeInd{4} = find(spk.tuned == 0);
    for ti = 1 : 4
        selectivity{mi,ti} = cellfun(@(x) max(cellfun(@mean, x)) - (sum(cellfun(@mean, x)) - max(cellfun(@mean, x))) / 6, spk.val(typeInd{ti}));
        propSelectivity(mi,ti) = length(typeInd{ti}) / length(spk.touchID);
    end
    for ai = 1 : length(angles)
        taDistribution(mi,ai) = length(find(spk.tunedAngle == angles(ai))) / length(find(spk.tuned));
    end
end

%%
%% S3F - Proportion of tuning type
%% 
figure('units','inch','pos',[0 0 3 5]), hold on
for i = 1 : 4
    bar(i, mean(propSelectivity(:,i)), 'facecolor', 'w', 'edgecolor', colors(i,:))
    errorbar(i, mean(propSelectivity(:,i)), sem(propSelectivity(:,i)), 'k')
end
xlim([0 5]), xticks(1:4),
xticklabels({'Specific', 'Broad', 'Complex', 'Non-selective'})
xtickangle(45)
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Proportion')

%%
%% S3G - Tuned angle distribution
%%
figure,
errorbar(angles, mean(taDistribution), sem(taDistribution), 'k--')
xlim([40, 140]), xticks(angles)
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Proportion')
xlabel('Tuned angle')
