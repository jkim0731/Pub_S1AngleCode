

baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% basic settings

calciumDir = [baseDir, 'Calcium\'];
whiskerModelDir = [baseDir, 'Whisker\Models\'];
colors = [248 171 66; 40 170 225] / 255;
featureNames = {'Push angle', 'Vertical displacement', 'Horizontal bending', 'Vertical bending', 'Slide distance', 'Touch duration', ...
    'Horizontal angle', 'Vertical angle', 'Horizontal curvature', 'Vertical curvature', 'Arc length', 'Touch count'};


%%
%% Figure 8A - Training with 2 angles
%%

%%
%% Figure 8B - Whisker feature distribution depending on 2 object angles
%%
%% Naive (top row)
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Naive'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

angles = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features{mi,ai,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

figure,

tempColorList = jet(7);
colorList = tempColorList([1,7],:);
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(angles));
    for ai = 1 : length(angles)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ai) = histcounts(tempCellFeature{mi,ai}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ai));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList(ai,:)})
    end
    ylim([0 0.3])
    
    title(featureNames{i})
    set(gca, 'fontname', 'Arial')
end
legend({'45\circ', '135\circ'})
sgtitle('Naive')

%% Expert (bottom row)

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

angles = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features{mi,ai,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

figure,

tempColorList = jet(7);
colorList = tempColorList([1,7],:);
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(angles));
    for ai = 1 : length(angles)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ai) = histcounts(tempCellFeature{mi,ai}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ai));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList(ai,:)})
    end
    ylim([0 0.3])
    
    title(featureNames{i})
    set(gca, 'fontname', 'Arial')
end
legend({'45\circ', '135\circ'})
sgtitle('Expert')

%%
%% 8C - Whisker feature distribution depending on 2 choices
%%
%% Naive (top row)
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Naive'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

angles = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features{mi,ai,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

figure,
colorList = {'c','m'};
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(angles));
    for ai = 1 : length(angles)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ai) = histcounts(tempCellFeature{mi,ai}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ai));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList{ai}})
    end
    ylim([0 0.3])
    
    title(featureNames{i})
    set(gca, 'fontname', 'Arial')
end
legend({'45\circ', '135\circ'})
sgtitle('Naive')

%% Expert (bottom row)

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

angles = sort(unique(groupMdl{1}.io.Y), 'descend');
features = cell(length(groupMdl), length(angles), size(groupMdl{1}.io.X,2));
for mi = 1 : length(groupMdl)
    for ai = 1 : length(angles)
        tempInd = find(groupMdl{mi}.io.Y == angles(ai));
        for i = 1 : size(features,3)
            features{mi,ai,i} = groupMdl{mi}.io.X(tempInd,i);
        end
    end
end

figure,
colorList = {'c','m'};
index = [1,2,5,6,9,10,3,4,7,8,11,12];
for i = 1 : size(features,3)
    subplot(3,4,index(i)), hold on
    tempCellFeature = features(:,:,i);
    histRange = -3:0.2:3.1; % standardized values
    hist = zeros(length(groupMdl), length(histRange)-1, length(angles));
    for ai = 1 : length(angles)
        for mi = 1 : length(groupMdl)            
            hist(mi,:,ai) = histcounts(tempCellFeature{mi,ai}, histRange, 'normalization', 'probability');
        end
        tempMat = squeeze(hist(:,:,ai));
        
        shadedErrorBar(histRange(1:end-1), mean(tempMat), std(tempMat)/sqrt(length(groupMdl)), 'lineprop', {'color', colorList{ai}})
    end
    ylim([0 0.3])
    
    title(featureNames{i})
    set(gca, 'fontname', 'Arial')
end
legend({'45\circ', '135\circ'})
sgtitle('Expert')




%%
%% 8D - Decoders performances
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
timing = 'lick'; % 'lick' or 'answer'


learnedGroup = {'Naive', 'Expert'}; % 'Naive', 'Expert', or ''
taskGroup = {'Two', 'Two'}; % 'Two', 'Discrete', or 'RadialDistance'

choicePerformance = zeros(6,2);
Yout = 'Choice'; % 'Touch' or 'Choice'
for i = 1 : length(learnedGroup)
    fn = ['mdl', taskGroup{i}, learnedGroup{i}, Xhow, Yout, '_new_', timing];
    data = load([whiskerModelDir, fn]);
    choicePerformance(:,i) = cellfun(@(x) mean(x.gof.mcc), data.groupMdl)';
end

anglePerformance = zeros(6,2);
Yout = 'Touch'; % 'Touch' or 'Choice'
for i = 1 : length(learnedGroup)
    fn = ['mdl', taskGroup{i}, learnedGroup{i}, Xhow, Yout, '_new_', timing];
    data = load([whiskerModelDir, fn]);
    anglePerformance(:,i) = cellfun(@(x) mean(x.gof.mcc), data.groupMdl)';
end

figure
hold on
plot(mean(anglePerformance), 'ko--')
plot(mean(choicePerformance), 'ko-', 'markerfacecolor', 'k')
errorbar(mean(anglePerformance), std(anglePerformance)/sqrt(6), 'k.')
errorbar(mean(choicePerformance), std(choicePerformance)/sqrt(6), 'k.')

ylabel('Performance (MCC)')
legend({'Angle decoder', 'Choice decoder'}, 'location', 'southeast')

xlim([0.5 2.5])
xticks(1:2)
yticks(0:0.2:1)
xticklabels({'Naive', 'Expert'})
set(gca, 'fontsize', 12, 'fontName', 'Arial')


%%
%% 8E - Feature importance for choice prediction (Expert 2-angle)
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Two'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

nGroup = length(groupMdl);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = mean(groupMdl{gi}.fitCoeffs,2);
    
    predMat = dataX * coeffs;    
    predY = exp(predMat)./(1+exp(predMat));
    LL(gi,1) = sum(log(binopdf(dataY,1,predY)));
    nuLL = sum(log(binopdf(dataY,1,0.5)));
    satLL = 0;
    DE(gi,1) = (LL(gi,1)-nuLL)/(satLL-nuLL);
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat)./(1+exp(predMat));
        
        LL(gi,fi+1) = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,fi+1) = (LL(gi,fi+1)-nuLL)/(satLL-nuLL);
    end
end

orderInd = [1,8:13, 2:5,7,6];

figure, 
tempMat = (DE(:,1) - DE(:,2:end)) ./ DE(:,1);
tempMat = tempMat(:,orderInd(2:length(orderInd))-1);
bar([1:6, 8:size(DE,2)], mean(tempMat), 'k'), hold on
errorbar([1:6, 8:size(DE,2)], mean(tempMat), std(tempMat)/sqrt(nGroup), 'k', 'linestyle', 'none')
xticklabels({featureNames{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylimVal = ylim();
ylim([0 ylimVal(2)])
ylabel('Feature importance for choice prediction')
box('off')
set(gca, 'fontsize', 12, 'fontname', 'Arial')








%%
%% S8A - Feature importance for choice prediction (Expert 7-angle)
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Choice'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or ''
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
load([whiskerModelDir, fn])

nGroup = length(groupMdl);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = mean(groupMdl{gi}.fitCoeffs,2);
    
    predMat = dataX * coeffs;    
    predY = exp(predMat)./(1+exp(predMat));
    LL(gi,1) = sum(log(binopdf(dataY,1,predY)));
    nuLL = sum(log(binopdf(dataY,1,0.5)));
    satLL = 0;
    DE(gi,1) = (LL(gi,1)-nuLL)/(satLL-nuLL);
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat)./(1+exp(predMat));
        
        LL(gi,fi+1) = sum(log(binopdf(dataY,1,predY)));
        nuLL = sum(log(binopdf(dataY,1,0.5)));
        satLL = 0;
        DE(gi,fi+1) = (LL(gi,fi+1)-nuLL)/(satLL-nuLL);
    end
end

orderInd = [1,8:13, 2:5,7,6];

figure, 
tempMat = (DE(:,1) - DE(:,2:end)) ./ DE(:,1);
tempMat = tempMat(:,orderInd(2:length(orderInd))-1);
bar([1:6, 8:size(DE,2)], mean(tempMat), 'k'), hold on
errorbar([1:6, 8:size(DE,2)], mean(tempMat), std(tempMat)/sqrt(nGroup), 'k', 'linestyle', 'none')
xticklabels({featureNames{orderInd(2:length(orderInd))-1}})
xtickangle(45)
ylimVal = ylim();
ylim([0 ylimVal(2)])
ylabel('Feature importance for choice prediction')
box('off')
set(gca, 'fontsize', 12, 'fontname', 'Arial')







%%
%% Figure 8F - Testing angle tuning in 7 angle sessions, before and after learning
%%




%%
%% 8G
%%

load([calciumDir, 'cellID_match_persAT_v9'], 'match', 'cellIDpersATnaive', 'cellIDpersATexpert')
load([calciumDir, 'wkv_angle_tuning_model_v9'], 'expert', 'naive')

colorsTransient = [248 171 66; 40 170 225] / 255;
colorsPersistent = [1 0 0; 0 0 1];
angles = 45:15:135;

learnerInd = [1:4,7,9];
naiveLearner = naive(learnerInd);
jk027inds = find(naiveLearner(2).cIDAll > 5000);
jk027AllLength = length(naiveLearner(2).cIDAll);

fn = fieldnames(naiveLearner);
for i = 1 : length(fn)
    if size(naiveLearner(2).(fn{i}),1) == jk027AllLength       
        naiveLearner(2).(fn{i}) = naiveLearner(2).(fn{i})(jk027inds,:);
    elseif size(naiveLearner(2).(fn{i}),1) == jk027AllLength
        naiveLearner(2).(fn{i}) = naiveLearner(2).(fn{i})(:,jk027inds);
    end
end

numMice = length(learnerInd);

%%
%% Fig 8G
%%
%% change of feature encoding (angle-tuned neurons)
%% pooled neurons

encodingPop = cell(numMice,2);
for mi = 1 : numMice
    % naive 
    indTuned = find(naiveLearner(mi).tuned);
    encodingPop{mi,1} = naiveLearner(mi).deFull(indTuned) - naiveLearner(mi).deWhiskerFeatures(indTuned,:);
    
    % expert
    indTuned = find(expert(mi).tuned);
    encodingPop{mi,2} = expert(mi).deFull(indTuned) - expert(mi).deWhiskerFeatures(indTuned,:);
    
end


wfEncodingNaive = cell2mat(encodingPop(:,1));
wfEncodingExpert = cell2mat(encodingPop(:,2));
wfEncodingExpert(isinf(wfEncodingExpert)) = nan;
barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, nanmean(wfEncodingNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, nanmean(wfEncodingExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, nanmean(wfEncodingNaive(:,1:6)), sem(wfEncodingNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, nanmean(wfEncodingExpert(:,1:6)), sem(wfEncodingExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
yticks(0:0.01:0.03)
ylabel({'Feature importance on fitting'; 'session-wide neuronal activity'})
title('Angle-tuned neurons')

p = nan(1,6);
m = cell(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(wfEncodingNaive(:,i), wfEncodingExpert(:,i));
end



%%
%% S8B - Show animal-averaged
%%
wfEncodingNaive = cell2mat(cellfun(@mean, encodingPop(:,1), 'un', 0));
wfEncodingExpert = cell2mat(cellfun(@(x) mean(x(isfinite(sum(x,2)),:)), encodingPop(:,2), 'un', 0));
barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, mean(wfEncodingNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, mean(wfEncodingExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, mean(wfEncodingNaive(:,1:6)), sem(wfEncodingNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, mean(wfEncodingExpert(:,1:6)), sem(wfEncodingExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
ylabel({'Mean feature importance';  'on fitting session-wide neuronal activity'})
title('Angle-tuned neurons')
xlim([0.5 6.5])
yticks(0:0.01:0.03)

p = nan(1,6);
m = cell(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(wfEncodingNaive(:,i), wfEncodingExpert(:,i));
end





%%
%% 8H - Feature impact on tuning
%%
% pooled neurons
featImpact = cell(numMice, 2);
for mi = 1 : numMice
    % naive 
    indTuned = find(naiveLearner(mi).tuned);
    featImpact{mi,1} = cell2mat(cellfun(@(x) x(1) - x(2:end), naiveLearner(mi).atCorrWhisker(indTuned), 'un', 0));
    
    % expert
    indTuned = find(expert(mi).tuned);
    featImpact{mi,2} = cell2mat(cellfun(@(x) x(1) - x(2:end), expert(mi).atCorrWhisker(indTuned), 'un', 0));
end

fiNaive = cell2mat(featImpact(:,1));
fiExpert = cell2mat(featImpact(:,2));
fiExpert(isinf(fiExpert)) = nan;
barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, nanmean(fiNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, nanmean(fiExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, nanmean(fiNaive(:,1:6)), sem(fiNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, nanmean(fiExpert(:,1:6)), sem(fiExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
yticks(-0.05:0.05:0.15)
ylabel({'Feature importance on fitting'; 'session-wide neuronal activity'})
title('Angle-tuned neurons')

p = nan(1,6);
m = cell(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(fiNaive(:,i), fiExpert(:,i));
end



%%
%% S8C - Mean feature impact on tuning
%%
fiNaive = cell2mat(cellfun(@(x) mean(x(isfinite(sum(x,2)),:)), featImpact(:,1), 'un', 0));
fiExpert = cell2mat(cellfun(@(x) mean(x(isfinite(sum(x,2)),:)), featImpact(:,2), 'un', 0));
barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, mean(fiNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, mean(fiExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, mean(fiNaive(:,1:6)), sem(fiNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, mean(fiExpert(:,1:6)), sem(fiExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
ylabel({'Mean feature importance';  'on fitting session-wide neuronal activity'})
title('Angle-tuned neurons')
xlim([0.5 6.5])
yticks(-0.05:0.05:0.15)

p = nan(1,6);
m = cell(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(fiNaive(:,i), fiExpert(:,i));
end




%%
%% S8D - Proportion of maximally impacting whisker features
%%
maxImpact = cell(numMice, 2);
for mi = 1 : numMice
    % naive
    indTuned = find(naiveLearner(mi).tuned);
    sigImpactFeat = cell2mat(cellfun(@(x) x(1)-x(2:end) == max(x(1)-x(2:end)), naiveLearner(mi).atCorrWhisker(indTuned), 'un', 0));
    maxImpact{mi,1} = mean(sigImpactFeat);
    
    % expert
    indTuned = find(expert(mi).tuned);
    sigImpactFeat = cell2mat(cellfun(@(x) x(1)-x(2:end) == max(x(1)-x(2:end)), expert(mi).atCorrWhisker(indTuned), 'un', 0));
    maxImpact{mi,2} = mean(sigImpactFeat);
end
miNaive = cell2mat(maxImpact(:,1));
miExpert = cell2mat(maxImpact(:,2));

barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, nanmean(miNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, nanmean(miExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, nanmean(miNaive(:,1:6)), sem(miNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, nanmean(miExpert(:,1:6)), sem(miExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
yticks(0:0.1:0.5)
ylabel('Proportion of angle-tuned neurons')
title('Maximally impacting object-angle tuning')

p = nan(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(miNaive(:,i), miExpert(:,i));
end



%%
%% S8E - Proportion of strongly impacting whisker features
%%
impactThreshold = 0.1;
sigImpactProp = cell(numMice, 2);
for mi = 1 : numMice
    % naive
    indTuned = find(naiveLearner(mi).tuned);
    sigImpactFeat = cell2mat(cellfun(@(x) x(1)-x(2:end) > impactThreshold, naiveLearner(mi).atCorrWhisker(indTuned), 'un', 0));
    sigImpactProp{mi,1} = mean(sigImpactFeat);
    
    % expert
    indTuned = find(expert(mi).tuned);
    sigImpactFeat = cell2mat(cellfun(@(x) x(1)-x(2:end) > impactThreshold, expert(mi).atCorrWhisker(indTuned), 'un', 0));
    sigImpactProp{mi,2} = mean(sigImpactFeat);
end
siNaive = cell2mat(sigImpactProp(:,1));
siExpert = cell2mat(sigImpactProp(:,2));

barOffset = 0.2;
barWidth = 0.4;
xpos = [1,2,3,4,6,5];
figure,
hold on
bar(xpos - barOffset, nanmean(siNaive(:,1:6)), barWidth, 'facecolor', colorsTransient(1,:))
bar(xpos + barOffset, nanmean(siExpert(:,1:6)), barWidth, 'facecolor', colorsTransient(2,:))
legend({'Naive', 'Expert'}, 'autoupdate', false, 'box', 'off', 'location', 'northwest')
errorbar(xpos - barOffset, nanmean(siNaive(:,1:6)), sem(siNaive(:,1:6)), 'k', 'lines', 'no')
errorbar(xpos + barOffset, nanmean(siExpert(:,1:6)), sem(siExpert(:,1:6)), 'k', 'lines', 'no')
xticks(min(xpos):max(xpos))
xticklabels(featureNames(xpos))
xtickangle(45)
yticks(0:0.1:0.4)
ylabel('Proportion of angle-tuned neurons')
title('Strongly impacting object-angle tuning')

p = nan(1,6);
for i = 1 : 6
    [~, p(xpos(i))] = ttest2(siNaive(:,i), siExpert(:,i));
end


