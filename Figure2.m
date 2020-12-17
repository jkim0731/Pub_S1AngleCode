
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')

%%
%% Figure 2A - 3D reconstruction of whisker, whisker pad, and the object

mouse = 'JK025';
session = 'S05';
sampleDir = [baseDir, 'Whisker\',mouse,session];
wfa = Whisker.WhiskerFinal_2padArray(sampleDir);
w3a = Whisker.Whisker3D_2padArray(sampleDir);

touchTrialInd = find(cellfun(@(x) length(x.protractionTFchunksByWhisking), wfa.trials));
targetAngle = 45;
angleTrialInd = find(cellfun(@(x) x.poleAngle == targetAngle, wfa.trials));
angleTouchInd = intersect(touchTrialInd, angleTrialInd);

[~,sortInd] = sort(cellfun(@(x) length(x.protractionTFchunksByWhisking{1}), wfa.trials(angleTouchInd)), 'descend');
angleTouchInd = angleTouchInd(sortInd);


%%
tIndInd = 2;
tInd = angleTouchInd(tIndInd);
w3Ind = find(w3a.trialNums == wfa.trialNums(tInd));

wf = wfa.trials{tInd};
w3 = w3a.trials{w3Ind};


allTouchFrames = cell2mat(cellfun(@(x) x', wf.protractionTFchunksByWhisking, 'uniformoutput', false));
allTouchTimes = wf.time(allTouchFrames);
allTouchFramesW3 = find(ismember(w3.time, allTouchTimes));
poleEdge = zeros(length(allTouchFramesW3),3);
for i = 1 : length(allTouchFramesW3)
    tempInd = allTouchFramesW3(i);
    [~,tempTouchInd] = min(sum((w3.trackerData{tempInd} - w3.intersectPoint(tempInd,:)).^2,2));
    poleEdge(i,:) = w3.trackerData{tempInd}(tempTouchInd,:);
end

polyX = polyfit(linspace(0,1, size(poleEdge,1)),poleEdge(:,1)', 1);
polyY = polyfit(linspace(0,1, size(poleEdge,1)),poleEdge(:,2)', 1);
polyZ = polyfit(linspace(0,1, size(poleEdge,1)),poleEdge(:,3)', 1);
numInitPoint = 2;
poleFitX = polyval(polyX,linspace(-2.5,2, numInitPoint));
poleFitY = polyval(polyY,linspace(-2.5,2, numInitPoint));
poleFitZ = polyval(polyZ,linspace(-2.5,2, numInitPoint));

circumPointNum = 100;
poleX = zeros(numInitPoint, circumPointNum);
poleY = zeros(numInitPoint, circumPointNum);
poleZ = zeros(numInitPoint, circumPointNum);

radius = 3; % the with of the pole is about 6 pixels, in case of JK025~041. 
for i = 1 : numInitPoint    
    poleX(i,:) = poleFitX(i);
    poleY(i,:) = poleFitY(i) + sin(linspace(0,2*pi,circumPointNum))*radius;
    poleZ(i,:) = poleFitZ(i) + cos(linspace(0,2*pi,circumPointNum))*radius;    
end


touchFramesWF = wf.protractionTFchunksByWhisking{1};
touchTimes = wf.time(touchFramesWF);
touchFramesW3 = find(ismember(w3.time, touchTimes));

wGradation = [0, linspace(0.8, 0.2, length(touchFramesW3)-2), 0];
pGradation = [0.7, linspace(0.6, 0.2, length(touchFramesW3)-2), 0];
whiskerColors = ones(length(touchFramesW3),3).*wGradation';

baseColors = [1-flip(pGradation'), zeros(length(touchFramesW3),2)];
touchColors = [zeros(length(touchFramesW3),2), 1-flip(pGradation')];

whiskerTracked = {w3.trackerData{touchFramesW3}};
touchPoints = zeros(length(touchFramesW3),3);
basePoints = w3.base(touchFramesW3,:);
mask3d = w3.mask3d;

anchorPoint = mean(basePoints);

isozpoints = length(find(mask3d(:,3) == mask3d(1,3)));
zlength = size(mask3d,1)/isozpoints;
maskx = zeros(isozpoints,zlength);
masky = zeros(isozpoints,zlength);
maskz = zeros(isozpoints,zlength);
for i = 1 : isozpoints
    maskx(i,:) = mask3d((i-1)*isozpoints+1:i*isozpoints,1);
    masky(i,:) = mask3d((i-1)*isozpoints+1:i*isozpoints,2);
    maskz(i,:) = mask3d((i-1)*isozpoints+1:i*isozpoints,3);
end
%
figure, hold on
subplot(211), hold on
surf((maskx-anchorPoint(1))/w3.pxPerMm, -(masky-anchorPoint(2))/w3.pxPerMm, (maskz-anchorPoint(3))/w3.pxPerMm, 'facecolor', 'k', 'linestyle', 'none', 'facealpha', 0.3)
for i = 1 : length(touchFramesW3)
    if i == 1 
        plot3((whiskerTracked{i}(:,1) - anchorPoint(1))/w3.pxPerMm, -(whiskerTracked{i}(:,2) - anchorPoint(2))/w3.pxPerMm, (whiskerTracked{i}(:,3) - anchorPoint(3))/w3.pxPerMm, ':', 'color', whiskerColors(i,:), 'linewidth', 2)
    elseif i == length(touchFramesW3)
        plot3((whiskerTracked{i}(:,1) - anchorPoint(1))/w3.pxPerMm, -(whiskerTracked{i}(:,2) - anchorPoint(2))/w3.pxPerMm, (whiskerTracked{i}(:,3) - anchorPoint(3))/w3.pxPerMm, '-', 'color', whiskerColors(i,:), 'linewidth', 2)
    else
        plot3((whiskerTracked{i}(:,1) - anchorPoint(1))/w3.pxPerMm, -(whiskerTracked{i}(:,2) - anchorPoint(2))/w3.pxPerMm, (whiskerTracked{i}(:,3) - anchorPoint(3))/w3.pxPerMm, '-', 'color', whiskerColors(i,:))
    end
    intersectPoint = w3.intersectPoint(touchFramesW3(i),:);
    [~,touchPointInd] = min(sum((whiskerTracked{i} - intersectPoint).^2,2));
    touchPoints(i,:) = whiskerTracked{i}(touchPointInd,:);
end

for i = 1 : length(touchFramesW3)
    plot3((touchPoints(i,1) - anchorPoint(1))/w3.pxPerMm, -(touchPoints(i,2) - anchorPoint(2))/w3.pxPerMm, (touchPoints(i,3) - anchorPoint(3))/w3.pxPerMm, '.', 'markersize', 20, 'color', touchColors(i,:))
    plot3((basePoints(i,1) - anchorPoint(1))/w3.pxPerMm, -(basePoints(i,2) - anchorPoint(2))/w3.pxPerMm, (basePoints(i,3) - anchorPoint(3))/w3.pxPerMm, '.', 'markersize', 20, 'color',baseColors(i,:))
end

poleX = (poleX - anchorPoint(1))/w3.pxPerMm;
poleY = (poleY - anchorPoint(2))/w3.pxPerMm;
poleZ = (poleZ - anchorPoint(3))/w3.pxPerMm;
surf(poleX,-poleY,poleZ + radius/w3.pxPerMm,'facecolor','k','LineStyle','none', 'facealpha', 0.3)
fill3(poleX(1,:), -poleY(1,:), poleZ(1,:) + radius/w3.pxPerMm, 'k', 'linestyle', 'none', 'facealpha', 0.3)
fill3(poleX(end,:), -poleY(end,:), poleZ(end,:) + radius/w3.pxPerMm, 'k', 'linestyle', 'none', 'facealpha', 0.3)
% 
xlim([min(min(poleX))-3 max(max(poleX))])
ypointsMin = min((cellfun(@(x) min(x(:,2)), whiskerTracked)) - anchorPoint(2))/w3.pxPerMm;
ypointsMax = max((cellfun(@(x) max(x(:,2)), whiskerTracked)) - anchorPoint(2))/w3.pxPerMm;
ylim([-ypointsMax -ypointsMin])
zlim([min(min(poleZ)) max(max(poleZ))+1])
xlabel('AP (mm)'), ylabel('ML (mm)'), zlabel('DV (mm)')
title([mouse, session, ' T#', num2str(w3.trialNum), ' ind ' num2str(tIndInd)])
axis equal

set(gca, 'fontsize', 14, 'fontname', 'Arial')

view([-55 -40])
xlim([min(min(poleX))-3 max(max(poleX))])


subplot(212), hold on

surf((maskx-anchorPoint(1))/w3.pxPerMm, -(masky-anchorPoint(2))/w3.pxPerMm, (maskz-anchorPoint(3))/w3.pxPerMm, 'facecolor', 'k', 'linestyle', 'none', 'facealpha', 0.3)
plot3((whiskerTracked{1}(:,1) - anchorPoint(1))/w3.pxPerMm, -(whiskerTracked{1}(:,2) - anchorPoint(2))/w3.pxPerMm, (whiskerTracked{1}(:,3) - anchorPoint(3))/w3.pxPerMm, ':', 'color', whiskerColors(1,:), 'linewidth', 2)
plot3((whiskerTracked{end}(:,1) - anchorPoint(1))/w3.pxPerMm, -(whiskerTracked{end}(:,2) - anchorPoint(2))/w3.pxPerMm, (whiskerTracked{end}(:,3) - anchorPoint(3))/w3.pxPerMm, '-', 'color', whiskerColors(end,:), 'linewidth', 2)
for i = 1 : length(touchFramesW3)
    plot3((touchPoints(i,1) - anchorPoint(1))/w3.pxPerMm, -(touchPoints(i,2) - anchorPoint(2))/w3.pxPerMm, (touchPoints(i,3) - anchorPoint(3))/w3.pxPerMm, '.', 'markersize', 20, 'color', touchColors(i,:))
end
plot3((basePoints(1,1) - anchorPoint(1))/w3.pxPerMm, -(basePoints(1,2) - anchorPoint(2))/w3.pxPerMm, (basePoints(1,3) - anchorPoint(3))/w3.pxPerMm, '.', 'markersize', 20, 'color',baseColors(1,:))
plot3((basePoints(end,1) - anchorPoint(1))/w3.pxPerMm, -(basePoints(end,2) - anchorPoint(2))/w3.pxPerMm, (basePoints(end,3) - anchorPoint(3))/w3.pxPerMm, '.', 'markersize', 20, 'color',baseColors(end,:))
surf(poleX,-poleY,poleZ + radius/w3.pxPerMm,'facecolor','k','LineStyle','none', 'facealpha', 0.3)
fill3(poleX(1,:), -poleY(1,:), poleZ(1,:) + radius/w3.pxPerMm, 'k', 'linestyle', 'none', 'facealpha', 0.3)
fill3(poleX(end,:), -poleY(end,:), poleZ(end,:) + radius/w3.pxPerMm, 'k', 'linestyle', 'none', 'facealpha', 0.3)
xlim([min(min(poleX))-3 max(max(poleX))])
ylim([-ypointsMax -ypointsMin])
zlim([min(min(poleZ)) max(max(poleZ))+1])
xlabel('AP (mm)'), ylabel('ML (mm)'), zlabel('DV (mm)')
axis equal

set(gca, 'fontsize', 14, 'fontname', 'Arial')
view([-55 -40])
xlim([min(min(poleX))-3 max(max(poleX))])

%%
%% 2B - Schematics to show 7-angle test
%%

%%
%% 2C - Feature discriminability for 7 angles
%% Compared between naive and expert
%% 

calciumBaseDir = [baseDir, 'Calcium\'];
mice = [25,27,30,36,39,52];
sessions = {[4,19], [3,10], [3,21], [1,17], [1,23], [3,21]};

colorsTransient = [248 171 66; 40 170 225] / 255;

% features: (1) maxDtheta, (2) maxDphi, (3) maxDkappaH, (4) maxDkappaV, (5) slide distance,
% (6) touch duration, (7) theta, (8) phi, (9) kappaH, (10) kappaV, (11) arc length, (12) touch counts
sFeatureNames = {'protractionTouchDThetaByWhisking', 'protractionTouchDPhiByWhisking', 'protractionTouchDKappaHByWhisking', ...
    'protractionTouchDKappaVByWhisking', 'protractionTouchSlideDistanceByWhisking', 'protractionTouchDurationByWhisking'};
mFeatureNames = {'theta', 'phi', 'kappaH', 'kappaV', 'arcLength'};
objectAngles = 45:15:135;

% first, gather feature distribution (absolute values)
featureDist = cell(12,length(mice), length(objectAngles), 2); % separated by sessions first, and then combine later
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
        
        for ai = 1 : length(objectAngles)
            angle = objectAngles(ai);
            angleInds = find(cellfun(@(x) x.angle == angle, u.trials(testInds)));
            for fi = 1 : length(sFeatureNames)            
                featureDist{fi,mi,ai,si} = cellfun(@(x,y) nanmean(x.(sFeatureNames{fi})(y)), u.trials(testInds(angleInds)), touchBeforeFirstLickInds(angleInds));
                % Nanmean in cellfun is to remove NaN effect from whisker feature calculation
                % Values can have NaN in non-touch trials
            end
            for fi = 1 : length(mFeatureNames)
                featureDist{fi+length(sFeatureNames), mi, ai, si} = cellfun(@(x,y) nanmean(x.(mFeatureNames{fi})(cellfun(@(z) z(1), x.protractionTouchChunksByWhisking(y)))), ...
                    u.trials(testInds(angleInds)), touchBeforeFirstLickInds(angleInds));
            end
            % Touch count can have 0
            featureDist{12,mi,ai,si} = cellfun(@length, touchBeforeFirstLickInds(angleInds));
        end
    end
end

%%
meanFeaturesAbs = zeros(12,length(mice),length(objectAngles),2);
% features standardized within each session: for separation between angles
meanFeaturesStdSession = zeros(12,length(mice),length(objectAngles),2);
% features standardized across 2 sessions from each mouse: for separation between sessions within each mouse
meanFeaturesStdMouse = zeros(12,length(mice),length(objectAngles),2);
for mi = 1 : length(mice)
    for si = 1 : length(sessions{mi})
        for fi = 1 : 12
            meanFeaturesAbs(fi,mi,:,si) = cellfun(@nanmean, squeeze(featureDist(fi,mi,:,si)));
            % Nanmean is to remove no-touch trials. Not affecting touch counts 
            
            featureMeanSession = nanmean(cell2mat(squeeze(featureDist(fi,mi,:,si))));
            featureStdSession = nanstd(cell2mat(squeeze(featureDist(fi,mi,:,si))));
            meanFeaturesStdSession(fi,mi,:,si) = cellfun(@nanmean, cellfun(@(x) (x-featureMeanSession)/featureStdSession, squeeze(featureDist(fi,mi,:,si)), 'un', 0));
            
            featureDistMouse = [cell2mat(squeeze(featureDist(fi,mi,:,1))); cell2mat(squeeze(featureDist(fi,mi,:,2)))];
            featureMeanMouse = nanmean(featureDistMouse);
            featureStdMouse = nanstd(featureDistMouse);
            meanFeaturesStdMouse(fi,mi,:,si) = cellfun(@nanmean, cellfun(@(x) (x-featureMeanMouse)/featureStdMouse, squeeze(featureDist(fi,mi,:,si)), 'un', 0));
        end
    end
end

%%

%%
subplotPos = [3,4,7,8,12,11,1,2,5,6,9,10];
titleStrs = {'Push angle (max\Delta\theta)','Vertical displacement (max\Delta\phi)','Horizontal bending (max\Delta\Kappa_H)','Vertical bending (max\Delta\Kappa_V)','Slide distance', 'Touch duration', ...
    'Azimuthal angle (\theta)', 'Vertical angle (\phi)', 'Horizontal curvature (\kappa_H)', 'Vertical curvature (\kappa_V)', 'Arc length', 'Touch count'};
figure, hold on
for fi = 1 : 12
    tempNaive = squeeze(meanFeaturesStdSession(fi,:,:,1));
    tempExpert = squeeze(meanFeaturesStdSession(fi,:,:,2));
%     if ismember(fi, [3,4,9,10])
%         tempNaive = -tempNaive;
%         tempExpert = -tempExpert;
%     end
    subplot(3,4,subplotPos(fi)), hold on
    plot(objectAngles, mean(tempNaive), 'color', colorsTransient(1,:));
    plot(objectAngles, mean(tempExpert), 'color', colorsTransient(2,:));
    if fi == 1
        legend({'Naive', 'Expert'}, 'autoupdate', 'off')
    end
    boundedline(objectAngles, mean(tempNaive), sem(tempNaive), 'cmap', colorsTransient(1,:))
    boundedline(objectAngles, mean(tempExpert), sem(tempExpert), 'cmap', colorsTransient(2,:))
    title(titleStrs{fi})
    xticks(45:45:135)
    ylim([-2 2])
end
sgtitle('Mean session-standardized feature values')



%%
%% 2D Contingency table
%% 2E Correct rate
%% 2F Error in angle prediction
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or '', or 'NaiveAll' for all 12 mice
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
dir = [baseDir, 'Whisker\Models\'];

expert = load([dir, fn]);

learned = 'Naive';
fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];

naive = load([dir, fn]);

nIter = length(naive.groupMdl{1}.fitCoeffs);
Ypairs = cell(6,2);
for gi = 1 : 6

    dataX = [ones(size(naive.groupMdl{gi}.io.X,1),1),naive.groupMdl{gi}.io.X];
    dataY = naive.groupMdl{gi}.io.Y;
    coeffs = zeros(size(naive.groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
        coeffs = coeffs + naive.groupMdl{gi}.fitCoeffs{ii}/nIter;
    end

    probMat = dataX * coeffs;
    probY = exp(probMat) ./ sum(exp(probMat),2);

    [~, maxind] = max(probY,[],2);
    predY = maxind * 15 + 30;

    Ypairs{gi, 1} = [dataY, predY];
    
    dataX = [ones(size(expert.groupMdl{gi}.io.X,1),1),expert.groupMdl{gi}.io.X];
    dataY = expert.groupMdl{gi}.io.Y;
    coeffs = zeros(size(expert.groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
        coeffs = coeffs + expert.groupMdl{gi}.fitCoeffs{ii}/nIter;
    end

    probMat = dataX * coeffs;
    probY = exp(probMat) ./ sum(exp(probMat),2);

    [~, maxind] = max(probY,[],2);
    predY = maxind * 15 + 30;

    Ypairs{gi, 2} = [dataY, predY];
end


%% Bootsrap to have 100 samples in each real angle (total 700 samples)
% And then calculate the performance (correct rate & error angle)
% Iterate 100 times and then take the mean
% Also calculate chance levels by shuffling matches 
angles = 45:15:135;
numSamplesBS = 100;
nIterBS = 100;
nShuffle = 100;
correctRate = zeros(6,2);
chanceCR = zeros(6,2);
errorAngle = zeros(6,2);
chanceEA = zeros(6,2);
confMat = cell(6,2);
for gi = 1 : 6
    for si = 1 : 2
        tempPair = Ypairs{gi,si};
        tempCR = zeros(nIterBS, 1);
        tempChanceCR = zeros(nIterBS,1);
        tempEA = zeros(nIterBS, 1);
        tempChanceEA = zeros(nIterBS,1);
        tempConfMat = zeros(length(angles), length(angles), nIterBS);
        angleInds = cell(length(angles),1);
        for ai = 1 : length(angles)
            angleInds{ai} = find(tempPair(:,1)==angles(ai));
        end
        for ii = 1 : nIterBS
            tempIterPair = zeros(numSamplesBS * length(angles), 2);
            for ai = 1 : length(angles)
                % bootstrapping
                tempInds = randi(length(angleInds{ai}),[numSamplesBS,1]);
                inds = angleInds{ai}(tempInds);
                tempIterPair( (ai-1)*numSamplesBS+1:ai*numSamplesBS, : ) = tempPair(inds,:);
            end
            tempCR(ii) = length(find(tempIterPair(:,2) - tempIterPair(:,1)==0)) / (numSamplesBS * length(angles));
            tempEA(ii) = mean(abs(tempIterPair(:,2) - tempIterPair(:,1)));
            
            tempTempCR = zeros(nShuffle,1);
            tempTempEA = zeros(nShuffle,1);
            for shuffi = 1 : nShuffle
                shuffledPair = [tempIterPair(:,1), tempIterPair(randperm(size(tempIterPair,1)),2)];
                tempTempCR(shuffi) = length(find(shuffledPair(:,2) - shuffledPair(:,1)==0)) / (numSamplesBS * length(angles));
                tempTempEA(shuffi) = mean(abs(shuffledPair(:,2) - shuffledPair(:,1)));
            end
            tempChanceCR(ii) = mean(tempTempCR);
            tempChanceEA(ii) = mean(tempTempEA);
            
            tempConfMat(:,:,ii) = confusionmat(tempIterPair(:,1), tempIterPair(:,2))/numSamplesBS;
        end
        correctRate(gi,si) = mean(tempCR);
        chanceCR(gi,si) = mean(tempChanceCR);
        errorAngle(gi,si) = mean(tempEA);
        chanceEA(gi,si) = mean(tempChanceEA);
        confMat{gi,si} = mean(tempConfMat,3);
    end
end

%%
%% 2D - Contingency tables
%%
contMatNaive = zeros(7,7,6);
for i = 1 : 6
    contMatNaive(:,:,i) = confMat{i,1};
end

figure, imagesc(mean(contMatNaive,3),[0 0.85]), axis square, colorbar
xticklabels(cellfun(@num2str, num2cell(45:15:135), 'un', 0))
yticklabels(cellfun(@num2str, num2cell(45:15:135), 'un', 0))
xlabel('Predicted object angle (\circ)')
ylabel('Object angle (\circ)')
title('Naive')

contMatExpert = zeros(7,7,6);
for i = 1 : 6
    contMatExpert(:,:,i) = confMat{i,2};
end

figure, imagesc(mean(contMatExpert,3),[0 0.85]), axis square, colorbar
xticklabels(cellfun(@num2str, num2cell(45:15:135), 'un', 0))
yticklabels(cellfun(@num2str, num2cell(45:15:135), 'un', 0))
xlabel('Predicted object angle (\circ)')
ylabel('Object angle (\circ)')
title('Expert')


%% 
%% 2E - Figure Classification performance - correct rate
%%

figure('units','norm','pos',[0.2 0.2 0.2 0.4]), hold on,
for i = 1 : 6
    plot(correctRate(i,:), 'ko-')
end
errorbar([1,2],mean(correctRate), sem(correctRate), 'ro', 'lines', 'no')
errorbar([1,2],mean(chanceCR), sem(chanceCR), 'o', 'lines', 'no', 'color', [0.6 0.6 0.6])

[~,p] = ttest(diff(correctRate,1,2))
% p = 0.1367
title('Correct rate')
xlim([0.5 2.5]), xticks([1,2]), xticklabels({'Naive', 'Expert'})
ylim([0 1])
yticks([0:0.1:1])
ylabel('Classification performance')
set(gca,'fontsize',12, 'fontname', 'Arial')

%%
%% 2F - Classification performance - prediction error
%%

figure('units','norm','pos',[0.2 0.2 0.2 0.4]), hold on,
for i = 1 : 6
    plot(errorAngle(i,:), 'ko-')
end
errorbar([1,2],mean(errorAngle), sem(errorAngle), 'ro', 'lines', 'no')
errorbar([1,2],mean(chanceEA), sem(chanceEA), 'o', 'lines', 'no', 'color', [0.6 0.6 0.6])

[~,p] = ttest(diff(errorAngle,1,2));

title({'Prediction error'; sprintf('%f',p)})
xlim([0.5 2.5]), xticks([1,2]), xticklabels({'Naive', 'Expert'})
ylim([0 35])
yticks([0:5:35])
ylabel('Prediction error (o)')
set(gca,'fontsize',12, 'fontname', 'Arial')


%%
%% 2G - Mouse performance VS Ideal observer
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or '', or 'NaiveAll' for all 12 mice
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
dir = [baseDir, 'Whisker\Models\'];

load([dir, fn])

%%
%% S2G - Individual mouse
%%
figure,
correctRateMouseExpert = zeros(6,1);
correctRateModelExpert = zeros(6,1);
for mi = 1 : 6
    temp = groupMdl{mi}.outcomes.matrix;
    inds = intersect(find(temp(7,:)), find(temp(5,:) > -1));
    trialAngles = temp(6,inds);
    Llicks = temp(4,inds);
    angles = unique(trialAngles);
    pLlickData = zeros(1,length(angles));
    model = cellfun(@(x) [ones(size(groupMdl{mi}.io.X,1),1), groupMdl{mi}.io.X] * x, groupMdl{mi}.fitCoeffs, 'un', 0);
    predAngleInd = zeros(length(inds), length(model));
    for i = 1 : length(model)
        [~, predAngleInd(:,i)] = max(model{i}, [], 2);
    end
    LlickModelBin = zeros(size(predAngleInd));
    LlickModelBin(find(predAngleInd > 4)) = 1;
    ind90d = find(predAngleInd(:) == 4);
    LlickModelBin(ind90d) = deal(0.5);
    pLlickModel = zeros(length(model), length(angles));

    for ai = 1 : length(angles)
        tempInds = find(trialAngles == angles(ai));
        pLlickData(ai) = sum(Llicks(tempInds)) / length(tempInds);
        pLlickModel(:,ai) = sum(LlickModelBin(tempInds,:)) / length(tempInds);
    end
    
    correctRateMouseExpert(mi) = (sum(1-pLlickData(1:3)) + sum(pLlickData(5:7))) / 6;
    correctRateModelExpert(mi) = (sum(sum(1-pLlickModel(:,1:3))) + sum(sum(pLlickModel(:,5:7)))) / 60;

    subplot(2,3,mi), hold on
    boundedline(angles, mean(pLlickModel), sem(pLlickModel), 'k')
    plot(angles, pLlickData, 'r')
    xticks([45:30:135])
    xlim([40,140])
end


%%
%% 2G
%%

figure, hold on
plot([1,2],[correctRateMouseExpert, correctRateModelExpert], 'ko-')
errorbar(1,mean(correctRateMouseExpert), sem(correctRateMouseExpert), 'ro')
errorbar(2,mean(correctRateModelExpert), sem(correctRateModelExpert), 'ro')

xlim([0.5 2.5])
xticks([1:2])
xticklabels({'Mouse', 'Ideal observer'})
ylim([0.5 1])
ylabel('Average correct rate')

[~,p] = paired_test(correctRateMouseExpert, correctRateModelExpert)




%%
%% 2H top - Variable importance of angle prediction from 7-angle sessions (Naive)
%%
Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Naive'; % 'Naive', 'Expert', or '', or 'NaiveAll' for all 12 mice
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
dir = [baseDir, 'Whisker\Models\'];
load([dir, fn])

angles = 45:15:135;
nGroup = length(groupMdl);
nIter = length(groupMdl{1}.fitCoeffs);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = zeros(size(groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
        coeffs = coeffs + groupMdl{gi}.fitCoeffs{ii}/nIter;
    end
    
    predMat = dataX * coeffs;
    predY = exp(predMat) ./ sum(exp(predMat),2);
    if find(abs(sum(predY,2)-1)>10^(-10))
        error('prob does not sum to 1')
    end
    
    patternY = zeros(length(dataY),length(listOfY));
    for pi = 1 : size(patternY,1)
        patternY(pi,find(angles==dataY(pi))) = 1;
    end
    
    nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
    satLL = 0;    
    LL(gi,1) = sum(log(mnpdf(patternY,predY)));
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat) ./ sum(exp(predMat),2);
        if find(abs(sum(predY,2)-1)>10^(-10))
            error('prob does not sum to 1')
        end
        LL(gi,fi+1) = sum(log(mnpdf(patternY,predY)));
    end
    DE(gi,:) = (LL(gi,:) - nuLL)/(satLL-nuLL);
end

tempMat = DE(:,1) - DE(:,2:end);
meanVal = DE(:,1) - mean(DE(:,2:end),2);
featureXpos = [8,9,10,11,13,12,1,2,3,4,5,6];

figure, hold on
bar(1, mean(meanVal), 'facecolor', [0.6 0.6 0.6])
bar(featureXpos+2, mean(tempMat), 'facecolor', colorsTransient(1,:))
errorbar(1, mean(meanVal), sem(meanVal), 'k')
errorbar(featureXpos+2, mean(tempMat), sem(tempMat), 'k', 'linestyle', 'none')
xlim([0 16])
xticks(sort([1, featureXpos+2]))
featureOrder = [7,8,9,10,11,12,1,2,3,4,6,5];
xticklabels(['Average', {groupMdl{1}.fitCoeffsFields{featureOrder}}])
xtickangle(45)
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Variable importance')
box('off')

pnaive = zeros(1,12);
for i = 1 : length(featureOrder)
    [~,pnaive(i)] = ttest(meanVal, tempMat(:,featureOrder(i)));
end






%%
%% 2H bottom - Variable importance of angle prediction from 7-angle sessions (Expert)
%%

Xhow = 'Mean'; %'Individual' or 'Mean'
Yout = 'Touch'; % 'Touch' or 'Choice'
learned = 'Expert'; % 'Naive', 'Expert', or '', or 'NaiveAll' for all 12 mice
task = 'Discrete'; % 'Two', 'Discrete', or 'RadialDistance'
timing = 'lick'; % 'lick' or 'answer'

fn = ['mdl', task, learned, Xhow, Yout, '_new_', timing];
dir = [baseDir, 'Whisker\Models\'];
load([dir, fn])

angles = 45:15:135;
nGroup = length(groupMdl);
nIter = length(groupMdl{1}.fitCoeffs);
nFeatures = length(groupMdl{1}.fitCoeffsFields);
LL = zeros(nGroup, 1+nFeatures);
DE = zeros(nGroup, 1+nFeatures);
for gi = 1 : nGroup
    dataX = [ones(size(groupMdl{gi}.io.X,1),1),groupMdl{gi}.io.X];
    dataY = groupMdl{gi}.io.Y;
    listOfY = unique(dataY);
    coeffs = zeros(size(groupMdl{gi}.fitCoeffs{1}));
    for ii = 1 : nIter
        coeffs = coeffs + groupMdl{gi}.fitCoeffs{ii}/nIter;
    end
    
    predMat = dataX * coeffs;
    predY = exp(predMat) ./ sum(exp(predMat),2);
    if find(abs(sum(predY,2)-1)>10^(-10))
        error('prob does not sum to 1')
    end
    
    patternY = zeros(length(dataY),length(listOfY));
    for pi = 1 : size(patternY,1)
        patternY(pi,find(angles==dataY(pi))) = 1;
    end
    
    nuLL = sum(log(mnpdf(patternY,ones(size(patternY))/length(listOfY))));
    satLL = 0;    
    LL(gi,1) = sum(log(mnpdf(patternY,predY)));    
    for fi = 1 : nFeatures
        tempDataX = dataX(:,setdiff(1:nFeatures+1, fi+1));
        tempCoeffs = coeffs(setdiff(1:nFeatures+1, fi+1),:);
        predMat = tempDataX * tempCoeffs;
        predY = exp(predMat) ./ sum(exp(predMat),2);
        if find(abs(sum(predY,2)-1)>10^(-10))
            error('prob does not sum to 1')
        end
        LL(gi,fi+1) = sum(log(mnpdf(patternY,predY)));
    end
    DE(gi,:) = (LL(gi,:) - nuLL)/(satLL-nuLL);
end

tempMat = DE(:,1) - DE(:,2:end);
meanVal = DE(:,1) - mean(DE(:,2:end),2);
featureXpos = [8,9,10,11,13,12,1,2,3,4,5,6];

figure, hold on
bar(1, mean(meanVal), 'facecolor', [0.6 0.6 0.6])
bar(featureXpos+2, mean(tempMat), 'facecolor', colorsTransient(2,:))
errorbar(1, mean(meanVal), sem(meanVal), 'k')
errorbar(featureXpos+2, mean(tempMat), sem(tempMat), 'k', 'linestyle', 'none')
xlim([0 16])
xticks(sort([1, featureXpos+2]))
featureOrder = [7,8,9,10,11,12,1,2,3,4,6,5];
xticklabels(['Average', {groupMdl{1}.fitCoeffsFields{featureOrder}}])
xtickangle(45)
ylim([0 0.4]), yticks(0:0.1:0.4)
ylabel('Variable importance')
box('off')

pexpert = zeros(1,12);
for i = 1 : length(featureOrder)
    [~,pexpert(i)] = ttest(meanVal, tempMat(:,featureOrder(i)));
end
    



%%
%%
%% Figure S2C - whisker 3D reconstruction
%%
mouse = 25;
session = 4;
whiskerDir = [baseDir, 'Whisker\'];
% Load whisker 3D files as an array, and find out the trial with biggest
% amplitude (roughly)

mouseName = sprintf('JK%03d',mouse);
sessionName = sprintf('S%02d',session);
sampleDir = [whiskerDir, mouseName, sessionName];
w3a = Whisker.Whisker3D_2padArray(sampleDir);
wfa = Whisker.WhiskerFinal_2padArray(sampleDir);

% show 3D plot, including -50 and +50 timepoints
trial45ind = find(cellfun(@(x) x.servoAngle == 45, w3a.trials));
maxAmpList = cellfun(@(x) max(smooth(x.theta(x.poleUpFrames))) - min(smooth(x.theta(x.poleUpFrames))), w3a.trials(trial45ind));
[~, wfind] = sort(maxAmpList, 'descend');

sortedi = 8;
w3 = w3a.trials{trial45ind(wfind(sortedi))};

wfind = find(cellfun(@(x) x.trialNum == w3.trialNum, wfa.trials));
wf = wfa.trials{wfind};
framesToUse = [1501:1800];
touchFrames = intersect(cell2mat(wf.protractionTFchunks'), framesToUse);

nonTouchFrames = setdiff(framesToUse, touchFrames);
refPoint = mean(w3.base);

figure, hold all
for i = 1 : length(nonTouchFrames)
    plot3((w3.trackerData{nonTouchFrames(i)}(:,1)-refPoint(1))/w3.pxPerMm, (w3.trackerData{nonTouchFrames(i)}(:,2)-refPoint(2))/w3.pxPerMm, ...
        (w3.trackerData{nonTouchFrames(i)}(:,3)-refPoint(3))/w3.pxPerMm, '-', 'color', [0.6 0.6 0.6])
    plot3((w3.base(nonTouchFrames(i),1)-refPoint(1))/w3.pxPerMm, (w3.base(nonTouchFrames(i),2)-refPoint(2))/w3.pxPerMm, ...
        (w3.base(nonTouchFrames(i),3)-refPoint(3))/w3.pxPerMm, 'r.', 'markersize', 20)
end
for i = 1 : length(touchFrames)
    plot3((w3.trackerData{touchFrames(i)}(:,1)-refPoint(1))/w3.pxPerMm, (w3.trackerData{touchFrames(i)}(:,2)-refPoint(2))/w3.pxPerMm, ...
        (w3.trackerData{touchFrames(i)}(:,3)-refPoint(3))/w3.pxPerMm, '-', 'color', [0 0 0])
    plot3((w3.base(touchFrames(i),1)-refPoint(1))/w3.pxPerMm, (w3.base(touchFrames(i),2)-refPoint(2))/w3.pxPerMm, ...
        (w3.base(touchFrames(i),3)-refPoint(3))/w3.pxPerMm, 'r.', 'markersize', 20)
end
ytickVals = yticks();
ytickVals = -ytickVals;
yticklabels(ytickVals)
xlabel('AP (mm)'), ylabel('ML (mm)'), zlabel('DV (mm)')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%%
%% S2D - Touch hyperplane
%%
mouse = 25;
session = 2;
mouseName = sprintf('JK%03d',mouse);
sessionName = sprintf('S%02d',session);
sampleDir = [whiskerDir, mouseName, sessionName];
wsa = Whisker.WhiskerSignalTrialArray_2pad(sampleDir);

trial45ind = find(cellfun(@(x) x.angle == 45, wsa.trials));
x = cell2mat(cellfun(@(x) x.whiskerEdgeCoord(:,1), wsa.trials(trial45ind)', 'uniformoutput', false));
y = cell2mat(cellfun(@(x) x.whiskerEdgeCoord(:,2), wsa.trials(trial45ind)', 'uniformoutput', false));
z = cell2mat(cellfun(@(x) ones(size(x.whiskerEdgeCoord,1),1) * x.apUpPosition, wsa.trials(trial45ind)', 'uniformoutput', false));
figure, plot3(x,y,z, 'k.', 'markersize', 0.01)



%%
%% S2E - 45 degrees example
%%

mouse = 25;
session = 4;
calciumDir = [baseDir, 'Calcium\'];
uberFn = sprintf('%s\\%03d\\UberJK%03dS%02d_NC', calciumDir, mouse, mouse, session);

load(uberFn)

touchInd = find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials));
correctInd = find(cellfun(@(x) length(x.drinkingOnsetTime), u.trials));
trial45ind = find(cellfun(@(x) x.angle == 45, u.trials));
exInd = intersect(intersect(touchInd, trial45ind), correctInd);
[~, max45ind] = sort(cellfun(@(x) nanstd(x.theta(ismember(x.whiskerTime,x.poleUpTime))), u.trials(exInd)), 'descend');
i = 8;
utrial = u.trials{exInd(max45ind(i))};

colorTemp = jet(25);
colorAngle = mat2cell(colorTemp([1,8],:),[1 1]);
colorCurve = mat2cell(colorTemp([20,25],:),[1 1]);
figure('unit','normalized','outerposition',[0 0.3 1 0.4]) 
h1 = subplot(211); hold on, 
ylimAngle = [min(union(utrial.theta, utrial.phi)), max(union(utrial.theta, utrial.phi))];
touchFrames = cell2mat(utrial.protractionTouchChunksByWhisking');
for i = 1 : length(touchFrames)
    timepoint = utrial.whiskerTime(touchFrames(i));
    plot([timepoint timepoint], [ylimAngle(1) ylimAngle(2)], '-', 'color', [0.8 0.8 0.8])
end
plot(utrial.whiskerTime, utrial.theta, '-', 'color', colorAngle{1})
ylabel('\theta')
yyaxis right,
plot(utrial.whiskerTime, utrial.phi, '-', 'color', colorAngle{2})
ylabel('\phi')
ax = gca;
ax.YAxis(1).Color = colorAngle{1};
ax.YAxis(2).Color = colorAngle{2};
ax.YAxis(1).Limits = ylimAngle;
ax.YAxis(2).Limits = ylimAngle;
xticks([])
set(gca, 'fontsize', 14)
xlim([min(utrial.whiskerTime), max(utrial.whiskerTime)])

dKh = zeros(1,length(utrial.whiskerTime));
dKv = zeros(1,length(utrial.whiskerTime));
for i = 1 : length(utrial.protractionTouchChunksByWhisking)
    tempInds = utrial.protractionTouchChunksByWhisking{i};
    for j = 2 : length(tempInds)
        dKh(tempInds(j)) = utrial.kappaH(tempInds(j)) - utrial.kappaH(tempInds(1));
        dKv(tempInds(j)) = utrial.kappaV(tempInds(j)) - utrial.kappaV(tempInds(1));
    end
end

subplot(212), hold on, 
% ylimCurve = [min(union(utrial.kappaH, utrial.kappaV)), max(union(utrial.kappaH, utrial.kappaV))];
ylimCurve = [min(union(dKh, dKv)), max(union(dKh, dKv))];
for i = 1 : length(touchFrames)
    timepoint = utrial.whiskerTime(touchFrames(i));
    plot([timepoint timepoint], [ylimCurve(1) ylimCurve(2)], '-', 'color', [0.8 0.8 0.8])
end
plot(utrial.whiskerTime, dKh, '-', 'color', colorCurve{1}),
ylabel('\Delta\kappa_H')
yyaxis right, plot(utrial.whiskerTime, dKv, '-', 'color', colorCurve{2}),
ylabel('\Delta\kappa_V')
ax = gca;
ax.YAxis(1).Color = colorCurve{1};
ax.YAxis(2).Color = colorCurve{2};
ax.YAxis(1).Limits = ylimCurve;
ax.YAxis(2).Limits = ylimCurve;
set(gca, 'fontsize', 14, 'fontname', 'Arial')
xlim([min(utrial.whiskerTime), max(utrial.whiskerTime)])
% xlim([1 4])
xlabel('Time (s)')

%%
%% S2F - 135 degrees example
%%
trial135ind = find(cellfun(@(x) x.angle == 135, u.trials));
[~, max135ind] = sort(cellfun(@(x) nanstd(x.theta(ismember(x.whiskerTime,x.poleUpTime))), u.trials(trial135ind)), 'descend');
i = 6; % JK025 S04 T473
utrial = u.trials{trial135ind(max135ind(i))};
colorTemp = jet(25);
colorAngle = mat2cell(colorTemp([1,8],:),[1 1]);
colorCurve = mat2cell(colorTemp([20,25],:),[1 1]);
figure('unit','normalized','outerposition',[0 0.3 1 0.4]) 
h1 = subplot(211); hold on, 
ylimAngle = [min(union(utrial.theta, utrial.phi)), max(union(utrial.theta, utrial.phi))];
touchFrames = cell2mat(utrial.protractionTouchChunksByWhisking');
for i = 1 : length(touchFrames)
    timepoint = utrial.whiskerTime(touchFrames(i));
    plot([timepoint timepoint], [ylimAngle(1) ylimAngle(2)], '-', 'color', [0.8 0.8 0.8])
end
plot(utrial.whiskerTime, utrial.theta, '-', 'color', colorAngle{1})
ylabel('\theta')
yyaxis right,
plot(utrial.whiskerTime, utrial.phi, '-', 'color', colorAngle{2})
ylabel('\phi')
ax = gca;
ax.YAxis(1).Color = colorAngle{1};
ax.YAxis(2).Color = colorAngle{2};
ax.YAxis(1).Limits = ylimAngle;
ax.YAxis(2).Limits = ylimAngle;
xticks([])
set(gca, 'fontsize', 14)
xlim([min(utrial.whiskerTime), max(utrial.whiskerTime)])

dKh = zeros(1,length(utrial.whiskerTime));
dKv = zeros(1,length(utrial.whiskerTime));
for i = 1 : length(utrial.protractionTouchChunksByWhisking)
    tempInds = utrial.protractionTouchChunksByWhisking{i};
    for j = 2 : length(tempInds)
        dKh(tempInds(j)) = utrial.kappaH(tempInds(j)) - utrial.kappaH(tempInds(1));
        dKv(tempInds(j)) = utrial.kappaV(tempInds(j)) - utrial.kappaV(tempInds(1));
    end
end

subplot(212), hold on, 
ylimCurve = [min(union(dKh, dKv)), max(union(dKh, dKv))];
for i = 1 : length(touchFrames)
    timepoint = utrial.whiskerTime(touchFrames(i));
    plot([timepoint timepoint], [ylimCurve(1) ylimCurve(2)], '-', 'color', [0.8 0.8 0.8])
end
plot(utrial.whiskerTime, dKh, '-', 'color', colorCurve{1}),
ylabel('\Delta\kappa_H')
yyaxis right, plot(utrial.whiskerTime, dKv, '-', 'color', colorCurve{2}),
ylabel('\Delta\kappa_V')
ax = gca;
ax.YAxis(1).Color = colorCurve{1};
ax.YAxis(2).Color = colorCurve{2};
ax.YAxis(1).Limits = ylimCurve;
ax.YAxis(2).Limits = ylimCurve;
set(gca, 'fontsize', 14, 'fontname', 'Arial')
xlim([min(utrial.whiskerTime), max(utrial.whiskerTime)])
xlabel('Time (s)')





