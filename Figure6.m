
baseDir = 'D:\TPM\JK\Pub_S1AngleCode\'; % The folder containing folders of data ('\Behavior', '\Calcium', '\Whisker') and dependent codes ('\MATLAB codes')
%% basic settings

calciumDir = [baseDir, 'Calcium\'];
%%
%% Fig 6A - Example FOV's for cell matching
%%
mouse = 39;
naiveSession = 1;
plane = 5;
load(sprintf('%s%03d\\cellIDmatch_JK%03d', calciumDir, mouse, mouse), 'us')
load(sprintf('%s%03d\\UberJK%03dS%02d_NC', calciumDir, mouse, mouse, naiveSession), 'u')
naiveIm = us.sessions(1).mimg{plane};
expertIm = us.sessions(2).mimg{plane};
tform = us.sessions(1).tform{plane};
movedIm = imwarp(naiveIm, tform, 'outputview', imref2d(size(expertIm)));
naiveCellmap = us.sessions(1).cellmap{plane};
expertCellmap = us.sessions(2).cellmap{plane};

naiveCells = setdiff(unique(naiveCellmap),0);
expertCells = setdiff(unique(expertCellmap),0);
movedCellmap = imwarp(naiveCellmap, tform, 'outputview', imref2d(size(expertCellmap)), 'interp', 'nearest');

% not-registered naive
% figure, 
% imshow((mat2gray((naiveIm))))
% hold on
% for ci = 1 : length(naiveCells)
%     cID = naiveCells(ci);
%     [ypix,xpix] = ind2sub(size(naiveCellmap),find(naiveCellmap == cID));
%     k = boundary(ypix, xpix);
%     patch(xpix(k), ypix(k), 'm', 'edgecolor', 'none')
% end
% 
% scaleBarPix = 100 / u.pixResolution;
% margins = 20; % 20 pixels each away from right bottom corner
% sizes = size(naiveIm);
% plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)
%% registered naive
figure, 
imshow(mat2gray(movedIm))
hold on
for ci = 1 : length(naiveCells)
    cID = naiveCells(ci);
    [ypix,xpix] = ind2sub(size(movedCellmap),find(movedCellmap == cID));
    k = boundary(ypix, xpix);
    patch(xpix(k), ypix(k), 'm', 'edgecolor', 'none')
end

scaleBarPix = 100 / u.pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(movedIm);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)

%% Expert FOV
figure,
imshow(mat2gray(expertIm))
hold on
for ci = 1 : length(expertCells)
    cID = expertCells(ci);
    [ypix,xpix] = ind2sub(size(expertCellmap),find(expertCellmap == cID));
    k = boundary(ypix, xpix);
    patch(xpix(k), ypix(k), 'c', 'edgecolor', 'none')
end

scaleBarPix = 100 / u.pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(expertIm);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)



%% Merged
figure, 
imshow(zeros(size(expertIm)))
hold on
for ci = 1 : length(naiveCells)
    cID = naiveCells(ci);
    [ypix,xpix] = ind2sub(size(movedCellmap),find(movedCellmap == cID));
    k = boundary(ypix, xpix);
    patch(xpix(k), ypix(k), 'm', 'edgecolor', 'none')
end
for ci = 1 : length(expertCells)
    cID = expertCells(ci);
    [ypix,xpix] = ind2sub(size(expertCellmap),find(expertCellmap == cID));
    k = boundary(ypix, xpix);
    patch(xpix(k), ypix(k), 'c', 'edgecolor', 'none')
end

overlapCellmap = expertCellmap .* imbinarize(movedCellmap);
for ci = 1 : length(expertCells)
    cID = expertCells(ci);
    [ypix,xpix] = ind2sub(size(overlapCellmap),find(overlapCellmap == cID));
    k = boundary(ypix, xpix);
    patch(xpix(k), ypix(k), 'w', 'edgecolor', 'none')
end



%% zoomed in images
before = movedIm(200:400, 200:400);
after = expertIm(200:400, 200:400);
% before = (before);
% after = imhistmatch(after, before);
figure, imshow(mat2gray(before)), hold on
scaleBarPix = 50 / u.pixResolution;
margins = 10; % 20 pixels each away from right bottom corner
sizes = size(before);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)

figure, imshow(mat2gray(after)), hold on
scaleBarPix = 50 / u.pixResolution;
margins = 10; % 20 pixels each away from right bottom corner
sizes = size(after);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)



%%
%% FOV matching methods and controls
%% Fig S6A-C
%%



%%

mouse = 39;
refSession = 1;
zSession = 998;

load(sprintf('%s%03d\\zstack_%03d_%03d', calciumDir, mouse, mouse, zSession), 'knobbyInfo', 'zstack')
load(sprintf('%s%03d\\regops_%03d_%03d', calciumDir, mouse, mouse, refSession), 'ops1')

zComp = zstack(ops1{1}.useY, ops1{1}.useX, :);
template = ops1{1}.mimg1;

%%
%% S6D - example correlation values
%%
figure, imshow(mat2gray(template))
hold on
pixResolution = 1.4 / str2double(ops1{1}.info.config.magnification_list(ops1{1}.info.config.magnification,:));
scaleBarPix = 100 / pixResolution;
margins = 20; % 20 pixels each away from right bottom corner
sizes = size(template);
plot([sizes(2)-margins-scaleBarPix, sizes(2)-margins], [sizes(1) - margins, sizes(1) - margins], 'w-', 'linewidth', 4)

repeat = 51;
interval = 1; % 1 is equal to 2 um
refPlane = 200;
startPlane = refPlane - interval * (repeat-1)/2;
corval = zeros(1,repeat);
zCompIms = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane + interval * (i-1);
    temp = zComp(:,:,pi);
    zCompIms{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
end
%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')





%%
%% 6E - z-axis span calculation
%%

mice = [25,27,30,36,39,52];
zstackFilenum = [1000, 1000, 2000, 997, 998, 995];
%% First, check zstackReg from interactive figure and select representative neurons

mi = 2; % for example
load(sprintf('%s%03d\\zstackReg_%03d_%d',calciumDir,mice(mi),mice(mi),zstackFilenum(mi)), 'zstackReg', 'zstackDepths')
load(sprintf('%s%03d\\%03d_%d_000',calciumDir,mice(mi),mice(mi),zstackFilenum(mi)),'info')

%% Manual selection of a chronically fluorescent neuron from visual inspection
mag=str2double(info.config.magnification_list(info.config.magnification,:));
xyumPerPix = 1.4 / mag;
depthDiff = abs(diff(zstackDepths));
zumPerPix = mean(depthDiff(find(depthDiff)));
numFrames = size(zstackReg,3);
%
%
frames = 1 : numFrames;
figure
while true
    currFrame = frames(1);
    imshow(zstackReg(:,:,currFrame))
    title(sprintf('Current frame num = %d', currFrame))
    drawnow
    if waitforbuttonpress
        value = double(get(gcf,'CurrentCharacter'));
        switch value
            case 28 % <-
                frames = circshift(frames,-1);
            case 29 % ->
                frames =circshift(frames,1);
            case 30 % up arrow
                frames = circshift(frames,-10);
            case 31 % down arrow
                frames = circshift(frames,10);
            case 93 % ]
                frames = circshift(frames,-50);
            case 91 % [
                frames = circshift(frames,50);

            case 27 % when pressing esc                    
                title(sprintf('Finished (current frame num = %d)', currFrame))
                break
        end
    end
end


%% Draw x-z iamge

offset = 0;
halfspan = 20;

% cellFrames = 104:125; % JK025
cellFrames = 206:227; % JK027
% cellFrames = 143:163; % JK052

estFrames = cellFrames(1)-20:cellFrames(end)+20;

% figure, imshow(mean(zstackReg(:,:,estFrames),3))
% [centerX, centerY] = ginput(1)

% centerX = 107; centerY = 350; % JK025
centerX = 131; centerY = 224; % JK027
% centerX = 421; centerY = 242; % JK052

xypix = 50;
scalebar = 20; % in um
xyscalebar = scalebar / xyumPerPix;
zscalebar = scalebar / zumPerPix;
useX = round(centerX-xypix/2):round(centerX+xypix/2);
useY = round(centerY-xypix/2):round(centerY+xypix/2);

subVolume = zstackReg(useY,useX,estFrames);
subFrames = length(estFrames);

figure,

% subplot(211)
% imshow(mean(subVolume(:,:,round(subFrames/4):round(subFrames*3/4)),3)), axis image, hold on
% plot([xypix-5-xyscalebar, xypix-5], [xypix-5 xypix-5], 'w-', 'linewidth', 2)
% subplot(212)


% imshow(flip(squeeze(max(subVolume,[],1))',1)), axis image, hold on
imshow(flip(squeeze(mean(subVolume,1))',1)), axis image, hold on
plot(xypix-5, round(subFrames/2)+offset + round(halfspan/zumPerPix), 'k*', 'linewidth', 2)
plot(xypix-5, round(subFrames/2)+offset, 'k*', 'linewidth', 2)
plot(xypix-5, round(subFrames/2)+offset - round(halfspan/zumPerPix), 'k*', 'linewidth', 2)
plot([xypix-5, xypix-5], [subFrames-5-zscalebar subFrames-5], 'w-', 'linewidth', 2)
 

daspect([1 xyumPerPix/zumPerPix 1])
title({sprintf('cell %d:%d, frames %d:%d', cellFrames(1), cellFrames(end), estFrames(1), estFrames(end));...
    sprintf('x center: %d, y center: %d', round(centerX), round(centerY))})


%% Draw x-y images
figure('units','norm','pos',[0.3 0 0.2 1]),
subplot(511)
imshow(subVolume(:,:, round(subFrames/2)+offset + round((halfspan+10)/zumPerPix) )), axis image
title(sprintf('Center - %d um', halfspan+10))

subplot(512)
imshow(subVolume(:,:, round(subFrames/2)+offset + round((halfspan)/zumPerPix) )), axis image
title(sprintf('Center - %d um', halfspan))

subplot(513)
imshow(subVolume(:,:,round(subFrames/2)+offset)), axis image
title('Center')

subplot(514)
imshow(subVolume(:,:, round(subFrames/2)+offset - round(halfspan/zumPerPix) )), axis image, hold on
% plot([xypix-5-xyscalebar, xypix-5], [xypix-5 xypix-5], 'w-', 'linewidth', 2)
title(sprintf('Center + %d um', halfspan))

subplot(515)
imshow(subVolume(:,:, round(subFrames/2)+offset - round((halfspan+10)/zumPerPix) )), axis image, hold on
plot([xypix-5-xyscalebar, xypix-5], [xypix-5 xypix-5], 'w-', 'linewidth', 2)
title(sprintf('Center + %d um', halfspan+10))



%%
%%  S6F - compare between tdTomato and GCaMP alignment
%%
% Manually collect peak correlation from each session, and calculate their
% difference across learning. Then, compare between tdTomato-aligned
% (JK052) and GCaMP-aligned (all other) mice

mouse = 52;
refSession = 3;
zSession = 995;

load(sprintf('%s%03d\\zstack_%03d_%03d', calciumDir, mouse, mouse, zSession), 'knobbyInfo', 'zstack')
load(sprintf('%s%03d\\regops_%03d_%03d', calciumDir, mouse, mouse, refSession), 'ops1')

%%
%% For GCaMP
%%
zComp = squeeze(zstack(ops1{1}.useY, ops1{1}.useX, :, 1)); % 1 for GREEN, 2 for RED
% %%
refPlane = 120;
plane = 2;
% figure, imshow(mat2gray(zComp(:,:,refPlane)))

template = ops1{plane}.mimg1; % mimg1 for GREEN, mimgRED for RED
% template = ops1{plane}.mimgRED;
% figure, imshow(mat2gray(template));
% %%
% template = ops1{plane}.mimgRED;

repeat = 51;
interval = 1; % 1 is equal to 2 um
% refPlane = 150;
startPlane = refPlane - interval * (repeat-1)/2;
corval = zeros(1,repeat);
zCompIms = cell(repeat,1);
for i = 1 : repeat
    pi = startPlane + interval * (i-1);
    temp = zComp(:,:,pi);
    zCompIms{i} = temp;
    [u, v] = jkfftalign(temp, template);
    temp = circshift(temp, [u,v]);
    if u >= 0
        if v >= 0
            tempIm = temp(u+1:end,v+1:end);
            tempRef = template(u+1:end,v+1:end);
        else
            tempIm = temp(u+1:end,1:end+v);
            tempRef = template(u+1:end,1:end+v);
        end
    else
        if v >= 0
            tempIm = temp(1:end+u,v+1:end);
            tempRef = template(1:end+u,v+1:end);
        else
            tempIm = temp(1:end+u,1:end+v);
            tempRef = template(1:end+u,1:end+v);
        end
    end
    temprho = corrcoef(double(tempIm(:)), double(tempRef(:)));
    corval(i) = temprho(1,2); 
end
%
figure
plot([knobbyInfo(startPlane:interval:startPlane+interval*(repeat-1)).z]/cosd(35), corval(1:i), 'k.', 'markersize', 20)
xlim([knobbyInfo(startPlane-1).z/cosd(35), knobbyInfo(startPlane+interval*(repeat)).z/cosd(35)])
% ylim([0.78 0.86])   
xlabel('Depth (\mum)'), ylabel('Spatial correlation')
set(gca,'box','off')






%%
%% Fig 6B - Data for Sankey plot
%%

matchFn = 'cellMatching_beforeNafter.mat';
load(sprintf('%s%s',calciumDir, matchFn), 'match');
angleTuningFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat';
tune = load(sprintf('%s%s',calciumDir, angleTuningFn), 'expert', 'naive');

colorsTransient = [248 171 66; 40 170 225] / 255;
colorsPersistent = [1 0 0; 0 0 1];

learnerInd = [1:4,7,9];

tune.naiveLearner = tune.naive(learnerInd);
tune = rmfield(tune, 'naive');

numMice = size(match,1);

jk027length = length(tune.naiveLearner(2).touchID);
jk027inds = find(tune.naiveLearner(2).touchID > 5000);
fn = fieldnames(tune.naiveLearner);
for i = 1 : length(fn)
    if size(tune.naiveLearner(2).(fn{i}), 1) == jk027length
        tune.naiveLearner(2).(fn{i}) = tune.naiveLearner(2).(fn{i})(jk027inds,:);
    elseif size(tune.naiveLearner(2).(fn{i}), 2) == jk027length
        tune.naiveLearner(2).(fn{i}) = tune.naiveLearner(2).(fn{i})(:,jk027inds);
    end
end

numAllChanges = zeros(numMice,4,4); % 1 not-found (nf), 2 non-touch (nt), 3 not-angle-tuned (nat), 4 angle-tuned (at)
numAllCells = zeros(numMice,1);
for mi = 1 : numMice
    numAllCell = length(match{mi,1}) + length(match{mi,3}) - length(find(match{mi,2}));
    
    % start with naive first (1st dim, or y axis), and then matching each to expert (2nd min, or x axis)
    expertIdNfNaive = setdiff(match{mi,3}, match{mi,2});    
    
    expertIdNtExpert = setdiff(match{mi,3}, tune.expert(mi).touchID);
    expertIdNatExpert = tune.expert(mi).touchID(tune.expert(mi).tuned == 0);
    expertIdAtExpert = tune.expert(mi).touchID(tune.expert(mi).tuned == 1);
    
    expertIdNf2Nt = intersect(expertIdNfNaive, expertIdNtExpert);
    expertIdNf2Nat = intersect(expertIdNfNaive, expertIdNatExpert);
    expertIdNf2At = intersect(expertIdNfNaive, expertIdAtExpert);
    
    naiveIdNtNaive = setdiff(match{mi,1}, tune.naiveLearner(mi).touchID);
    naiveIdNatNaive = tune.naiveLearner(mi).touchID(tune.naiveLearner(mi).tuned == 0);
    naiveIdAtNaive = tune.naiveLearner(mi).touchID(tune.naiveLearner(mi).tuned == 1);
    
    naiveIdNfExpert = match{mi,1}(find(match{mi,2} == 0));
    
    naiveIdNt2Nf = intersect(naiveIdNtNaive, naiveIdNfExpert);
    naiveIdNat2Nf = intersect(naiveIdNatNaive, naiveIdNfExpert);
    naiveIdAt2Nf = intersect(naiveIdAtNaive, naiveIdNfExpert);
    
    % now only the matched ones are left. Stick to expert ID's. 
    naiveIndNtNaive = find(1 - ismember(match{mi,1}, tune.naiveLearner(mi).touchID));
    naiveIndNatNaive = find(ismember(match{mi,1}, tune.naiveLearner(mi).touchID(tune.naiveLearner(mi).tuned == 0)));
    naiveIndAtNaive = find(ismember(match{mi,1}, tune.naiveLearner(mi).touchID(tune.naiveLearner(mi).tuned == 1)));
    
    expertIdNtNaive = setdiff(match{mi,2}(naiveIndNtNaive), 0);
    expertIdNatNaive = setdiff(match{mi,2}(naiveIndNatNaive), 0);
    expertIdAtNaive = setdiff(match{mi,2}(naiveIndAtNaive), 0);
    
    expertIdNt2Nt = intersect(expertIdNtNaive, expertIdNtExpert);
    expertIdNt2Nat = intersect(expertIdNtNaive, expertIdNatExpert);
    expertIdNt2At = intersect(expertIdNtNaive, expertIdAtExpert);
    
    expertIdNat2Nt = intersect(expertIdNatNaive, expertIdNtExpert);
    expertIdNat2Nat = intersect(expertIdNatNaive, expertIdNatExpert);
    expertIdNat2At = intersect(expertIdNatNaive, expertIdAtExpert);
    
    expertIdAt2Nt = intersect(expertIdAtNaive, expertIdNtExpert);
    expertIdAt2Nat = intersect(expertIdAtNaive, expertIdNatExpert);
    expertIdAt2At = intersect(expertIdAtNaive, expertIdAtExpert);
    
    numAllChanges(mi,1,2) = length(expertIdNf2Nt);
    numAllChanges(mi,1,3) = length(expertIdNf2Nat);
    numAllChanges(mi,1,4) = length(expertIdNf2At);
        
    numAllChanges(mi,2,1) = length(naiveIdNt2Nf);
    numAllChanges(mi,2,2) = length(expertIdNt2Nt);
    numAllChanges(mi,2,3) = length(expertIdNt2Nat);
    numAllChanges(mi,2,4) = length(expertIdNt2At);

    numAllChanges(mi,3,1) = length(naiveIdNat2Nf);
    numAllChanges(mi,3,2) = length(expertIdNat2Nt);
    numAllChanges(mi,3,3) = length(expertIdNat2Nat);
    numAllChanges(mi,3,4) = length(expertIdNat2At);
    
    numAllChanges(mi,4,1) = length(naiveIdAt2Nf);
    numAllChanges(mi,4,2) = length(expertIdAt2Nt);
    numAllChanges(mi,4,3) = length(expertIdAt2Nat);
    numAllChanges(mi,4,4) = length(expertIdAt2At);
    
    if sum(sum(numAllChanges(mi,:,:))) ~= numAllCell
        error('total num cell does not match')
    end
    numAllCells(mi) = numAllCell;
end
%% Sanity check
allChanges = numAllChanges ./ sum(sum(numAllChanges,3),2);
sum(sum(allChanges(3,:,:)))
sum(numAllCells)




%%
%% S7A - Quantification of the overlap
%%


mice = [25,27,30,36,39,52];
numMice = length(mice);

refi = 2;
moviList = [1,3];
    
similarity = cell(numMice,1);
app = cell(numMice,1);
disapp = cell(numMice,1);

for mi = 1 : numMice
    mouse = mice(mi);
    load(sprintf('%s%03d\\cellIDmatch_JK%03d',calciumDir,mouse,mouse))
    numPlane = length(us.sessions(1).cellmap);
    
    similarity{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    app{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    disapp{mi} = nan(numPlane,length(moviList)); % 1 for 1st to 2nd, 2 for 3rd to 2nd
    
    for planei = 1 : numPlane
        for movii = 1 : length(moviList)
            movi = moviList(movii);
            cm = us.sessions(refi).cellmap{planei};
            if ~isempty(cm)
                refx(1) = min(find(sum(cm)));
                refx(2) = max(find(sum(cm)));
                refy(1) = min(find(sum(cm,2)));
                refy(2) = max(find(sum(cm,2)));
                refTemplate = zeros(size(us.sessions(refi).cellmap{planei}));
                refTemplate(refy(1):refy(2), refx(1):refx(2)) =1;

                cm = us.sessions(movi).cellmap{planei};
                movx(1) = min(find(sum(cm)));
                movx(2) = max(find(sum(cm)));
                movy(1) = min(find(sum(cm,2)));
                movy(2) = max(find(sum(cm,2)));

                movTemplate = zeros(size(us.sessions(movi).cellmap{planei}));
                movTemplate(movy(1):movy(2), movx(1):movx(2)) =1;

                % figure, imshowpair(refTemplate,movTemplate)
                movedTemplate = imwarp(movTemplate, us.sessions(movi).tform{planei}, 'outputview', imref2d(size(refTemplate)));
                movedTemplate(find(movedTemplate)) = 1;
                movedCellMapTemp = imwarp(us.sessions(movi).cellmap{planei}, us.sessions(movi).tform{planei}, 'outputview', imref2d(size(refTemplate)));
                movedCellMap = zeros(size(movedCellMapTemp));
                integerInds = find(movedCellMapTemp - floor(round(movedCellMapTemp)) == 0);
                nonZeroInds = find(movedCellMapTemp);
                finalInds = intersect(integerInds, nonZeroInds);
                movedCellMap(finalInds) = movedCellMapTemp(finalInds);
                % figure, imshowpair(refTemplate,movedTemplate)
                jointTerritory = refTemplate .* movedTemplate;
                refCells = setdiff(unique(us.sessions(refi).cellmap{planei} .* jointTerritory), 0);
                movCellTemp = setdiff(unique(movedCellMap .* jointTerritory), 0);
                id = unique(round(movCellTemp/1000));
                freq = zeros(length(id),1);
                for idi = 1 : length(id)
                    freq(idi) = length(find(floor(movCellTemp/1000) == id(idi)));
                end
                [~, maxind] = max(freq);
                maxid = id(maxind);

                inds = find(floor(movCellTemp/1000) == id(idi));
                movCells = movCellTemp(inds);

                movCi = find(ismember(us.sessions(movi).cellID, movCells));
                candCi = us.sessions(movi).matchedRefCellID(movCi);
                trackedCimov = movCi(find(candCi));
                trackedCIDmov = us.sessions(movi).cellID(trackedCimov);
                trackedCIDref = candCi(find(candCi));

                similarity{mi}(planei,movii) = length(find(candCi)) / (length(refCells) + length(movCells) - length(find(candCi)));
                app{mi}(planei,movii) = length(setdiff(refCells,trackedCIDref)) / length(refCells);
                disapp{mi}(planei,movii) = length(setdiff(movCells, trackedCIDmov)) / length(movCells);
            end
        end
    end
end

% draw
figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), similarity, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), similarity, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), similarity, 'un', 0)));

bar(1, nanmean(mat1), 'w')
bar(2, nanmean(mat2), 'w')
bar(3, nanmean(mat3), 'w')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Similarity')
ylim([0 0.6])
set(gca, 'fontsize', 12, 'fontname', 'Arial')

figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), disapp, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), disapp, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), disapp, 'un', 0)));
bar(1, nanmean(mat1), 'm')
bar(2, nanmean(mat2), 'm')
bar(3, nanmean(mat3), 'm')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Disappearred rate')
ylim([0 0.6])
set(gca, 'fontsize', 12, 'fontname', 'Arial')

figure, hold on
mat1 = (cell2mat(cellfun(@(x) nanmean(x(1:2,1)), app, 'un', 0)));
mat2 = (cell2mat(cellfun(@(x) nanmean(x(3:6,1)), app, 'un', 0)));
mat3 = (cell2mat(cellfun(@(x) nanmean(x(7:8,1)), app, 'un', 0)));
bar(1, nanmean(mat1), 'c')
bar(2, nanmean(mat2), 'c')
bar(3, nanmean(mat3), 'c')
errorbar(1, nanmean(mat1), sem(mat1), 'k')
errorbar(2, nanmean(mat2), sem(mat2), 'k')
errorbar(3, nanmean(mat3), sem(mat3), 'k')
xticks(1:3)
xticklabels({'1~2', '3~6', '7~8'})
xlabel('Imaging planes')
ylabel('Appearred rate')
ylim([0 0.6])
set(gca, 'fontsize', 12, 'fontname', 'Arial')


%% Total similarity, disappearrance, and appearrance

matSim = (cell2mat(cellfun(@nanmean, similarity, 'un', 0)));
matDis = (cell2mat(cellfun(@nanmean, disapp, 'un', 0)));
matApp = (cell2mat(cellfun(@nanmean, app, 'un', 0)));

mean(matSim)
sem(matSim)

mean(matDis)
sem(matDis)

mean(matApp)
sem(matApp)






%%
%% S7B - Within each category, how much of them fell silent after learning?
%%

changesToSilent = zeros(numMice,3); % 1 non-touch -> non-touch, 2 not-tuned -> not-tuned, 3 tuned -> tuned 
for mi = 1 : numMice
    naiveIdSilentExpert = match{mi,1}(find(match{mi,2} == 0));
    
    naiveIdTouch = tune.naiveLearner(mi).touchID;
    naiveIdNT = setdiff(match{mi,1}, naiveIdTouch);
    naiveIdNotTuned = tune.naiveLearner(mi).touchID(find(tune.naiveLearner(mi).tuned == 0));
    naiveIdTuned = tune.naiveLearner(mi).touchID(find(tune.naiveLearner(mi).tuned));
    
    changesToSilent(mi,1) = length(intersect(naiveIdNT, naiveIdSilentExpert)) / length(naiveIdNT);
    changesToSilent(mi,2) = length(intersect(naiveIdNotTuned, naiveIdSilentExpert)) / length(naiveIdNotTuned);
    changesToSilent(mi,3) = length(intersect(naiveIdTuned, naiveIdSilentExpert)) / length(naiveIdTuned);
end

colors = [ones(1,3); ones(1,3) * 0.5; zeros(1,3)];
colorsTransient = [248 171 66; 40 170 225] / 255;
figure, hold on
for i = 1 : 3
    bar(i, mean(changesToSilent(:,i)), 'facecolor', colors(i,:), 'edgecolor', colorsTransient(1,:))
    errorbar(i, mean(changesToSilent(:,i)), sem(changesToSilent(:,i)), 'k')
end

xticks(1:3)
xticklabels({'Non-touch -> silent', 'Not-tuned -> silent', 'Tuned -> silent'})
xtickangle(45)
ylabel('Proportion of naive categories')

%%
%% S7C - Within each category in expert, how much of them newly recruited?
%%

changesFromSilent = zeros(numMice,3); % 1 non-touch -> non-touch, 2 not-tuned -> not-tuned, 3 tuned -> tuned 
for mi = 1 : numMice
    expertIdSilentNaive = match{mi,3}(find(ismember(match{mi,3}, match{mi,2})==0));
    
    idTouchExpert = tune.expert(mi).touchID;
    
    expertIdNTExpert = setdiff(match{mi,3}, idTouchExpert);
    expertIdNotTunedExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned == 0));
    expertIdTunedExpert = tune.expert(mi).touchID(find(tune.expert(mi).tuned));
    
    changesFromSilent(mi,1) = length(intersect(expertIdNTExpert, expertIdSilentNaive)) / length(expertIdNTExpert);
    changesFromSilent(mi,2) = length(intersect(expertIdNotTunedExpert, expertIdSilentNaive)) / length(expertIdNotTunedExpert);
    changesFromSilent(mi,3) = length(intersect(expertIdTunedExpert, expertIdSilentNaive)) / length(expertIdTunedExpert);
end

colors = [ones(1,3); ones(1,3) * 0.5; zeros(1,3)];
figure, hold on
for i = 1 : 3
    bar(i, mean(changesFromSilent(:,i)), 'facecolor', colors(i,:), 'edgecolor', colorsTransient(2,:))
    errorbar(i, mean(changesFromSilent(:,i)), sem(changesFromSilent(:,i)), 'k')
end

xticks(1:3)
xticklabels({'Silent -> non-touch', 'Silent -> not-tuned', 'Silent -> tuned'})
xtickangle(45)
ylabel('Proportion of expert categories')



%%
%% S7D and extra - Cumulative distribution of # of inferred spikes / session
%%


%%
mice = [25,27,30,36,39,52];
numMice = length(mice);
sessions = {[4,19],[3,10],[3,21],[1,17],[1,23],[3,21]};

load([calciumDir, 'cellMatching_beforeNafter'])

colors = [248, 171, 66; ...
    40, 164, 221] / 255;

%%
persSR = cell(numMice,2); % persistent neurons' spike rates
transSR = cell(numMice,2);
persSN = cell(numMice,2); % persistent neurons' # of spikes
transSN = cell(numMice,2);
for mi = 1 : numMice
    mouse = mice(mi);
    naive = load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumDir, mouse, mouse, sessions{mi}(1)), 'u');
    expert = load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumDir, mouse, mouse, sessions{mi}(2)), 'u');
    
    % naive first
    upperTrialInds = find(cellfun(@(x) ismember(1,x.planes), naive.u.trials));
    lowerTrialInds = find(cellfun(@(x) ismember(5,x.planes), naive.u.trials));
    upperSpikes = cell2mat(cellfun(@(x) x.spk, naive.u.trials(upperTrialInds)', 'un', 0));
    lowerSpikes = cell2mat(cellfun(@(x) x.spk, naive.u.trials(lowerTrialInds)', 'un', 0));
    upperSpikeRate = mean(upperSpikes,2) * naive.u.frameRate;
    lowerSpikeRate = mean(lowerSpikes,2) * naive.u.frameRate;
    upperSpikeNum = sum(upperSpikes,2);
    lowerSpikeNum = sum(lowerSpikes,2);
    if mi == 2
        naiveSR = lowerSpikeRate;
        naiveSN = lowerSpikeNum;
    else
        naiveSR = [upperSpikeRate; lowerSpikeRate];
        naiveSN = [upperSpikeNum; lowerSpikeNum];
    end
    if length(naiveSR) ~= length(match{mi,1})
        error('Naive # of neurons does not match')
    end
    
    naivePersInd = find(match{mi,2});
    naiveTransInd = find(match{mi,2}==0);
    persSR{mi,1} = naiveSR(naivePersInd);
    transSR{mi,1} = naiveSR(naiveTransInd);
    persSN{mi,1} = naiveSN(naivePersInd);
    transSN{mi,1} = naiveSN(naiveTransInd);
    
    
    % and then expert
    upperTrialInds = find(cellfun(@(x) ismember(1,x.planes), expert.u.trials));
    lowerTrialInds = find(cellfun(@(x) ismember(5,x.planes), expert.u.trials));
    upperSpikes = cell2mat(cellfun(@(x) x.spk, expert.u.trials(upperTrialInds)', 'un', 0));
    lowerSpikes = cell2mat(cellfun(@(x) x.spk, expert.u.trials(lowerTrialInds)', 'un', 0));
    upperSpikeRate = mean(upperSpikes,2) * expert.u.frameRate;
    lowerSpikeRate = mean(lowerSpikes,2) * expert.u.frameRate;
    upperSpikeNum = sum(upperSpikes,2);
    lowerSpikeNum = sum(lowerSpikes,2);
    if mi == 2
        expertSR = lowerSpikeRate;
        expertSN = lowerSpikeNum;
    else
        expertSR =[upperSpikeRate; lowerSpikeRate];
        expertSN =[upperSpikeNum; lowerSpikeNum];
    end
    if length(expertSR) ~= length(match{mi,3})
        error('Expert # of neurons does not match')
    end
    
    expertPersID = setdiff(match{mi,2},0);
    expertTransID = setdiff(match{mi,3}, match{mi,2});
    expertPersInd = find(ismember(match{mi,3}, expertPersID));
    expertTransInd = find(ismember(match{mi,3}, expertTransID));
    persSR{mi,2} = expertSR(expertPersInd);
    transSR{mi,2} = expertSR(expertTransInd);
    persSN{mi,2} = expertSN(expertPersInd);
    transSN{mi,2} = expertSN(expertTransInd);
end


%%
%% Comparing between persistent and transient neurons
%%
histEdge = [0:0.02:1];
figure, hold on
plot(histEdge(1:end-1), histcounts(transSR{1,1}, histEdge) /(length(find(transSR{1,1}<1)) + length(find(persSR{1,1}<1))) );
plot(histEdge(1:end-1), histcounts(persSR{1,1}, histEdge) /(length(find(transSR{1,1}<1)) + length(find(persSR{1,1}<1))) );
legend({'Transient', 'Persistent'})


%%
%% S7D
%%
countEdge = 0:10:1010;
countCell = zeros(numMice,length(countEdge)-1);
for mi = 1 : numMice
    tempMat = [persSN{mi,1}; transSN{mi,1}];
    countCell(mi,:) = histcounts(tempMat, countEdge, 'norm', 'cdf');
end
figure
boundedline(countEdge(1:end-1), mean(countCell), sem(countCell), 'k')
xlabel('Number of inferred spikes / session')
ylabel('Cumulative proportion')
set(gca, 'fontsize', 12, 'fontname', 'Arial')

%%
%% S7D inset
%%
countEdge = 0:1:101;
countCell = zeros(numMice,length(countEdge)-1);
for mi = 1 : numMice
    tempMat = [persSN{mi,1}; transSN{mi,1}];
    countCell(mi,:) = histcounts(tempMat, countEdge, 'norm', 'cdf');
end
figure, hold on
boundedline(countEdge(1:end-1), mean(countCell), sem(countCell), 'k')
plot([55, 55], [0, mean(countCell(:,56))], 'k--')
plot([0, 55], [mean(countCell(:,56)), mean(countCell(:,56))], 'k--')
ylim([0 0.08])
set(gca, 'fontsize', 12, 'fontname', 'Arial')



%%
%% 6C
%%
% comparison of spike rates between persistent and transient neurons
% from all neurons (mouse-averaged; inset)
% and cumulative distribution

act = load([calciumDir, 'allActiveInfo']);
load([calciumDir, 'cellID_match_persAT_v9'], 'match');

% cleaning up nonlearners and upper layer of JK027
learnerInd = [1:4,7,9];
act.naiveLearner = act.naive(learnerInd);

jk027inds = find(act.naiveLearner(2).cellID > 5000);
jk027length = length(act.naiveLearner(2).cellID);

fn = fieldnames(act.naiveLearner);
for i = 1 : length(fn)
    if size(act.naiveLearner(2).(fn{i}),1) == jk027length
        act.naiveLearner(2).(fn{i}) = act.naiveLearner(2).(fn{i})(jk027inds,:);
    elseif length(act.naiveLearner(2).(fn{i})) == jk027length
        act.naiveLearner(2).(fn{i}) = act.naiveLearner(2).(fn{i})(:,jk027inds);
    end
end

act = rmfield(act, 'naive');

colorsTransient = [248 171 66; 40 170 225] / 255;
colorsPersistent = [1 0 0; 0 0 1];

numMice = length(learnerInd);

spikeRateNaive = cell(numMice, 2); % (:,1) persistent, (:,2) transient
spikeRateExpert = cell(numMice, 2); % (:,1) persistent, (:,2) transient
for mi = 1 : numMice
    % naive
    idPers = match{mi,1}(find(match{mi,2}));
    indPers = find(ismember(act.naiveLearner(mi).cellID, idPers));
    indTrans = setdiff(1:length(act.naiveLearner(mi).cellID), indPers);
    spikeRateNaive{mi,1} = act.naiveLearner(mi).spikeRates(indPers);
    spikeRateNaive{mi,2} = act.naiveLearner(mi).spikeRates(indTrans);
    
    % expert
    idPers = setdiff(match{mi,2}, 0);
    indPers = find(ismember(act.expert(mi).cellID, idPers));
    indTrans = setdiff(1:length(act.expert(mi).cellID), indPers);
    spikeRateExpert{mi,1} = act.expert(mi).spikeRates(indPers);
    spikeRateExpert{mi,2} = act.expert(mi).spikeRates(indTrans);
end

figure, hold on
naivePersSR = cellfun(@mean, spikeRateNaive(:,1));
naiveTransSR = cellfun(@mean, spikeRateNaive(:,2));
expertPersSR = cellfun(@mean, spikeRateExpert(:,1));
expertTransSR = cellfun(@mean, spikeRateExpert(:,2));
p = zeros(1,2); m = cell(1,2);
[~, p(1), m{1}] = paired_test(naivePersSR, naiveTransSR);
[~, p(2), m{2}] = paired_test(expertPersSR, expertTransSR);
barOffset = 0.2;
barWidth = 0.4;
bar(1-barOffset, mean(naivePersSR), barWidth, 'facecolor', 'w', 'edgecolor', colorsPersistent(1,:))
errorbar(1-barOffset, mean(naivePersSR), sem(naivePersSR), 'k')
bar(1+barOffset, mean(naiveTransSR), barWidth, 'facecolor', colorsTransient(1,:), 'edgecolor', colorsTransient(1,:))
errorbar(1+barOffset, mean(naiveTransSR), sem(naiveTransSR), 'k')

bar(2-barOffset, mean(expertPersSR), barWidth, 'facecolor', 'w', 'edgecolor', colorsPersistent(2,:))
errorbar(2-barOffset, mean(expertPersSR), sem(expertPersSR), 'k')
bar(2+barOffset, mean(expertTransSR), barWidth, 'facecolor', colorsTransient(2,:), 'edgecolor', colorsTransient(2,:))
errorbar(2+barOffset, mean(expertTransSR), sem(expertTransSR), 'k')

xlim([0.5 2.5])
xticks([1,2]), xticklabels({'Naive', 'Expert'})
ylabel('Rate of inferred spike (Hz)')
title(sprintf('p = %s (m = %s),   p = %s (m = %s)', num2str(p(1),3), m{1}, num2str(p(2),3), m{2}))


%% Cumulative proportion
histRange = 0:0.02:3;
cumNaivePers = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'cdf'), spikeRateNaive(:,1), 'un', 0));
cumNaiveTrans = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'cdf'), spikeRateNaive(:,2), 'un', 0));
cumExpertPers = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'cdf'), spikeRateExpert(:,1), 'un', 0));
cumExpertTrans = cell2mat(cellfun(@(x) histcounts(x, histRange, 'norm', 'cdf'), spikeRateExpert(:,2), 'un', 0));

figure, hold on
boundedline(histRange(2:end), mean(cumNaivePers), sem(cumNaivePers), 'cmap', colorsPersistent(1,:))
boundedline(histRange(2:end), mean(cumNaiveTrans), sem(cumNaiveTrans), 'cmap', colorsTransient(1,:))
boundedline(histRange(2:end), mean(cumExpertPers), sem(cumExpertPers), 'cmap', colorsPersistent(2,:))
boundedline(histRange(2:end), mean(cumExpertTrans), sem(cumExpertTrans), 'cmap', colorsTransient(2,:))
xlabel('Rate of inferred spikes (Hz)')
ylabel('Cumulative proportion')
set(gca, 'fontName', 'Arial', 'fontsize', 15, 'box', 'off')






%%
%% 6D - Proportion angle-tuned
%%
matchFn = 'cellMatching_beforeNafter.mat';
angleTuningFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat';

load(sprintf('%s%s',calciumDir, matchFn), 'match');
tune = load(sprintf('%s%s',calciumDir, angleTuningFn), 'expert', 'naive');

learnerInd = [1:4,7,9];
tune.naive = tune.naive(learnerInd);
numMice = length(learnerInd);


propTuned = zeros(numMice,2); % 1: naive 2: expert
for mi = 1 : numMice
    if mi == 2
        tempInd = find(tune.naive(mi).touchID > 5000);
        propTuned(mi,1) = sum(tune.naive(mi).tuned(tempInd)) / length(match{mi,1});
    else
        propTuned(mi,1) = sum(tune.naive(mi).tuned) / length(match{mi,1});
    end
    propTuned(mi,2) = sum(tune.expert(mi).tuned) / length(match{mi,3});
end

figure, hold on
for i = 1 : size(propTuned,1)
    plot(propTuned(i,:), 'ko-')
end

[~,p] = ttest(propTuned(:,1) - propTuned(:,2));

errorbar(mean(propTuned), sem(propTuned), 'ro', 'lines', 'no')
xticks([1:2])
xticklabels({'Naive', 'Expert'})
xtickangle(45)
xlim([0.5 2.5])
ylim([0 0.45]), yticks(0:0.1:0.3)
ylabel('Proportion (/touch cells)')
title({'Tuned cells'; sprintf('p = %f', p)})







%%
%% 6E and Extra
%%
%% 2020/10/26 Comparing sorting between naive and expert
% Matched neurons change in angle tuning
% Follow-up from idPersTuned_naive, because this means these neurons are
% tuned in both sessions. Find in each 'tune' to ensure they are correctly
% matched between sessions.

stretchFactor = 50;
naiveTC = [];
tunedAngleNaive = [];
sharpnessNaive = [];
expertTC = [];
for mi = 1 : numMice
    idNaivePers = match{mi,1}(find(match{mi,2}));
    idExpertPers = match{mi,2}(find(match{mi,2}));
    
    % Starting from expert
    idExpertPersTunedExpert = intersect(idExpertPers, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    indAllNaivePersTunedExpert = find(ismember(match{mi,2}, idExpertPersTunedExpert)); % indAll for index in 'match' / ind is for index in 'tune'
    idNaivePersTunedExpert = match{mi,1}(indAllNaivePersTunedExpert);
    
    % back to naive
    idNaivePersTunedNaive = intersect(idNaivePers, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    
    % ID's for naive
    idPersTuned_naive = intersect(idNaivePersTunedNaive, idNaivePersTunedExpert); % for the last process, _naive

    for ci = 1 : length(idPersTuned_naive)
        naiveTC = [naiveTC; cellfun(@mean, tune.naive(mi).val{find(tune.naive(mi).touchID == idPersTuned_naive(ci))})'];
        indAllPersTuned_naive = find(match{mi,1} == idPersTuned_naive(ci));
        idPersTuned_expert = match{mi,2}(indAllPersTuned_naive);
        expertTC = [expertTC; cellfun(@mean, tune.expert(mi).val{find(tune.expert(mi).touchID == idPersTuned_expert)})'];
    end
    tunedAngleNaive = [tunedAngleNaive; tune.naive(mi).tunedAngle(find(ismember(tune.naive(mi).touchID, idPersTuned_naive)))];
    sharpnessNaive = [sharpnessNaive; tune.naive(mi).sharpness(find(ismember(tune.naive(mi).touchID, idPersTuned_naive)))];
end

[~, indsortSharpness] = sort(sharpnessNaive, 'descend');
tempTunedAngle = tunedAngleNaive(indsortSharpness);
[~, indsortAngle] = sort(tempTunedAngle);
indsort = indsortSharpness(indsortAngle);

tilingMapNaive = zeros(length(tunedAngleNaive), 7 * stretchFactor);
tilingMapExpert = zeros(length(tunedAngleNaive), 7 * stretchFactor);
for i = 1 : length(tunedAngleNaive)
    tempMap = naiveTC(i,:);
    normTempMap = min_max_normalization(tempMap);
    for j = 1 : 7
        tilingMapNaive(i,(j-1)*stretchFactor+1:j*stretchFactor) = deal(normTempMap(j));
    end
    tempMap = expertTC(i,:);
    normTempMap = min_max_normalization(tempMap);
    for j = 1 : 7
        tilingMapExpert(i,(j-1)*stretchFactor+1:j*stretchFactor) = deal(normTempMap(j));
    end
end
%%
tempMapNaive = tilingMapNaive(:,[1,stretchFactor*1+1, stretchFactor*2+1, stretchFactor*3+1, stretchFactor*4+1, stretchFactor*5+1, stretchFactor*6+1]);
inds = [1:size(tempMapNaive,1)]';
numGroup = size(tempMapNaive,2);
sortiCell = cell(numGroup,1);
for i = 1 : numGroup
    row2sort = find(tempMapNaive(:,i)==1); % since it's min-max normalized.
    if isempty(row2sort)
        sortiCell{i} = [];
    else
        col2sort = setdiff(1:numGroup,i);
        tempsorti = nested_sorting(tempMapNaive(row2sort,col2sort),inds(row2sort));
        if size(tempsorti,2) ~= 1
            tempsorti = tempsorti';
        end
        sortiCell{i} = tempsorti;
    end
end
sorti = cell2mat(sortiCell);
    
newMapNaive = tilingMapNaive(sorti,:);
newMapExpertNaive = tilingMapExpert(sorti,:);

tempMapExpert = tilingMapExpert(:,[1,stretchFactor*1+1, stretchFactor*2+1, stretchFactor*3+1, stretchFactor*4+1, stretchFactor*5+1, stretchFactor*6+1]);
inds = [1:size(tempMapExpert,1)]';
numGroup = size(tempMapExpert,2);
sortiCell = cell(numGroup,1);
for i = 1 : numGroup
    row2sort = find(tempMapExpert(:,i)==1); % since it's min-max normalized.
    if isempty(row2sort)
        sortiCell{i} = [];
    else
        col2sort = setdiff(1:numGroup,i);
        tempsorti = nested_sorting(tempMapExpert(row2sort,col2sort),inds(row2sort));
        if size(tempsorti,2) ~= 1
            tempsorti = tempsorti';
        end
        sortiCell{i} = tempsorti;
    end
end
sorti = cell2mat(sortiCell);
    
newMapExpert = tilingMapExpert(sorti,:);
%%
%% 6E
%%
figure
subplot(121)
imshow(newMapNaive), hold on, 
axis tight, title('Naive')


subplot(122)
imshow(newMapExpert), hold on, 
axis tight, title('Expert')

%% Extra - Expert sorted by Naive
figure
imshow(newMapExpertNaive), hold on, 
axis tight, title('Expert - sorted by Naive')







%%
%% S7E - Naive VS Expert angle tuning curve, in L2/3 and L4
%% Extra - compare between L2/3 and L4 in Naive and Expert
%%
% load data and run calculation
tune = load([calciumDir, 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat']);
touch = load([calciumDir, 'glmResults_devExp_touch_NC'], 'naive', 'expert');
mice = [25,27,30,36,39,52];
sessions = {[4,19], [3,10], [3,21], [1,17], [1,23], [3,21]};
learnerInd = [1,2,3,4,7,9];
tune.naive = tune.naive(learnerInd);
depthThresh = 350;
angles = 45:15:135;
tuningNaive = zeros(length(learnerInd), length(angles), 2); % (:,:,1) - L2/3, (:,:,2) - L4
tuningExpert = zeros(length(learnerInd), length(angles), 2); % (:,:,1) - L2/3, (:,:,2) - L4
numL4 = zeros(6,2); % (:,1) naive, (:,2) expert
propL4 = zeros(6,3,2); % (:,1,:) touch prop, (:,2,:) tuned prop, (:,3,:) not-tuned prop
absnumL4 = zeros(6,3,2); 
for i = 1 : length(tune.naive)
    temp = tune.naive(i);
    L23ind = find(temp.depth < depthThresh);
    L4ind = find(temp.depth >= depthThresh);
    for j = 1 : length(angles)
        tuningNaive(i,j,1) = length(find(temp.tunedAngle(L23ind) == angles(j))) / sum(temp.tuned(L23ind));
        tuningNaive(i,j,2) = length(find(temp.tunedAngle(L4ind) == angles(j))) / sum(temp.tuned(L4ind));
    end
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumDir,mice(i),mice(i),sessions{i}(1)));
    cellID = u.cellNums;
    cellDepths = u.cellDepths;
    numL4(i,1) = length(find(cellDepths >= depthThresh));
    absnumL4(i,1,1) = length(L4ind);
    absnumL4(i,2,1) = sum(temp.tuned(L4ind));
    absnumL4(i,3,1) = absnumL4(i,1,1) - absnumL4(i,2,1);
    propL4(i,1,1) = absnumL4(i,1,1) / numL4(i,1);
    propL4(i,2,1) = absnumL4(i,2,1) / absnumL4(i,1,1);
    propL4(i,3,1) = absnumL4(i,3,1) / absnumL4(i,1,1);
    
    temp = tune.expert(i);
    L23ind = find(temp.depth < depthThresh);
    L4ind = find(temp.depth >= depthThresh);
    for j = 1 : length(angles)
        tuningExpert(i,j,1) = length(find(temp.tunedAngle(L23ind) == angles(j))) / sum(temp.tuned(L23ind));
        tuningExpert(i,j,2) = length(find(temp.tunedAngle(L4ind) == angles(j))) / sum(temp.tuned(L4ind));
    end
    
    load(sprintf('%s%03d\\UberJK%03dS%02d_NC',calciumDir,mice(i),mice(i),sessions{i}(2)));
    cellDepths = u.cellDepths;
    numL4(i,2) = length(find(cellDepths >= depthThresh));
    absnumL4(i,1,2) = length(L4ind);
    absnumL4(i,2,2) = sum(temp.tuned(L4ind));
    absnumL4(i,3,2) = absnumL4(i,1,2) - absnumL4(i,2,2);
    propL4(i,1,2) = absnumL4(i,1,2) / numL4(i,2);
    propL4(i,2,2) = absnumL4(i,2,2) / absnumL4(i,1,2);
    propL4(i,3,2) = absnumL4(i,3,2) / absnumL4(i,1,2);
end


%%
%% Extra - Naive L2/3 vs L4
%%
validInd = 1:6;
naiveL23 = squeeze(tuningNaive(validInd,:,1));
naiveL4 = squeeze(tuningNaive(validInd,:,2));
expertL23 = squeeze(tuningExpert(validInd,:,1));
expertL4 = squeeze(tuningExpert(validInd,:,2));


figure('units', 'normalized', 'position', [0.1 0.1 0.5 0.3])
subplot(1,2,1), hold on
errorbar(angles, mean(naiveL23), sem(naiveL23), 'k*-')
errorbar(angles, mean(naiveL4), sem(naiveL4), 'ko--')
legend({'L2/3', 'L4'})
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([0:0.1:0.5]), ylim([0 0.5]), xlim([40 140])
ylabel('Proportion (/all tuned cells)')
title('Naive L2/3 vs L4')

subplot(1,2,2), hold on
naiveDiff = naiveL23 - naiveL4;
errorbar(angles, mean(naiveDiff), sem(naiveDiff), 'k-')
plot([40 140], [0 0], '--', 'color', [0.6 0.6 0.6])
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([-0.2:0.1:0.2]), ylim([-0.25 0.2]), xlim([40 140])
ylabel('Diff in proportion')
title('Naive L2/3 - L4')


% Bonferroni-Holm correction of multiple t-test
alpha = 0.05;
pvals = zeros(7,1);
for i = 1 : 7
    [~, pvals(i)] = ttest(naiveDiff(:,i));
end
[pvalsorted, sorti] = sort(pvals);
passed = [];
for i = 1 : 7
    if pvalsorted(i) < alpha/(7-i+1)
        passed = [passed, sorti(i)];
    else
        break
    end
end
passed



%%
%% Extra - Expert L2/3 vs L4
%%
figure('units', 'normalized', 'position', [0.15 0.15 0.5 0.3])
subplot(1,2,1), hold on
errorbar(angles, mean(expertL23), sem(expertL23), 'k*-')
errorbar(angles, mean(expertL4), sem(expertL4), 'ko--')
legend({'L2/3', 'L4'})
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([0:0.1:0.5]), ylim([0 0.5]), xlim([40 140])
ylabel('Proportion (/all tuned cells)')
title('Expert L2/3 vs L4')

subplot(1,2,2), hold on
expertDiff = expertL23 - expertL4;
errorbar(angles, mean(expertDiff), sem(expertDiff), 'k-')
plot([40 140], [0 0], '--', 'color', [0.6 0.6 0.6])
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([-0.2:0.1:0.2]), ylim([-0.25 0.2]), xlim([40 140])
ylabel('Diff in proportion')
title('Expert L2/3 - L4')

% Bonferroni-Holm correction of multiple t-test
alpha = 0.05;
pvals = zeros(7,1);
for i = 1 : 7
    [~, pvals(i)] = ttest(expertDiff(:,i));
end
[pvalsorted, sorti] = sort(pvals);
passed = [];
for i = 1 : 7
    if pvalsorted(i) < alpha/(7-i+1)
        passed = [passed, sorti(i)];
    else
        break
    end
end
passed


%%
%% S7E Left - L2/3 Naive VS Expert
%%
figure('units', 'normalized', 'position', [0.2 0.2 0.5 0.3])
subplot(1,2,1), hold on
errorbar(angles, mean(naiveL23), sem(naiveL23), 'k*-')
errorbar(angles, mean(expertL23), sem(expertL23), 'ko--')
legend({'Naive', 'Expert'})
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([0:0.1:0.5]), ylim([0 0.5]), xlim([40 140])
ylabel('Proportion (/all tuned cells)')
title('L2/3 naive vs expert')

subplot(1,2,2), hold on
L23Diff = naiveL23 - expertL23;
errorbar(angles, mean(L23Diff), sem(L23Diff), 'k-')
plot([40 140], [0 0], '--', 'color', [0.6 0.6 0.6])
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([-0.2:0.1:0.2]), ylim([-0.25 0.2]), xlim([40 140])
ylabel('Diff in proportion')
title('L2/3 naive - expert')

% Bonferroni-Holm correction of multiple t-test
alpha = 0.05;
pvals = zeros(7,1);
for i = 1 : 7
    [~, pvals(i)] = ttest(L23Diff(:,i));
end
[pvalsorted, sorti] = sort(pvals);
passed = [];
for i = 1 : 7
    if pvalsorted(i) < alpha/(7-i+1)
        passed = [passed, sorti(i)];
    else
        break
    end
end
passed


%%
%% S7E Right - L4 Naive VS Expert
%%
figure('units', 'normalized', 'position', [0.25 0.25 0.5 0.3])
subplot(1,2,1), hold on
errorbar(angles, mean(naiveL4), sem(naiveL4), 'k*-')
errorbar(angles, mean(expertL4), sem(expertL4), 'ko--')
legend({'Naive', 'Expert'})
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([0:0.1:0.5]), ylim([0 0.5]), xlim([40 140])
ylabel('Proportion (/all tuned cells)')
title('L4 naive vs expert')

subplot(1,2,2), hold on
L4Diff = naiveL4 - expertL4;
errorbar(angles, mean(L4Diff), sem(L4Diff), 'k-')
plot([40 140], [0 0], '--', 'color', [0.6 0.6 0.6])
xlabel('Object angle (\circ)')
xticks([45:15:135]), yticks([-0.2:0.1:0.2]), ylim([-0.25 0.2]), xlim([40 140])
ylabel('Diff in proportion')
title('L4 naive - expert')

% Bonferroni-Holm correction of multiple t-test
alpha = 0.05;
pvals = zeros(7,1);
for i = 1 : 7
    [~, pvals(i)] = ttest(L4Diff(:,i));
end
[pvalsorted, sorti] = sort(pvals);
passed = [];
for i = 1 : 7
    if pvalsorted(i) < alpha/(7-i+1)
        passed = [passed, sorti(i)];
    else
        break
    end
end
passed





%% Figure 6F
% persistent neurons angle-tunind distribution
%% basic setting

matchFn = 'cellMatching_beforeNafter.mat';
angleTuningFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat';

load(sprintf('%s%s',calciumDir, matchFn), 'match');
tune = load(sprintf('%s%s',calciumDir, angleTuningFn), 'expert', 'naive');


learnerInd = [1,2,3,4,7,9];
numMice = size(match,1);

% cleaning up nonlearners and upper layer of JK027
tune.naive = tune.naive(learnerInd);
jk027inds = find(tune.naive(2).touchID > 5000);

fn = fieldnames(tune.naive);
for i = 1 : length(fn)
    tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(jk027inds);
end

angles = 45:15:135;

colorsTransient = [248 171 66; 40 170 225] / 255;
colorsPersistent = [1 0 0; 0 0 1];

%% angle-tuning distribution in different groups
% divide into:
% (1) persistently angle-tuned neurons (can follow up across sessions) - Figure 6F
% (2) persistently active transiently angle-tuned neurons
% (3) transient & angle-tuned neurons - Figure S7F Left
% (4) persistent & angle-tuned (including (1) and (2))  - Figure S7F Right

persTuned = zeros(numMice, length(angles), 2); % (:,:,1) naive, (:,:,2) expert
persActiveTransTuned = zeros(numMice, length(angles), 2);
transTuned = zeros(numMice, length(angles), 2);
persAndTuned = zeros(numMice, length(angles), 2);

for mi = 1 : numMice
    % basic ID's
    idNaivePers = match{mi,1}(find(match{mi,2}));
    idExpertPers = match{mi,2}(find(match{mi,2}));
    idNaiveTrans = match{mi,1}(find(match{mi,2} == 0));
    idExpertTrans = setdiff(match{mi,3}, match{mi,2});
    
    % Starting from expert
    idExpertPersTunedExpert = intersect(idExpertPers, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    idExpertPersNottunedExpert = setdiff(idExpertPers, tune.expert(mi).touchID(find(tune.expert(mi).tuned))); % include both non-selective touch and non-touch neurons
    indAllNaivePersTunedExpert = find(ismember(match{mi,2}, idExpertPersTunedExpert)); % indAll for index in 'match' / ind is for index in 'tune'
    indAllNaivePersNottunedExpert = find(ismember(match{mi,2}, idExpertPersNottunedExpert));
    idNaivePersTunedExpert = match{mi,1}(indAllNaivePersTunedExpert);
    idNaivePersNottunedExpert = match{mi,1}(indAllNaivePersNottunedExpert);
    
    % back to naive
    idNaivePersTunedNaive = intersect(idNaivePers, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    idNaivePersNottunedNaive = setdiff(idNaivePers, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    
    % ID's for naive
    idPersTuned_naive = intersect(idNaivePersTunedNaive, idNaivePersTunedExpert); % for the last process, _naive
    idPersActiveTransTuned_naive = intersect(idNaivePersTunedNaive, idNaivePersNottunedExpert);
    idTransTuned_naive = intersect(idNaiveTrans, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    idPersAndTuned_naive = intersect(idNaivePers, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    
    % Indices for tune.naive(mi)
    indPersTuned_naive = find(ismember(tune.naive(mi).touchID, idPersTuned_naive));
    indPersActiveTransTuned_naive = find(ismember(tune.naive(mi).touchID, idPersActiveTransTuned_naive));
    indTransTuned_naive = find(ismember(tune.naive(mi).touchID, idTransTuned_naive));
    indPersAndTuned_naive = find(ismember(tune.naive(mi).touchID, idPersAndTuned_naive));
    % sanity check
    if ~all(ismember(indPersAndTuned_naive, union(indPersTuned_naive, indPersActiveTransTuned_naive)))
        error('Error in matching naive session')
    end
    
    % Indices for match{mi,1} (naive)
    indAllNaivePersTunedBoth = find(ismember(match{mi,1}, idPersTuned_naive));
    indAllNaivePersNottunedNaive = find(ismember(match{mi,1}, idNaivePersNottunedNaive));
    
    % Going back to expert (ID's for expert)
    idPersTuned_expert = match{mi,2}(indAllNaivePersTunedBoth);
    idPersActiveTransTuned_expert = intersect(match{mi,2}(indAllNaivePersNottunedNaive), idExpertPersTunedExpert);
    idTransTuned_expert = intersect(idExpertTrans, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    idPersAndTuned_expert = intersect(idExpertPers, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    
    % Indices for tune.expert(mi)
    indPersTuned_expert = find(ismember(tune.expert(mi).touchID, idPersTuned_expert));
    indPersActiveTransTuned_expert = find(ismember(tune.expert(mi).touchID, idPersActiveTransTuned_expert));
    indTransTuned_expert = find(ismember(tune.expert(mi).touchID, idTransTuned_expert));
    indPersAndTuned_expert = find(ismember(tune.expert(mi).touchID, idPersAndTuned_expert));
    % sanity check
    if ~all(ismember(indPersAndTuned_expert, union(indPersTuned_expert, indPersActiveTransTuned_expert)))
        error('Error in matching expert session')
    end
    
    for ai = 1 : length(angles)
        persTuned(mi,ai,1) = length(find(tune.naive(mi).tunedAngle(indPersTuned_naive) == angles(ai))) / length(indPersTuned_naive);
        persTuned(mi,ai,2) = length(find(tune.expert(mi).tunedAngle(indPersTuned_expert) == angles(ai))) / length(indPersTuned_expert);
        
        persActiveTransTuned(mi,ai,1) = length(find(tune.naive(mi).tunedAngle(indPersActiveTransTuned_naive) == angles(ai))) / length(indPersActiveTransTuned_naive);
        persActiveTransTuned(mi,ai,2) = length(find(tune.expert(mi).tunedAngle(indPersActiveTransTuned_expert) == angles(ai))) / length(indPersActiveTransTuned_expert);
        
        transTuned(mi,ai,1) = length(find(tune.naive(mi).tunedAngle(indTransTuned_naive) == angles(ai))) / length(indTransTuned_naive);
        transTuned(mi,ai,2) = length(find(tune.expert(mi).tunedAngle(indTransTuned_expert) == angles(ai))) / length(indTransTuned_expert);
        
        persAndTuned(mi,ai,1) = length(find(tune.naive(mi).tunedAngle(indPersAndTuned_naive) == angles(ai))) / length(indPersAndTuned_naive);
        persAndTuned(mi,ai,2) = length(find(tune.expert(mi).tunedAngle(indPersAndTuned_expert) == angles(ai))) / length(indPersAndTuned_expert);
        
    end
end

%%
%% 6F
%%
figure, hold on
for i = 1 : 2
    tempMat = squeeze(persTuned(:,:,i));
    errorbar(angles, mean(tempMat), sem(tempMat), 'o-', 'markersize', 10, 'color', colorsPersistent(i,:), 'capsize', 10)
end
legend({'Naive', 'Expert'}, 'box', 'off')
xlim([40 140]), xticks(angles), xlabel('Tuned angle(\circ)')
ylabel('Proportion (persistently tuned)')
ylim([0 0.4]), yticks(0:0.1:0.4)





%% Extra

figure, hold on
for i = 1 : 2
    tempMat = squeeze(persActiveTransTuned(:,:,i));
    errorbar(angles, mean(tempMat), sem(tempMat), 'o--', 'markersize', 10, 'color', colorsPersistent(i,:), 'capsize', 10)
end
legend({'Naive', 'Expert'}, 'box', 'off')
xlim([40 140]), xticks(angles), xlabel('Tuned angle(\circ)')
ylabel('Proportion (persistently active transiently tuned)')
ylim([0 0.4]), yticks(0:0.1:0.4)



%%
%% S7F left
%%
figure, hold on
for i = 1 : 2
    tempMat = squeeze(transTuned(:,:,i));
    errorbar(angles, mean(tempMat), sem(tempMat), 'o-', 'markersize', 10, 'color', colorsTransient(i,:), 'capsize', 10)
end
legend({'Naive', 'Expert'}, 'box', 'off')
xlim([40 140]), xticks(angles), xlabel('Tuned angle(\circ)')
ylabel('Proportion (transient & tuned)')
ylim([0 0.4]), yticks(0:0.1:0.4)


%%
%% S7F right
%%
figure, hold on
for i = 1 : 2
    tempMat = squeeze(persAndTuned(:,:,i));
    errorbar(angles, mean(tempMat), sem(tempMat), 'o-', 'markersize', 10, 'markerfacecolor', colorsPersistent(i,:), 'color', colorsPersistent(i,:), 'capsize', 10)
end
legend({'Naive', 'Expert'}, 'box', 'off')
xlim([40 140]), xticks(angles), xlabel('Tuned angle()')
ylabel('Proportion (persistent & tuned)')
ylim([0 0.4]), yticks(0:0.1:0.4)




%% Extra - Showing the difference only
% persistently tuned neurons and transient & tuned neurons
persTunedDiff = squeeze(persTuned(:,:,2)) - squeeze(persTuned(:,:,1));
transTunedDiff = squeeze(transTuned(:,:,2)) - squeeze(transTuned(:,:,1));

figure, hold on
errorbar(angles, mean(persTunedDiff), sem(persTunedDiff), 'o-', 'markersize', 10, 'color', colorsPersistent(2,:), 'capsize', 10);
errorbar(angles, mean(transTunedDiff), sem(transTunedDiff), 'o-', 'markersize', 10, 'color', colorsTransient(2,:), 'capsize', 10);
legend({'Persistently tuned', 'Transient & tuned'}, 'autoupdate', false, 'box', 'off')
plot([40 140], [0 0], '--', 'color', [0.6 0.6 0.6])
xlim([40 140]), xticks(angles), xlabel('Tuned angle ()')
ylim([-0.2 0.2]), yticks([-0.2:0.1:0.2])
ylabel('DProportion (Expert - Naive)')

%% statistical test
p = zeros(2, length(angles));
m = cell(2, length(angles));
for ai = 1 : length(angles)
    [~, p(1,ai), m{1,ai}] = paired_test(persTunedDiff(:,ai));
    [~, p(2,ai), m{2,ai}] = paired_test(transTunedDiff(:,ai));
end










%% Figure 6G
% Matched neurons change in angle tuning
% Follow-up from idPersTuned_naive, because this means these neurons are
% tuned in both sessions. Find in each 'tune' to ensure they are correctly
% matched between sessions.

matchFn = 'cellMatching_beforeNafter.mat';
angleTuningFn = 'angle_tuning_summary_preAnswer_perTouch_NC_PTC.mat';

load(sprintf('%s%s',calciumDir, matchFn), 'match');
tune = load(sprintf('%s%s',calciumDir, angleTuningFn), 'expert', 'naive');


learnerInd = [1,2,3,4,7,9];
numMice = size(match,1);

% cleaning up nonlearners and upper layer of JK027
tune.naive = tune.naive(learnerInd);
jk027inds = find(tune.naive(2).touchID > 5000);

fn = fieldnames(tune.naive);
for i = 1 : length(fn)
    tune.naive(2).(fn{i}) = tune.naive(2).(fn{i})(jk027inds);
end

angles = 45:15:135;


histRange = 0:15:105;
matchedAngleDiff = cell(numMice, 1);
angleDiffDist = zeros(numMice, length(histRange)-1);
angleDiffDistShuffle = zeros(numMice, length(histRange)-1);
numShuffle = 1000;
for mi = 1 : numMice
    idNaivePers = match{mi,1}(find(match{mi,2}));
    idExpertPers = match{mi,2}(find(match{mi,2}));
    
    % Starting from expert
    idExpertPersTunedExpert = intersect(idExpertPers, tune.expert(mi).touchID(find(tune.expert(mi).tuned)));
    indAllNaivePersTunedExpert = find(ismember(match{mi,2}, idExpertPersTunedExpert)); % indAll for index in 'match' / ind is for index in 'tune'
    idNaivePersTunedExpert = match{mi,1}(indAllNaivePersTunedExpert);
    
    % back to naive
    idNaivePersTunedNaive = intersect(idNaivePers, tune.naive(mi).touchID(find(tune.naive(mi).tuned)));
    
    % ID's for naive
    idPersTuned_naive = intersect(idNaivePersTunedNaive, idNaivePersTunedExpert); % for the last process, _naive
    
    matchedAngleDiff{mi} = zeros(length(idPersTuned_naive),1);
    for ci = 1 : length(idPersTuned_naive)
        tunedAngleNaive = tune.naive(mi).tunedAngle(find(tune.naive(mi).touchID == idPersTuned_naive(ci)));
        indAllPersTuned_naive = find(match{mi,1} == idPersTuned_naive(ci));
        idPersTuned_expert = match{mi,2}(indAllPersTuned_naive);
        tunedAngleExpert = tune.expert(mi).tunedAngle(find(tune.expert(mi).touchID == idPersTuned_expert));
        matchedAngleDiff{mi}(ci) = abs(tunedAngleExpert - tunedAngleNaive);
    end

    angleDiffDist(mi,:) = histcounts(matchedAngleDiff{mi}, histRange, 'norm', 'prob');
    
    % For shuffled data
    tempDist = zeros(numShuffle,length(histRange)-1);
    indPersTuned_naive = find(ismember(tune.naive(mi).touchID, idPersTuned_naive)); % indices of persistently tuned neurons in naive session, from touchID
    indAllNaive = find(ismember(match{mi,1}, idPersTuned_naive)); % indices of persistently tuned neurons in naive session, from ALL NEURONS (from match{mi,1})
    idAllExpert = match{mi,2}(indAllNaive); % indices of persistently tuned neurons in expert session, from ALL NEURONS (from match{mi,2})
    % sanity check
    if ~isempty(find(idAllExpert == 0))
        error('Error matching tuned neurons in expert sessions')
    end
    indPersTuned_expert = find(ismember(tune.expert(mi).touchID, idAllExpert)); % indices of persistently tuned neurons in expert session, from touchID
    tunedAngleNaiveAll = tune.naive(mi).tunedAngle(indPersTuned_naive);
    tunedAngleExpertAll = tune.expert(mi).tunedAngle(indPersTuned_expert);
    for shi = 1 : numShuffle
        tunedAngleExpertShuffle = tunedAngleExpertAll(randperm(length(tunedAngleExpertAll)));
        tempDiff = abs(tunedAngleExpertShuffle - tunedAngleNaiveAll);
        tempDist(shi,:) = histcounts(tempDiff, histRange, 'norm', 'prob');
    end
    angleDiffDistShuffle(mi,:) = mean(tempDist);
end

figure, hold on
bar(histRange(1:2), mean(angleDiffDist(:,1:2)), 'facecolor', 'k')
bar(histRange(1:2), mean(angleDiffDistShuffle(:,1:2)), 'facecolor', [0.6 0.6 0.6])
legend({'Data', 'Shuffle'}, 'autoupdate', 'off')
bar(histRange(3:7), mean(angleDiffDistShuffle(:,3:7)), 'facecolor', [0.6 0.6 0.6])
bar(histRange(3:7), mean(angleDiffDist(:,3:7)), 'facecolor', 'k')

errorbar(histRange(1:end-1), mean(angleDiffDist), sem(angleDiffDist), 'k', 'lines', 'no')
errorbar(histRange(1:end-1), mean(angleDiffDistShuffle), sem(angleDiffDistShuffle), 'color', [0.6 0.6 0.6], 'lines', 'no')
xticks(histRange(1:end-1))
xlabel('\DeltaTuned angle (\circ)')
ylabel('Proportion')






