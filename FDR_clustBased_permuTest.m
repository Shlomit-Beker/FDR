% This function gets time(rows)*frequency(columns) maps as inputs, calcultes z-scored p values, 
% and correct for False Discovery Rate (FDR) using cluster-based permutation test (Mike X Cohen,
% 2014: Analyzing Neural Time Series Data (MIT Press)). 
% No channel dinension for this version
% Inputs: MAT1 MAT2 - participants (rows)*time (columns)*frequency (3rd dim) matrices (could be ITPC or TFR)
%           timeoi - vector with times-of-interest 
%           freqoi - vector with frequencies-of-interest 
% Outputs: cluster_thresh - threshold for cluster size (to which real
%           clusters are compared); figures 
%
%  Shlomit Beker 2021 <shlomitbeker@gmail.com>

function cluster_thresh = FDR_clustBased_permuTest(MAT1,MAT2,timeoi,freqoi)

N = 10000; % number of permutation - defauls
pval = 0.05; % p value - default

% convert p-value to Z value
zval = abs(norminv(pval));
timeReduction = [100:length(timeoi)-100]; 
freqReduction = [1:length(freqoi)];

numTime = length(timeReduction);
numFreq = length(freqReduction);
clim = [0 50];                                                        % color scale
map1 = MAT1(:,freqReduction,timeReduction); 
map2 = MAT2(:,freqReduction,timeReduction); 

realDiffItc = squeeze(mean(map1,1) - mean(map2,1)); % the real difference in ITC between group1 and group1
ITCall = cat(1, map1, map2);

nASD = size(map2,1);                                             % number of subjects in the 1st group
nTD = size(map1,1);                                                 % number of subjects in the 2nd group
randItcDiff = zeros(N,numFreq,numTime);

% loop through permutations

% randomly assign subjects to Group1 or Group1 
for n = 1:N
    randSubj = randperm(size(ITCall,1));                       
    Group1_rand = ITCall(randSubj(1:nTD),:,:);
    Group2_rand = ITCall(randSubj(nASD+1:end),:,:);
    randItcDiff(n,:,:) = squeeze(mean(Group1_rand,1) - mean(Group2_rand,1));
end

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(randItcDiff));
std_h0  = squeeze(std(randItcDiff));

% now threshold real data with Z-score
zmap = (realDiffItc-mean_h0) ./ std_h0;
zmap(isnan(zmap)==1) = 0; 

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;

%plot
figure(3), clf
subplot(221)
imagesc(timeoi(timeReduction),freqoi(freqReduction),realDiffItc);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
title('TF map of real power values')

subplot(222)
imagesc(timeoi(timeReduction),freqoi(freqReduction),realDiffItc);
hold on
contour(timeoi(timeReduction),freqoi(freqReduction),logical(zmap),1,'linecolor','k');
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
title('Power values and outlined significance regions')

subplot(223)
imagesc(timeoi(timeReduction),freqoi(freqReduction),zmap);
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
title('Thresholded TF map of Z-values')

%%

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,N);

% for maximum-pixel based correction
max_val = zeros(N,2); % "2" for min/max

for permi = 1:N
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(randItcDiff(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    threshimg(isnan(threshimg)==1)=0;
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
 
    % get extreme values (smallest and largest)
    temp = sort( reshape(randItcDiff(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end

%% show histograph of maximum cluster sizes

figure(3), clf
hist(max_cluster_sizes,20);
xlabel('Maximum cluster sizes'), ylabel('Number of observations')
title('Expected cluster sizes under the null hypothesis')


% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

% plot tresholded results
figure(4), clf
subplot(221)
imagesc(timeoi(timeReduction),freqoi(freqReduction),realDiffItc)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power, no thresholding') 
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(222)
imagesc(timeoi(timeReduction),freqoi(freqReduction),realDiffItc)
hold on
contour(timeoi(timeReduction),freqoi(freqReduction),logical(zmap),1,'linecolor','k')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF power with contour')
set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')

subplot(223)
imagesc(timeoi(timeReduction),freqoi(freqReduction),zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded')
set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal')




