% This function correct TimeFreq maps for multiple comparisons using the Guthrie-Buchwald (1991) 
%
% Inputs: map1, map2 - time(rows)*frequency(column) matrices (could be ITPC or Power)
%           timeoi - vector with times-of-interest 
%           freqoi - vector with frequencies-of-interest 
%           GBnum - critical number of consecutive time points (see table 1
%           in Guthrie and Buchwald (1991) 
%
% Outputs: cluster_thresh - threshold for cluster size (to which real
%           clusters are compared); figures.
%
%  Shlomit Beker 2021 <shlomitbeker@gmail.com>


function GBmap = FDR_GB(map1,map2,timeoi,freqoi,timeLine,timeWindow,GBnum)
% calculate zmap

N = 10000; % number of permutation - defauls
pval = 0.05;

% convert p-value to Z value
zval = abs(norminv(pval));
timeReduction = [100:length(timeoi)-100]; 
freqReduction = [1:length(freqoi)];

numTime = length(timeReduction);
numFreq = length(freqReduction);
clim = [0 0.3];                                                           % color scale
map1 = tfrTD(:,freqReduction,timeReduction); 
map2 = tfrASD(:,freqReduction,timeReduction); 

realDiffItc = squeeze(mean(map1,1) - mean(map2,1)); % the real difference in ITC between group1 and group1
ITCall = cat(1, map1, map2);

nTD = size(map1,1); 
nASD = size(map2,1);                                             % number of subjects in the 1st group
                                                % number of subjects in the 2nd group
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


%% define the areas of interest

logMap = logical(zmap);                                                 %logical 0/1 map for the whole array
timeLine = [259;426;592;759];                                     
timeWindow = 150;
sampWindow = floor(timeWindow/1000*256);                    %add +- ms in samples
timeWin = [timeLine-sampWindow,timeLine+sampWindow];
GBnum = 8;

%% This loops check N consecutive values in the logical map

clear column_inds
sigInds = [];

for m = 1:length(timeWin)
    logMapWin = logMap;
    logRows = 1:size(logMap,2);
    z = find(logRows<timeWin(m,1)|logRows>timeWin(m,2));
    logMapWin(:,z) = 0;
    [row,column] = find(logMapWin);
    rows_sig = unique(row);
    for i = 1:length(rows_sig)
        row_check_consec = find(row==rows_sig(i));
        column_check_consec = column(row_check_consec);
        check_consec = diff(column_check_consec);
        gaps = find(check_consec>1);
        if isempty(gaps) && length(check_consec)>GBnum
            sigInds = cat(2,sigInds,[rows_sig(i);min(column_check_consec); max(column_check_consec)]);
        else
            j = [1,gaps'+1,length(column_check_consec)];
            long_consec = find(diff(j)>GBnum);
            for k= long_consec
                sigInds = cat(2,sigInds,[rows_sig(i);column_check_consec(j(k));column_check_consec(j(k+1)-1)]);
            end
        end
    end
end

%% Plot the Guthrie-Buchwald significance Map

GBmap = zeros(size(logMap));     
for j = 1:length(sigInds)
            GBmap(sigInds(1,j),sigInds(2,j):sigInds(3,j))=1;
end
%figure; imagesc(GBmap)
%set(gca,'xlim',xlim,'ydir','norm')


%
Nint = 500;
[x,y] = meshgrid(Timeoi,freqoi); % low-res grid
[x2,y2] = meshgrid(Timeoi(1):1/Nint/5:Timeoi(end),freqoi(1):.01:freqoi(end));  %high-res grid
dataInterp_logic = interp2(x,y,GBmap, x2,y2, 'linear'); %interpolate up 

% plot significat areas as white contours

pSig_GB =dataInterp_logic>0;

xx=Timeoi(1):1/Nint/5:Timeoi(end);
yy=[freqoi(1):.01:freqoi(end)]';
hold on;
contour3(xx,yy,pSig_GB+100,'k', 'LineWidth',1); %the +100 is to plot the contour on top of the highest surf value


