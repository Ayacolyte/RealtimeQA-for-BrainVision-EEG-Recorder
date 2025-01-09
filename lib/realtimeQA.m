function [data,OverallBadMask,NoSignalMask,HighAmpliMask,FrequencyNoiseMask,LowCorrelationMask,NS_track_New] = realtimeQA(data,props,fig,HighPassband,selechanns,badWindowThreshold,...
          robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,NS_track,chann_names)
% realtime_QA - Performs real time quality assessments on input EEG data from RDA (Remote Data Access) protocal of BrainVision Recorder.
% 
% This work is adapted from previous work done by Li Dong of UESTC (DRT) 
%
% Author: Changyuan Ni
% Contact: changyuan.d.ni@outlook.com
% Affiliation: ShanghaiTech University
%
% Date of Creation: January 5th, 2024
%
% Notes: This function was written on MATLAB R2023a, update: tested on R2024b.
% Needs Image Proccessing Toobox and Statistics and Machine Learning Toolbox
%% USAGE %%
% To test the output of this work with pre-existing data, run test.m
% To deploy the 
%%---------------------------Version History----------------------------%%
% 2024.01.05 - Testing using sliding 10s window. Conclusion:insufficient data to produce meaningful QA
% 2024.01.09 - Try using by-second analysis plus cumulative documentation of noisiness and deviation
% 2024.01.11 - Using persistant to declare local static variables that hold its value between function calls.
% 2024.01.31 - Written and tested the UI generation script createLamps
% 2024.02.02 - Visualisation completed, however filter is not robust enough for 1-window-input
% 2024.12.31 - Adds dynamic UI positioning for varying electrode number
% 2025.01.02 - Optimise blank space in UI
% 2025.01.05 - Adds a global scaler, compensating EEGlab's built-in scaling (not seen in RDA)
% 2025.01.07 - Adds a column displaying bad channel causes, and a legend
% 2025.01.08 - Adds Always-on-top feature
%
%


% handle insufficient inputs
% In Jan 9th this was nargin < 3, being data, props and channelDeviation_cuml
%persistent channelDeviation_cuml noisiness_cuml 


if nargin < 3
    error ('Three inputs (data,props and fog) are reqiured at least!!!!!');
elseif nargin == 3
    % If only two inputs are present, initiate a GUI to prompt the user for the
    % parameters needed. 

    geometry = {[3, 1], [3, 1], [3, 1], [3, 1], [3, 1], [3, 1], [3, 1], 1};
    geomvert = [1 1 1 1 1 1 1 1];

    uilist = {{'style', 'text', 'string', 'High Pass Band (Hz)'} ...
              {'style', 'edit', 'string', '1'} ...
              {'style', 'text', 'string', 'Selected Channals ("all" if want to select all)'} ...
              {'style', 'edit', 'string', 'all'} ...
              {'style', 'text', 'string', 'Bad Window Threshold'} ...
              {'style', 'edit', 'string', '0.6'} ...
              {'style', 'text', 'string', 'Robust Deviation Threshold'} ...
              {'style', 'edit', 'string', '5'} ...
              {'style', 'text', 'string', 'Power Frequency (Hz)'} ...
              {'style', 'edit', 'string', '50'} ...
              {'style', 'text', 'string', 'Frequency Noise Threshold'} ...
              {'style', 'edit', 'string', '3'} ...
              {'style', 'text', 'string', 'Correlation Threshold'} ...
              {'style', 'edit', 'string', '0.6'} ...
              {'style', 'checkbox', 'string', 'Notch filter the data instead of pass band', 'value', 0} ...
              };

        result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, 'title', 'Realtime_QA', 'helpcom', 'pophelp(''pop_eegfiltnew'')');

    if isempty(result), return; end

    HighPassband = str2double(result{1});
    selechanns = result{2};
    badWindowThreshold = str2double(result{3});
    robustDeviationThreshold = str2double(result{4});
    PowerFrequency = str2double(result{5});
    FrequencyNoiseThreshold = str2double(result{6});
    correlationThreshold = str2double(result{7});
    flagNotchFilter = str2double(result{8});
    disp(result{2})
elseif nargin == 4
    selechanns = 'all';      
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
elseif nargin == 5     
    badWindowThreshold = 0.4;         
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
elseif nargin == 6      
    robustDeviationThreshold = 5;     
    PowerFrequency = 50;              
    FrequencyNoiseThreshold = 3;      
    flagNotchFilter = 0;              
    correlationThreshold = 0.6;       
elseif nargin == 7
    PowerFrequency = 50;
    FrequencyNoiseThreshold = 3;
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
elseif nargin == 8
    FrequencyNoiseThreshold = 3;
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
elseif nargin == 9
    flagNotchFilter = 0;
    correlationThreshold = 0.6;
else
    correlationThreshold = 0.6;
end
%----------------------------Check Input Validity-------------------------%
if isempty(HighPassband)
    HighPassband = 1;
end
if isempty(selechanns)
    selechanns= 'all';
end
if isempty(badWindowThreshold)
    badWindowThreshold = 0.4;
end
if isempty(robustDeviationThreshold)
    robustDeviationThreshold = 5;
end
if isempty(PowerFrequency)
    PowerFrequency = 50;
end
if isempty(FrequencyNoiseThreshold)
    FrequencyNoiseThreshold = 3;
end
if isempty(flagNotchFilter)
    flagNotchFilter = 0;
end
if isempty(correlationThreshold)
    correlationThreshold = 0.6;
end

% check EEG data
try
    if isempty(data)
        error('EEG.data is empty!!!!!');
    end
catch
    error('EEG.data does not exist!!!!');
end

%----------------------------Initialising Constants-----------------------%
% get sampling rate
try 
    srate = props.samplingInterval; % sampling rate
    if isfinite(srate)
        disp(['sampling rate = ',num2str(srate)]);
    else
        disp('sampling rate is invalide');
        error('sampling rate is invalide');
    end
catch
    disp('sampling rate is not found');
    error('sampling rate is not found');
end
%
% assignin('base',"selechanns",selechanns)
% get selected number of channels and number of channels Nchan
if isequal(selechanns,'all')
    selechanns = 1:props.channelCount;
end
Nchan = length(selechanns); % No. of selected channs
disp(['No. of selected channels: ',num2str(Nchan)]);



 
%--------------------------------Filtering--------------------------------%
% High Pass Filtering
disp('High Pass filtering...');
data = pop_eegfiltnew_raw(data,props,HighPassband,[],[],0);
% notch filtering
disp('Notch filtering for power frequency...');
data = pop_eegfiltnew_raw(data,props,PowerFrequency - 5, PowerFrequency + 5,[],1); % notch filtering for power frequency
%assignin("base","filtereddatausingrawversion",data)
% %------------------------Cut 10s Data into 1s windows---------------------%
% WindowSeconds = 1;
% N_timepoints= size(data,2);
% winlenth = srate * WindowSeconds; % window length (time points)
% Nwin = floor(N_timepoints/winlenth); % No. of windows 
% selectdata = data(selechanns,1:winlenth*Nwin);
% reshapedata = reshape(selectdata,Nchan,winlenth,Nwin);
% disp(['No. of time points: ',num2str(N_timepoints)]);
% disp(['No. of windows: ',num2str(Nwin)]);

data = data(:,1:(end-1)); % excludes the last point as it consistantly displays unrealistic values.
Nwin = 1;

%%
%--------------------------------------------------------------------------------%
% Method 1:
% Detect constant or NaN/Inf signals in each window
disp('------------');
disp('Detecting constant, Inf, or NaN channels......');
median1 = reshape(mad(data, 1, 2),Nchan,1);
std1 = reshape(std(data, 1, 2),Nchan,1);
NanSignalMask = reshape(sum(~isfinite(data),2),Nchan,1);
%assignin('base','NanSignalMaskinworkspace',NanSignalMask)
NoSignalMask = double( median1 < 10e-10 | std1 < 10e-10) + NanSignalMask;
% ---------------------------------------
% Method 2:
% Detect unusually high or low amplitude using robust STD
disp('------------');
disp('Detecting unusually high or low amplitude using robust STD......');
disp(['Robust deviation threshold: ',num2str(robustDeviationThreshold)]);
index1 = abs(data) > 150;  % absolute amplitude > 150 μV
high1 = reshape(sum(index1,2),Nchan,1) > 0;
% assignin("base","high1",high1)
% assignin("base","index1",index1)
channelDeviation = reshape(0.7413 *iqr(data,2),Nchan,1); % Robust estimate of SD
%channelDeviation_cuml = [channelDeviation_cuml channelDeviation];
channelDeviationSD = 0.7413*iqr(channelDeviation(:));
channelDeviationMedian = nanmedian(channelDeviation(:),1);
%assignin("base","channelDeviationMedian",channelDeviationMedian)
robustChannelDeviation = (channelDeviation - channelDeviationMedian) / channelDeviationSD;
HighAmpliMask = abs(robustChannelDeviation) > robustDeviationThreshold | isnan(robustChannelDeviation) | high1;
% ----------------------------------------
% Method 3: Compute the noise-to-signal ratio (based on Christian Kothe's clean_channels)
% Note: RANSAC and global correaltion uses the filtered values X of the data
disp('------------');
if flagNotchFilter == 1
    disp('Detecting high frequency noise and power frequency noise using noise-to-signal ratio......');
else
    disp('Detecting high frequency noise using noise-to-signal ratio......');
end

disp(['Frequency noise threshold: ',num2str(FrequencyNoiseThreshold)]);
% FrequencyNoiseMask = zeros(Nchan,Nwin);
if srate > 2*PowerFrequency
    % Remove signal content above 40Hz/50Hz and below 1 Hz
    disp('low Pass filtering...');
    data_alt = pop_eegfiltnew_raw(data,props,[],PowerFrequency - 10,[],0); % In China, the power frequency is 50Hz, so set the high pass frequency as 50-10=40.
    if flagNotchFilter == 1
        disp('Notch filtering for 0.5*power frequency...');
        data_alt = pop_eegfiltnew_raw(data_alt,props,0.5*PowerFrequency - 5, 0.5*PowerFrequency + 5,[],1); % notch filtering for 0.5* power frequency
    end
    % checking high frequency noise 
    % -------------- 
    winlenth = 1000;
    X = data_alt;
    %X = data_alt(selechanns,1:winlenth*Nwin);
    %X = reshape(X,Nchan,winlenth,Nwin);
    % Determine z-scored level of EM noise-to-signal ratio for each window
    noisiness = mad(data - X, 1, 2)./mad(X, 1, 2);
    noisiness = reshape(noisiness,Nchan,Nwin);
    %noisiness_cuml = [noisiness_cuml,noisiness];
    noisinessMedian = nanmedian(noisiness(:));
    noisinessSD = mad(noisiness(:), 1)*1.4826; % median absolute deviation
    zscoreFreNoiseTemp = (noisiness - noisinessMedian) ./ noisinessSD;
        
%     noisinessMedian = nanmedian(noisiness);
%     noisinessSD = mad(noisiness, 1)*1.4826;
%     zscoreFreNoiseTemp = bsxfun ( @minus, noisiness, noisinessMedian);
%     zscoreFreNoiseTemp = bsxfun ( @rdivide, zscoreFreNoiseTemp,noisinessSD);
    
    FrequencyNoiseMask = (abs(zscoreFreNoiseTemp) > FrequencyNoiseThreshold) | isnan(zscoreFreNoiseTemp) | (abs(noisiness) > 0.5); % or the absolute noise-to-signal ratio > 0.5
    FrequencyNoiseMask = FrequencyNoiseMask .* (abs(noisiness) > 0.0075);   % the error between signals of twice low pass fitering is about 0.0075, so the noisiness < 0.0075 may be indistinguishable.
    %a = (abs(noisiness) > 0.0075) == 1;
    %assignin('base','sevenfive',a);
else
    warning('The sampling rate is below 2*PowerFrequency (too low), detecting high frequency noise is skipped');
    X = EEG.data(selechanns,1:winlenth*Nwin);
    X = reshape(X,Nchan,winlenth,Nwin);
        FrequencyNoiseMask = zeros(Nchan,Nwin);
end
    
% -----------------------------------------
% Method 4: Global correlation criteria in time domain (from Nima Bigdely-Shamlo)
disp('------------');
disp('Detecting low correlation with other channels......');
disp(['correlation threshold: ',num2str(correlationThreshold)]);

channelCorrelations = zeros(Nchan,Nwin);
for k1 = 1:Nwin 
    eegPortion = squeeze(X(:, :, k1))'; % using filtered data X
    %assignin('base','eegPortion',eegPortion)
    windowCorrelation = corrcoef(eegPortion);
    %assignin('base','windowCorrelation',windowCorrelation)
    abs_corr = abs(windowCorrelation - eye(Nchan,Nchan));
    %assignin('base','abs_corr',abs_corr)
    channelCorrelations(:,k1)  = quantile(abs_corr, 0.98); % approximate maximal correlation: quantile of 98%
end

dropOuts = isnan(channelCorrelations) | isnan(noisiness);
channelCorrelations(dropOuts) = 0;
% noisiness(dropOuts) = 0;
LowCorrelationMask = channelCorrelations < correlationThreshold;

OverallBadMask = NoSignalMask + HighAmpliMask + FrequencyNoiseMask + LowCorrelationMask; % considering all methods

for i = 1 : Nchan
    clr = "red"; 
    msg = "";
    res = find([NoSignalMask(i), HighAmpliMask(i), FrequencyNoiseMask(i), LowCorrelationMask(i)] == 1);
    %disp(res==3)
    if length(res) > 1
        msg = "MU";
    elseif isempty(res)
        clr = "green"; 
    else
        switch res
            case 1
                msg = "NS";
                NS_track(1,i)  = NS_track(1,i) + 1;
            case 2
                msg = "HA";
                NS_track(1,i)  = NS_track(1,i) + 1;
            case 3
                msg = "FN";
                NS_track(1,i)  = 0;
            case 4
                msg = "LC";
                NS_track(1,i)  = 0;
        end
    end
    if NS_track(1,i) >= 15
        local_msg = sprintf("%d号电极(%s)可能脱落",i,chann_names{i});
        warndlg(local_msg,"Warning")
        NS_track = zeros(1,Nchan);
    end
    NS_track_New = NS_track;
    set(fig.UserData.lamp_handles(i), 'BackgroundColor', clr);
    set(fig.UserData.msg_handles(i), 'String', msg);
end
% figure;
% subplot(1,4,1)
% plot(NoSignalMask)
% subplot(1,4,2)
% plot(HighAmpliMask)
% subplot(1,4,3)
% plot(FrequencyNoiseMask)
% subplot(1,4,4)
% plot(LowCorrelationMask)
% disp(NoSignalMask)
% disp(HighAmpliMask)
% disp(FrequencyNoiseMask)
% disp(LowCorrelationMask)
% disp(channelDeviation_cuml)
% disp(noisiness_cuml)

end

