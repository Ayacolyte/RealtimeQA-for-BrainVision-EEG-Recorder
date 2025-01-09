clear
close all
clc
% initilisation with loaded simulated data
current_dir= pwd;
lib_path1 = fullfile(current_dir, 'DRT_lib');
lib_path2 = fullfile(current_dir, 'lib');
addpath(lib_path1);
addpath(lib_path2);

load("reveng_data_filtering.mat")
load("1225_RDAtest_workspace.mat")
str = 1001;
ed = 2000;
filter_test_data_before = EEGinworkspace_beforefilter.data;

num_channels = double(props.channelCount);
chann_names = props.channelNames;

HighPassband = 1;
seleChanns = 'all';
badWindowThreshold = 0.6;
robustDeviationThreshold = 5;
PowerFrequency = 50;
FrequencyNoiseThreshold = 3;
flagNotchFilter = 0;
correlationThreshold = 0.6;


%for i = 1 : floor(size(filter_test_data_before,2)/1000)
idx1 = 1;
idx2 = 20;
data1s = [];
NoSignalMask_mat = [];
HighAmpliMask_mat = [];
FrequencyNoiseMask_mat = [];
LowCorrelationMask_mat = [];
OverallBadMask_mat = [];

NS_track = zeros(1,num_channels);


close all
fig1 = createLamps(num_channels,chann_names);
WinOnTop( fig1, true ); % set figure to be always on top, script by Igor
%set(fig1, 'WindowStyle', 'alwaysontop');
%drawnow;
while idx2 <= length(filter_test_data_before)
    tic;
    data = filter_test_data_before(:,idx1:idx2);
    EEGData = data;
    %EEGData = reshape(data, props.channelCount, length(data) / props.channelCount);
    data1s = [data1s EEGData];
    dims = size(data1s);
    elapsedTime = toc;
    remainingTime = 0.020 - elapsedTime;
    if remainingTime > 0
        pause(remainingTime);
    end
    if dims(2) > 1000000 / props.samplingInterval
        data1s = data1s(:, dims(2) - 1000000 / props.samplingInterval : dims(2));
        disp(size(data1s))
        %avg = mean(mean(data1s.*data1s));
        %disp(['Average power: ' num2str(avg)]);
        [~,OverallBadMask,NoSignalMask,HighAmpliMask,FrequencyNoiseMask,LowCorrelationMask,NS_track_new]=realtimeQA(data1s,props,fig1,HighPassband,seleChanns,badWindowThreshold,...
          robustDeviationThreshold,PowerFrequency,FrequencyNoiseThreshold,flagNotchFilter,correlationThreshold,NS_track,chann_names);
        NS_track = NS_track_new;
        % set data buffer to empty for next full second
        data1s = [];
        OverallBadMask(OverallBadMask >= 1) = 1;
        OverallBadMask_mat = [OverallBadMask_mat,OverallBadMask];
        NoSignalMask_mat = [NoSignalMask_mat,NoSignalMask];
        HighAmpliMask_mat = [HighAmpliMask_mat,HighAmpliMask];
        FrequencyNoiseMask_mat = [FrequencyNoiseMask_mat,FrequencyNoiseMask];
        LowCorrelationMask_mat = [LowCorrelationMask_mat,LowCorrelationMask];
    end

    idx1 = idx1 + 20;
    idx2 = idx2 + 20;
end
figure(),clf
subplot(5,3,1)
imagesc(NoSignalMask_mat)
subplot(5,3,4)
imagesc(HighAmpliMask_mat)
subplot(5,3,7)
imagesc(FrequencyNoiseMask_mat)
subplot(5,3,10)
imagesc(LowCorrelationMask_mat)
subplot(5,3,2)
imagesc(OverallBadMask_mat)