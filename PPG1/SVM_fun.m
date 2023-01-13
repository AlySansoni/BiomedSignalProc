%% Function that adds the noise and computes all the STDs and relative labels
% plot_lb is the boolean parameters to the optional print of some plots

function[STDhr,STDamp,STDsd,STDwav,labels] = SVM_fun(SyntData,plot_lb)
%Plotting the original syntetic signal
% SyntData = load('tmp.mat').tmp.v;
samples = 1:105000;
Fs = 500;
t = samples/Fs;
    
% Adding noise to the whole signal
SNR = 25;
VarSig = var(SyntData);
VarNoise = VarSig/power(10,SNR/10);
NoiseData = SyntData + sqrt(VarNoise) * randn(1,105000);

% Adding more noise to a random part of the signal with a fixed length
SNR2 = -25;
VarNoise2 = VarSig/power(10,SNR2/10);
noise_length = 30000;
noise_init = randi(length(SyntData)-noise_length);
noise_end = noise_init + noise_length;
NoiseData2 = NoiseData(noise_init+1:(noise_end)) + sqrt(VarNoise2) * randn(1,noise_length);
NoiseData = [NoiseData(1:(noise_init)) NoiseData2 NoiseData((noise_end+1):end)];

%% Preprocessing the data

% 3th order IIR band pass filter with cut-off frequencies of 0.5 and 12 Hz
fc1 = 0.5;
fc2 = 12;
fs = 500;
[b,a] = butter(3,[fc1/(fs/2) fc2/(fs/2)], 'bandpass');

% Zero-phase forward and reverse filtering
zeroPhase = filtfilt(b,a, NoiseData);

% Plotting pre-processing steps
if plot_lb
    subplot(3,1,1);
    plot(t,SyntData);
    title('Original Synthetic Signal');
    axis([max((noise_init/500)-10,0) min((noise_end/500)+10,210) min(SyntData)-1 max(SyntData)+1]);
    %axis([max((noise_init/500)-10,0) (noise_init/500)+10 min(SyntData)-1 max(SyntData)+1]);
    xlabel('Seconds');
    ylabel('Signal amplitude [au]');
    subplot(3,1,2)
    plot(t,NoiseData);
    title('Signal with noise');
    axis([max((noise_init/500)-10,0) min((noise_end/500)+10,210) min(NoiseData)-1 max(NoiseData)+1]);
    %axis([max((noise_init/500)-10,0) (noise_init/500)+10 min(NoiseData)-1 max(NoiseData)+1]);
    xlabel('Seconds');
    ylabel('Signal amplitude [au]');
    subplot(3,1,3)
    plot(t,zeroPhase);
    title('Filtered Signal');
    axis([max((noise_init/500)-10,0) min((noise_end/500)+10,210) min(zeroPhase)-1 max(zeroPhase)+1]);
    %axis([max((noise_init/500)-10,0) (noise_init/500)+10 min(zeroPhase)-1 max(zeroPhase)+1]);
    xlabel('Seconds');
    ylabel('Signal amplitude [au]');
    %save('Signal','zeroPhase'); %Saving the signal for the 2-nd paper
end
%% Parameters

% Initializing the four parameters matrixes

STDhr = zeros(1,29);
STDamp = zeros(1,29);
STDsd = zeros(1,29);
STDwav = zeros(1,29);

%Computing peaks and troughs of the signal
[peaks,locpk] = findpeaks(zeroPhase);
[troughs, loctr] = findpeaks(-zeroPhase);
troughs = -troughs;

% Splitting the signal in 7-s segments
for n=0:29
    % Computing the right indexes of peaks and troughs inside the segmen
    segment1 = find(locpk >= (1 + n*3500));
    segment2 = find(locpk <= 3500*(n+1));
    segmentPk = intersect(segment2,segment1);
    
    segment1 = find(loctr >= (1 + n*3500));
    segment2 = find(loctr <= 3500*(n+1));
    segmentTr = intersect(segment2,segment1);
    
    % Choosing the right starting point in the segment in order to avoid
    % partial pulses
    if peaks(1,segmentPk(1,2)) > peaks(1,segmentPk(1,1))
        segmentPk = segmentPk(1,2:end); 
    end 
    while locpk(1,segmentPk(1,1)) > loctr(1,segmentTr(1,1))
        segmentTr = segmentTr(1,2:end);
    end
    segmentTr = segmentTr(1,2:end);

    % Standard Deviation of Peak-to-Peak Interval
    pp_interval = zeros(1,length(segmentPk)-1);
    j = 1;
    % Step of two in order to avoid false peaks
    for i=3:2:length(segmentPk)
        pp_interval(j) = (locpk(1,segmentPk(i)) - locpk(1,segmentPk(i-2))) / 500; % Difference in seconds
        j = j+1;
    end
    pp_interval = pp_interval(1,pp_interval~=0);
    mean_ppi = mean(pp_interval(1,:));
    sum_int = 0;
    for k=1:length(pp_interval)
        diff_int = (pp_interval(k)-mean_ppi)^2;
        sum_int = sum_int + diff_int;
    end
    STDhr(n+1) = sqrt(sum_int/length(pp_interval)); % std p-p interval
  
    % Standard Deviation of Peak-to-Peak Amplitude
    % Choosing the same number of peaks and troughs in the segment
    if length(segmentTr) < length(segmentPk)
        segmentPk=segmentPk(1,1:length(segmentTr));
    else
        segmentTr=segmentTr(1,1:length(segmentPk));
    end
    pp_amplitude = zeros(1,length(segmentTr)-1);
    
    j = 1;
    % Step of two in order to avoid false peaks and troughs
    for i=1:2:(length(segmentTr)-2)
        pp_amplitude(j) = (peaks(1,segmentPk(i)) - troughs(1,segmentTr(i+2))) / 500; % Difference in seconds
        j = j+1;
    end
    
    pp_amplitude = pp_amplitude(1,pp_amplitude~=0);
    mean_ppa = mean(pp_amplitude(1,:));
    sum_int = 0;
    for k=1:length(pp_amplitude)
        diff_int = (pp_amplitude(k)-mean_ppa)^2;
        sum_int = sum_int + diff_int;
    end
    STDamp(n+1) = sqrt(sum_int/length(pp_amplitude)); % std p-p amplitude
    
    %Standard Deviation of Systolic and Diastolic Ratio
    j = 1;
    sd_ratio = zeros(1,length(segmentTr));
    % Step of two in order to avoid false peaks and troughs
    for i=3:2:length(segmentTr)
        sd_ratio(j) = (loctr(1,segmentTr(i)) - locpk(1,segmentPk(i))) / (locpk(1,segmentPk(i)) - loctr(1,segmentTr(i-2)));
        sd_ratio(j) = sd_ratio(j)/500; % Difference in seconds
        j = j+1;
    end
    
    sd_ratio = sd_ratio(1,sd_ratio~=0);
    mean_sd_ratio = mean(sd_ratio(1,:));
    sum_int = 0;
    for k=1:length(sd_ratio)
        diff_int = (sd_ratio(k)-mean_sd_ratio)^2;
        sum_int = sum_int + diff_int;
    end
    sd_ratio = sd_ratio(2:end);
    STDsd(n+1) = sqrt(sum_int/length(sd_ratio)); % std p-p amplitude
    
    %Mean-Standard Deviation of Pulse Shape (STDwav)
    seg_length = length(segmentPk);
    step = fix((3500/length(segmentPk))/16); % step by which the samples are taken
    STDwav_m = zeros(1,4); % Choosing 4 central representative peaks of the segment
    for m = -2:2
        % Choosing the reference peak from the centre of the segment
        ref_peak = locpk(segmentPk(fix(seg_length/2))-m);
        % Taking a number of samples around the reference one
        samples = zeroPhase(1,ref_peak-8:step:ref_peak+8);

        mean_samples = mean(samples(1,:));
        sum_int = 0;
        for k=1:length(samples)
            diff_int = (samples(k)-mean_samples)^2;
            sum_int = sum_int + diff_int;
        end      
       STDwav_m(m+3) = sqrt(sum_int/length(samples)); % std p-p amplitude for a pulse
    end
        STDwav(n+1) = sum(STDwav_m)/4; % E[std] which is the mean of the previous result
end
%% Computing labels for training or for the confusion matrix
    
% Declaring clean and corrupted labels
labels = zeros(30,1);
for n = 0:28
    seg_init = 1+n*3500;
    seg_end = 3500*(n+1);
    % Checking if there is intersection with the noisy part
    if (seg_init < noise_init && seg_end < noise_init) || (seg_init > noise_end && seg_end > noise_end)
        continue;
    % Checking if the segment is totally contained in the noise interval
    elseif seg_init >= noise_init && seg_end <= noise_end
        labels(n+1,1) = 1;
    % Checking if the partial intersection is enought to be considered corrupted   
    elseif seg_init < noise_init
        if noise_init - seg_init < (0.85*3500)
            labels(n+1,1) = 1;
        end
    elseif seg_end > noise_end
        if seg_end - noise_end < (0.85*3500)
            labels(n+1,1) = 1;
        end
    end
end
    
%Second way to save labels for the 2-nd
%      if plot_lb
%          save('Labels','labels'); 
%          paper
%      end 
end