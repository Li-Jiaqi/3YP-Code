clear 
close all


load('dryproject_keyboard.mat')
RA = dryproject_Ch7.values;
Oblique = dryproject_Ch8.values;
time_stamp = dryproject_Ch31.times;
load('dryproject.mat')

[r,c] = size(RA);
fs = SmrData.SR;
datapt = 1:1:r;
datapt = datapt.';
time = datapt/fs;


%% filtering, rectification 



window = hamming(size(Oblique,1)); 
[pxx_normal,f_normal]=periodogram(Oblique, window,[],fs,'psd');
figure(8)
plot(f_normal,pxx_normal)
xlim([20,500])

f_low= 300; % attenuate low-frequency noise such as breathing
f_high = 40; 
norder = 4;
d_bp=design(fdesign.bandpass('N,F3dB1,F3dB2',norder,f_high,f_low,fs),'butter');
[b_bp,a_bp] = tf(d_bp);

RA_filtered = filtfilt(b_bp,a_bp,RA);
RA_filtered = abs(RA_filtered);
Oblique_filtered = filtfilt(b_bp,a_bp,Oblique);
Oblique_filtered = abs(Oblique_filtered);

window = hamming(size(Oblique_filtered,1)); 
[pxx_normal,f_normal]=periodogram(Oblique_filtered, window,[],fs,'psd');
plot(f_normal,pxx_normal,'b')
xlim([10,500])
xlabel('Frequency (Hz)')
ylabel('PSD')
title('PSD of raw sEMG recording of EO channel')

figure(1)
subplot(2,1,1)
plot(time,RA)
subplot(2,1,2)
plot(time,RA_filtered)

figure(2)
subplot(2,1,1)
plot(datapt,Oblique)
subplot(2,1,2)
plot(datapt,Oblique_filtered)

%% sliding window average


n_start = 1;
window_size = 1000;
overlap_size = 800;
window_step = window_size - overlap_size;
Oblique_window = zeros(1,floor(r/(window_size -overlap_size)));
RA_window = zeros(1,floor(r/(window_size -overlap_size)));
i = 1;
while (n_start+window_size<=r)
    Oblique_window(i) =sqrt( mean(Oblique_filtered(n_start:n_start+window_size-1,:).^2));
    RA_window(i) =sqrt( mean(RA_filtered(n_start:n_start+window_size-1,:).^2));
    n_start = n_start + window_size - overlap_size;
    i = i+1;
end
Oblique_window=Oblique_window(1:end-2);
RA_window=RA_window(1:end-2);
datapt_window = window_size/2:(window_size -overlap_size):r;
time_window = datapt_window./fs;
figure(3)
subplot(2,1,1)
plot(time,Oblique_filtered,'Color', [150, 150, 150]./255)
hold on
plot(time_window, Oblique_window,'')
xlim([0,220])
legend('before smoothing by sliding window','after smoothing by sliding window')
xlabel('time (s)')
title('Oblique')
subplot(2,1,2)
plot(time,RA_filtered,'Color', [150, 150, 150]./255)
hold on
plot(time_window, RA_window,'')
xlim([0,220])
legend('before smoothing by sliding window','after smoothing by sliding window')
xlabel('time (s)')
title('RA')


%% time stamp segmentation
[standing_cough_time,standing_cough] = extractSeg(Oblique_window,11,time_stamp, time_window );

figure(4)
plot(standing_cough_time,standing_cough)

[squat_time,squat] = extractSeg(Oblique_window,6,time_stamp, time_window );
[legrise_time,legrise] = extractSeg(Oblique_window,5,time_stamp, time_window );
[sit_time,sit] = extractSeg(Oblique_window,1,time_stamp, time_window );
[stand_time,stand] = extractSeg(Oblique_window,2,time_stamp, time_window );




%% extract features - envelop
% low pass filtering for envelop detection 

f_low = 3; 
norder = 4;
d_bp=design(fdesign.lowpass('N,F3db',norder,f_low,50),'butter');
[b_bp,a_bp] = tf(d_bp);

squat_envelop = filtfilt(b_bp,a_bp,squat);
standcough_envelop = filtfilt(b_bp,a_bp,standing_cough);
legrise_envelop = filtfilt(b_bp,a_bp,legrise);
sit_envelop = filtfilt(b_bp,a_bp,sit);
stand_envelop = filtfilt(b_bp,a_bp,stand);

figure(5)
subplot(3,1,1)
plot(squat_time,squat,'Color', [150, 150, 150]./255)
hold on
plot(squat_time,squat_envelop,'b')
title('squat')
ylim([0,0.2])

subplot(3,1,2)
plot(standing_cough_time,standing_cough,'Color', [150, 150, 150]./255)
hold on
plot(standing_cough_time,standcough_envelop,'b')
ylim([0,0.2])
title('standing cough')

subplot(3,1,3)
plot(legrise_time,legrise,'Color', [150, 150, 150]./255)
hold on
plot(legrise_time,legrise_envelop,'b')
ylim([0,0.2])
title('standing leg rise')

figure(6)
subplot(2,1,1)
plot(sit_time,sit,'Color', [150, 150, 150]./255)
hold on
plot(sit_time,sit_envelop,'b')
ylim([0,0.15])
title('sit')
subplot(2,1,2)
plot(stand_time,stand,'Color', [150, 150, 150]./255)
hold on
plot(stand_time,stand_envelop,'b')
ylim([0,0.15])
title('stand')

rest = find(sit_time>16&sit_time<24);
rest_level = mean(sit_envelop(rest));
fprintf('the rest EMG level is %5.3f\n',rest_level)

rest_stand = find(stand_time>25&stand_time<29);
reststand_level = mean(sit_envelop(rest_stand));
fprintf('the stand EMG level is %5.3f\n',reststand_level)


%% First derivative of EMG envelop
standcough_diff = diff(standcough_envelop);
legrise_diff = diff(legrise_envelop);
squat_diff = diff(squat_envelop);

standcough_diff2 = diff(standcough_diff);
legrise_diff2 = diff(legrise_diff);
squat_diff2 = diff(squat_diff);

figure(7)
subplot(3,1,1)
plot(standing_cough_time,standing_cough,'Color', [150, 150, 150]./255)
hold on
plot(standing_cough_time,standcough_envelop,'b')
plot(standing_cough_time(1:end-1),standcough_diff,'r')
plot(standing_cough_time(1:end-2),standcough_diff2,'g')
yline(0.011)
title('standing cough')
ylim([0,0.15])

subplot(3,1,2)
plot(squat_time,squat,'Color', [150, 150, 150]./255)
hold on
plot(squat_time,squat_envelop,'b')
plot(squat_time(1:end-1),squat_diff,'r')
plot(squat_time(1:end-2),squat_diff2,'g')
yline(0.011)
title('squat')
ylim([0,0.15])

subplot(3,1,3)
plot(legrise_time,legrise,'Color', [150, 150, 150]./255)
hold on
plot(legrise_time,legrise_envelop,'b')
plot(legrise_time(1:end-1),legrise_diff,'r')
plot(legrise_time(1:end-2),legrise_diff2,'g')
yline(0.011)
title('standing leg rise')
ylim([0,0.15])



%% amplitude mormalisation & merging 


EO_envelop = filtfilt(b_bp,a_bp,Oblique_window);
RA_envelop = filtfilt(b_bp,a_bp,RA_window);
[MVC_EO,IEO] = max(EO_envelop);
t_max_EO = time_window(IEO);
[MVC_RA,IRA] = max(RA_envelop);
t_max_RA = time_window(IRA);
% amplitude modulation - express as a percentage to MVC

EO_modulated = EO_envelop/MVC_EO*100;
RA_modulated = RA_envelop/MVC_RA*100;
EMG_merge = EO_modulated+RA_modulated;
[MVC,I] = max(EMG_merge);
EMG_modulated = EMG_merge/MVC*100;


figure(9)
subplot(3,1,1)
plot(time_window, EO_modulated,'b')
hold on
plot(t_max_EO,EO_modulated(IEO),'ro')
xlim([0,220])
xlabel('time(s)')
ylim([0,105])
ylabel('percentage of MVC (%)')
title('finding MVC')
subplot(3,1,2)
plot(time_window,RA_modulated,'b')
hold on
plot(t_max_RA,RA_modulated(IRA),'ro')
xlim([0,220])
xlabel('time(s)')
ylim([0,105])
ylabel('percentage of MVC (%)')
title('finding MVC')
subplot(3,1,3)
plot(time_window, EMG_modulated)
hold on 
plot(time_window(I), EMG_modulated(I),'ro')
xlim([0,220])



%% pre-processing demonstration
[standcough_datapt, standcough] = extractSeg(Oblique,1,[185.942,188.227]*fs,datapt);

f_low= 300; % attenuate low-frequency noise such as breathing
f_high = 40; 
norder = 4;
d_bp=design(fdesign.bandpass('N,F3dB1,F3dB2',norder,f_high,f_low,fs),'butter');
[b_bp,a_bp] = tf(d_bp);

standcough_filtered = abs(filtfilt(b_bp,a_bp,standcough));

n_start = 1;
window_size = 800;
overlap_size = 700;
window_step = window_size - overlap_size;
standcough_window = zeros(1,ceil(size(standcough_filtered,1)/(window_size -overlap_size)));

i = 1;
while (n_start+window_size<=size(standcough_filtered,1))
    standcough_window(i) =sqrt( mean(standcough_filtered(n_start:n_start+window_size-1,:).^2));
    n_start = n_start + window_size - overlap_size;
    i = i+1;
end
standcough_window=standcough_window(1:end-4);

datapt_window = window_size/2:(window_size -overlap_size):r;
standcough_time = (185.942*fs+window_size/2:(window_size -overlap_size):188.227*fs)/fs;

f_low = 3; 
norder = 4;
d_bp=design(fdesign.lowpass('N,F3db',norder,f_low,50),'butter');
[b_bp,a_bp] = tf(d_bp);

scsegment_envelop = filtfilt(b_bp,a_bp,standcough_window);
scsegment_modulated = scsegment_envelop/0.1636*100;


figure(18)
subplot(3,1,1)
plot(standcough_datapt/fs,standcough,'k')
xlabel('time (s)')
ylabel('raw sEMG (mV)')
legend('raw sEMG')
title('raw sEMG')
set(gca,'FontSize',20)
subplot(3,1,2)
plot(standcough_datapt/fs,standcough_filtered,'Color',[150 150 150]./255)
hold on 
plot(standcough_time,standcough_window,'b','LineWidth',2)
plot(standcough_time,scsegment_envelop,'r','LineWidth',1.5)
xlabel('time (s)')
ylabel('sEMG (mV)')
ylim([0,0.35])
legend('filtered and rectified sEMG','RMS sliding window smoothing','Envelop extration')
title('filtered and rectified, sliding window smoothed, and envelop of sEMG')
set(gca,'FontSize',20)
subplot(3,1,3)
plot(standcough_time,scsegment_modulated,'r','LineWidth',1.5)
xlabel('time (s)')
ylabel('% of MVC')
ylim([0,100])
legend('Amplitude normalised sEMG')
title('Pre-processed sEMG')
set(gca,'FontSize',20)


%% preprocessing demonstration for ppt
[stand_datapt, stand] = extractSeg(Oblique,1,[26,28]*fs,datapt);
[standcough_datapt, standcough] = extractSeg(Oblique,1,[165,167]*fs,datapt);

%standcough_datapt = cat(1,stand_datapt, standcough_datapt,stand_datapt);
standcough = cat(1,stand(1:end/4*3),standcough,stand, stand);
seg_length = size(standcough,1);
standcough_datapt = 1:1:(seg_length);


f_low= 300; % attenuate low-frequency noise such as breathing
f_high = 40; 
norder = 4;
d_bp=design(fdesign.bandpass('N,F3dB1,F3dB2',norder,f_high,f_low,fs),'butter');
[b_bp,a_bp] = tf(d_bp);

standcough_filtered = abs(filtfilt(b_bp,a_bp,standcough));

n_start = 1;
window_size = 800;
overlap_size = 700;
window_step = window_size - overlap_size;
standcough_window = zeros(1,ceil(size(standcough_filtered,1)/(window_size -overlap_size)));

i = 1;
while (n_start+window_size<=size(standcough_filtered,1))
    standcough_window(i) =sqrt( mean(standcough_filtered(n_start:n_start+window_size-1,:).^2));
    n_start = n_start + window_size - overlap_size;
    i = i+1;
end
standcough_window=standcough_window(1:end-4);

datapt_window = window_size/2:(window_size -overlap_size):r;
standcough_time = (0*fs+window_size/2:(window_size -overlap_size):7.5*fs)/fs;

f_low = 3; 
norder = 4;
d_bp=design(fdesign.lowpass('N,F3db',norder,f_low,50),'butter');
[b_bp,a_bp] = tf(d_bp);

scsegment_envelop = filtfilt(b_bp,a_bp,standcough_window);
scsegment_modulated = scsegment_envelop/0.1636*100;


figure(18)
subplot(3,1,1)
rectangle('Position',[3.2,-0.5,0.9,5],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
hold on
plot(standcough_datapt/fs,standcough,'k')
xlabel('time (s)')
%ylabel('raw sEMG (mV)')
%legend('raw sEMG')
title('Raw sEMG')
text([1.5,3.4,5.8],[0.4,0.4,0.4],{'Stand','Valsalva','Stand'},'FontSize',20,'Color','k')
set(gca,'FontSize',20)
xlim([0,7.5])
ylim([-0.5,0.5])

subplot(3,1,2)
rectangle('Position',[3.2,-0.5,0.9,5],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
hold on
plot(standcough_datapt/fs,standcough_filtered,'Color',[150 150 150]./255)
plot(standcough_time,standcough_window,'b','LineWidth',2)
plot(standcough_time,scsegment_envelop,'r','LineWidth',1.5)

xlabel('time (s)')
%ylabel('sEMG (mV)')
ylim([0,0.2])
%legend('filtered and rectified sEMG','RMS sliding window smoothing','Envelop extration')
title('After Filtering')
set(gca,'FontSize',20)
xlim([0,7.5])

subplot(3,1,3)
rectangle('Position',[3.2,0,0.9,70],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
hold on
plot(standcough_time,scsegment_modulated,'r','LineWidth',1.5)
xlabel('time (s)')
ylabel('% of MVC')
ylim([0,60])
%legend('Amplitude normalised sEMG')
title('After Amplitude Modulation')
set(gca,'FontSize',20)
xlim([0,7.5])
%% following



%% combining EMG channels demonstration
[legriseseg_time,legriseseg_EO] = extractSeg(EO_modulated,1,[88.9248,92.496], time_window );
[legriseseg_time,legriseseg_RA] = extractSeg(RA_modulated,1,[88.9248,92.496], time_window );
[legriseseg_time,legriseseg] = extractSeg(EMG_modulated,1,[88.9248,92.496], time_window );
[jumpseg_time,jumpseg_EO] = extractSeg(EO_modulated,1,[201.091,202.55], time_window );
[jumpseg_time,jumpseg_RA] = extractSeg(RA_modulated,1,[201.091,202.55], time_window );
[jumpseg_time,jumpseg] = extractSeg(EMG_modulated,1,[201.091,202.55], time_window );
figure(19)
subplot(2,1,1)
plot(legriseseg_time,legriseseg_EO,'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
hold on 
plot(legriseseg_time,legriseseg_RA,'Color',[0.3010 0.7450 0.9330],'LineWidth',1)
plot(legriseseg_time,legriseseg,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
ylim([0,60])
xlim([88.75,92.75])
xlabel('time (s)')
ylabel('sEMG (% of MVC)')
legend('External Oblique','Rectus Abdominis','Merged Signal')
title('One Right Leg Raise')
subplot(2,1,2)
plot(jumpseg_time,jumpseg_EO,'Color',[0.9290 0.6940 0.1250],'LineWidth',1)
hold on
plot(jumpseg_time,jumpseg_RA,'Color',[0.3010 0.7450 0.9330],'LineWidth',1)
plot(jumpseg_time,jumpseg,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
ylim([0,60])
xlabel('time (s)')
ylabel('sEMG (% of MVC)')
legend('External Oblique','Rectus Abdominis','Merged Signal')
title('One Jump')


%% Relating EMG amplitude to IAP
[sit_time,sit] = extractSeg(EMG_envelop,1,time_stamp, time_window );
[stand_time,stand] = extractSeg(EMG_envelop,2,time_stamp, time_window );
[bendknee_time,bendkee] = extractSeg(EMG_envelop,3,time_stamp, time_window );
[bendwaist_time,bendwaist] = extractSeg(Oblique_window,4,time_stamp, time_window );
[legrise_time,legrise] = extractSeg(EMG_envelop,5,time_stamp, time_window );
[squat_time,squat] = extractSeg(EMG_envelop,6,time_stamp, time_window );
[sitstand_time,sitstand] = extractSeg(EMG_envelop,7,time_stamp, time_window );
[sitvalsalva_time,sitvalsalva] = extractSeg(Oblique_window,8,time_stamp, time_window );
[standvalsalva_time,standvalsalva] = extractSeg(Oblique_window,9,time_stamp, time_window );
[sitcought_time,sitcough] = extractSeg(Oblique_window,8,time_stamp, time_window );
[standing_cough_time,standing_cough] = extractSeg(EMG_envelop,11,time_stamp, time_window );
[jump_time,jump] = extractSeg(Oblique_window,12,time_stamp, time_window );

[peak,I] = findpeaks(EMG_envelop, 'MinPeakHeight',0.02,'MinPeakDistance',1);




%% first and second derivative
EMG_diff = diff(EMG_modulated);
EMG_diff2 = diff(EMG_diff);

[peak1,I1] = findpeaks(EMG_diff, 'MinPeakHeight',3,'MinPeakDistance',1);
[peak2,I2] = findpeaks(EMG_diff2, 'MinPeakHeight',1,'MinPeakDistance',1);

figure(10)
subplot(3,1,1)
plot(time_window, EMG_modulated,'Color', [150, 150, 150]./255)
hold on
plot(time_window(I), peak,'ro')
xline(time_stamp)

subplot(3,1,2)
plot(time_window(1:end-1), EMG_diff,'Color', [150, 150, 150]./255)
hold on
plot(time_window(I1), peak1,'ro')
xline(time_stamp)

subplot(3,1,3)
plot(time_window(1:end-2), EMG_diff2,'Color', [150, 150, 150]./255)
hold on
plot(time_window(I2), peak2,'ro')
xline(time_stamp)


figure(11)
subplot(3,1,1)
plot(time_window, EMG_envelop)
hold on
xline(time_stamp)
xlim([0,215])

subplot(3,1,2)
plot(time_window(1:end-1), EMG_diff)
hold on
xline(time_stamp)
xlim([0,215])

subplot(3,1,3)
plot(time_window(1:end-2), EMG_diff2)
hold on
xline(time_stamp)
xlim([0,215])

%% building the test segment

seg1 = find(time_window>16&time_window<24);  % sit
seg2 = find(time_window>147.6 & time_window<150);   % sit to stand
seg12 = find(time_window>26 & time_window<28);   % stand
seg11 = find(time_window>28 & time_window<29); 
seg3 = cat(2,seg12,seg12);
seg4 = find(time_window>126 & time_window<131.7);      % squat
seg5 = find(time_window>111.6 & time_window<114);      % standing leg rise
seg6 = find(time_window>166 & time_window<167);        % sitting valsalva
seg7 = find(time_window>186.6 & time_window<188.3);    % stand cough
seg8 = find(time_window>205 & time_window<206.3);      % jump

testI = cat(2,seg1, seg2,seg3,seg3,seg4,seg3,seg5,seg3,seg6,seg3,seg7,seg7,seg3,seg8,seg3,seg3);

length = size(testI, 2);
test_time = 0:30e-3:(length-1)*30e-3;

figure(12)
plot(test_time,EMG_modulated(testI))

%% method1: single threshold - abitrarily determined
EMG_amplitude = xlsread('EMG recording.xlsx','B2:F13');
IAP = xlsread('EMG recording.xlsx','G2:K13');

EMG_amplitude_modulated = xlsread('EMG recording.xlsx','E20:I31');


[xData, yData] = prepareCurveData( EMG_amplitude_modulated, IAP );

% Set up fittype and options.
ft = fittype( '280/(1+exp(-(x+b)/c))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [100,80];

% Fit model to data.
[EMG_gaussian, gof] = fit( xData, yData, ft, opts );

% Fit to linear model
ft = fittype( 'poly1' );

% Fit model to data.
[EMG_linear, gof] = fit( xData, yData, ft );

% Set up fittype and options.
ft = fittype( 'a+b.*log(x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [280 10];

% Fit model to data.
[EMG_log, gof] = fit( xData, yData, ft, opts );

%Goodness of fit:
  %SSE: 1.852e+05
  %R-square: 0.5338
  %Adjusted R-square: 0.5256
  %RMSE: 57

% Goodness of fit: Linear
  %SSE: 1.834e+05
  %R-square: 0.5459
  %Adjusted R-square: 0.5381
  %RMSE: 56.24  
 
  
x = linspace(0,100,1000);
y_linear = 2.576*x+4.883;
y_log = -206.6+ 89.86 .*log(x);


% Plot fit with data.
figure(13)
subplot(2,1,1)
scatter(EMG_amplitude_modulated(1:2,:), IAP(1:2),250,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',0.7);
hold on
scatter(EMG_amplitude_modulated(3:5,:), IAP(3:5),250,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',0.7);
scatter(EMG_amplitude_modulated(6:8,:), IAP(6:8),250,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.7);
scatter(EMG_amplitude_modulated(9:12,:), IAP(9:12),250,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerFaceAlpha',0.7);

p1 = plot (EMG_gaussian,'k','predfunc',0.95);

%legend([p1,p2],{'Gaussian model','Linear Model'})

set(p1, 'LineWidth',1.5)
xlabel('percenrage of MVC (%)')
ylabel('Estimated IAP rise (cmH_2O)')
grid on
title('Best fit sigmoidal curve')
ylim([0,300])

%patch([0,100,100,0],[0,0,40,40],[0.4660 0.6740 0.1880],'EdgeColor','none','FaceAlpha',0.3)
%patch([0,100,100,0],[40,40,65,65],[0.9290 0.6940 0.1250],'EdgeColor','none','FaceAlpha',0.3)
%patch([0,100,100,0],[65,65,120,120],[0.8500 0.3250 0.0980],'EdgeColor','none','FaceAlpha',0.3)
%patch([0,100,100,0],[120,120,300,300],[0.6350 0.0780 0.1840],'EdgeColor','none','FaceAlpha',0.3)
l1=yline(280,'r--','LineWidth',1);

legend(p1,{'Best-fit sigmoidal model','95% confidence interval'})
text(2,260,'Estimated Maximum IAP during Exercise','Color','r','FontSize',20)
set(gca,'FontSize',20)

subplot(2,1,2)
scatter(EMG_amplitude_modulated(1:2,:), IAP(1:2),250,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',0.7);
hold on
scatter(EMG_amplitude_modulated(3:5,:), IAP(3:5),250,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',0.7);
scatter(EMG_amplitude_modulated(6:8,:), IAP(6:8),250,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.7);
scatter(EMG_amplitude_modulated(9:12,:), IAP(9:12),250,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerFaceAlpha',0.7);
scatter(EMG_amplitude_modulated(1:2,:), IAP(1:2),250,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',0.7);

p2 = plot(x,y_linear,'k-.');
p3 = plot(x,y_log,'k-');


legend([p2,p3],{'best-fit Linear Model','best-fit Logarithmic Model'})

xlabel('percenrage of MVC (%)')
ylabel('Estimated IAP rise (cmH_2O)')
grid on
title('Best fit linear and logarithmic curves')
ylim([0,300])
set(gca,'FontSize',20)


%% fit the sigmoid model to the first and second derivative 
% Set up fittype and options.
ft = fittype( '280/(1+exp(-(x+b)/c))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [100,80];

diff1_amplitude = xlsread('EMG recording.xlsx','M2:Q13');
IAP = xlsread('EMG recording.xlsx','G2:K13');

diff2_amplitude= xlsread('EMG recording.xlsx','S2:W13');

[diff1Data, yData] = prepareCurveData( diff1_amplitude, IAP );
[diff1_sigmoid, gof] = fit( diff1Data, yData, ft, opts );

[diff2Data, yData] = prepareCurveData( diff2_amplitude, IAP );
[diff2_sigmoid, gof] = fit( diff2Data, yData, ft, opts );




% Plot fit with data.
figure(14)

subplot(2,1,1)
scatter(diff1_amplitude(1:2,:), IAP(1:2),250,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',0.5);
hold on
scatter(diff1_amplitude(3:5,:), IAP(3:5),250,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',0.5);
scatter(diff1_amplitude(6:8,:), IAP(6:8),250,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.5);
scatter(diff1_amplitude(9:12,:), IAP(9:12),250,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerFaceAlpha',0.5);

p1 = plot (diff1_sigmoid,'k','predfunc',0.95);


set(p1, 'LineWidth',1.5)
xlabel('1^{st} derivate of EMG amplitude')
ylabel('Estimated IAP rise (cmH_2O)')
grid on
title(' Estimated IAP against the First Derivative of sEMG')
legend(p1,{'fitted Sigmoidal model','95% confidence interval'})
set(gca,'FontSize',20)

subplot(2,1,2)
scatter(diff2_amplitude(1:2,:), IAP(1:2),250,'filled','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerFaceAlpha',0.5);
hold on
scatter(diff2_amplitude(3:5,:), IAP(3:5),250,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerFaceAlpha',0.5);
scatter(diff2_amplitude(6:8,:), IAP(6:8),250,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerFaceAlpha',0.5);
scatter(diff2_amplitude(9:12,:), IAP(9:12),250,'filled','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerFaceAlpha',0.5);

p1 = plot (diff2_sigmoid,'k','predfunc',0.95);
set(p1, 'LineWidth',1.5)
xlabel('2^{nd} derivate of sEMG amplitude')
ylabel('Estimated IAP rise (cmH_2O)')
grid on
set(gca,'FontSize',20)

title('Estimated IAP against the Second Derivative of EMG')
legend(p1,{'fitted Sigmoidal model','95% confidence interval'})

%Goodness of fit:
 % SSE: 1.091e+05
 % R-square: 0.73
 % Adjusted R-square: 0.7254
 % RMSE: 43.36

%Goodness of fit: 2nd derivative
  %SSE: 9.531e+04
  %R-square: 0.7641
  %Adjusted R-square: 0.76
  %RMSE: 40.54
%% thresholding

% single threshold on EMG
% sliding window of 3 samples

% calculating mu and sigma using resting EMG
seg9 =find(time_window>132 & time_window<138.5);
seg10 = find(time_window>196 & time_window<200);


rest_sampleI = cat(2,seg1,seg3, seg9,seg10);
mu = mean(EMG_modulated(rest_sampleI));
sigma = std(EMG_modulated(rest_sampleI));

% window length
W = 10;       % time increment between each sample is 20 ms, hence 3 samples cover the range of 60 ms
length = size(testI,2);
i =1;      % pointer
test_sample = EMG_modulated(testI);
threshold = 10.2;
FES = zeros([1,length-(W-1)]);
H = zeros([1,length-(W-1)]);
test_bar = zeros([1,length-(W-1)]);

% apply the threshold
while((i+W-1)<length)
    test_bar(i) = sum(test_sample(i:(i+W-1)))/W;

    if test_bar(i)>threshold
        FES(i)=1;
    end
    i=i+1;

end

% actuation target
target = zeros([1,2915]);
target(find(test_time>13&test_time<18))=0.5;
target(find(test_time>30&test_time<37.35))=0.5;
target(find(test_time>44.2&test_time<49.14))=0.5;
target(find(test_time>54.8&test_time<55.47))=0.5;
target(find(test_time>62.13&test_time<66.12))=0.5;
target(find(test_time>73&test_time<75.45))=0.5;

x = linspace(0,87.42,2915);

%% creating IAP simulation 
peak = 223;
rest = 32.5;

tao1 = 0.1;
tao2 = 1;

x1= linspace(0,2.525,100);
y1 = exp((x1-2)/tao1)+rest;
x2 = linspace(2.54,4.3,100);
y2 = peak* ones(1,size(x2,2));
x3 = linspace(4.3, 10, 500);
y3 = (peak-rest)*(exp(-1*(x3-4.3)/tao2))+rest;
x3 = (x3-4.3)/5.7;

xIAP_sample = cat(2,x1,x2);
xIAP_sample = xIAP_sample(1:124)/2.98;
yIAP_sample = cat(2,y1,y2);
yIAP_sample = yIAP_sample(77:200);
time_sample = linspace(0,10,100);
IAP_sit = ones(1,100)*32;
IAP_stand = ones(1,100)*33;


inc1= (yIAP_sample-33)/190*48.9+33;
y31 = (y3-33)/190*48.9+33;
inc2 = (yIAP_sample-33)/190*42+33;
y32 = (y3-33)/190*42+33;
inc3 = (yIAP_sample-33)/190*17.5+33;
y33 = (y3-33)/190*17.5+33;
inc4 = (yIAP_sample-33)/190*154.6+33;
y34 = (y3-33)/190*154.6+33;
inc5 = (yIAP_sample-33)/190*190+33;
y35 = (y3-33)/190*190+33;
inc6 =  (yIAP_sample-33)/190*242.4+33;
y36 = (y3-33)/190*242.4+33;
%simux = cat(2,time_sample/10*13,xIAP_sample(1:200)/10*6+13,time_sample/10*11+19,xIAP_sample/10*7.35+30,time_sample/10*6.85+37.35,xIAP_sample/10*4.94+44.2,time_sample/10*4.66+49.14,xIAP_sample/10*7+53.8,time_sample/10*0.33+60.8,xIAP_sample+61.13,time_sample/10*0+71.13,xIAP_sample/10*14+71.13,time_sample/10*10.43+85.13);
%simuy = cat(2,IAP_sit,inc1,IAP_stand,inc2,IAP_stand,inc3,IAP_stand,inc4,IAP_stand,inc5,IAP_stand,inc6,IAP_stand);

simux = cat(2,time_sample/10*13,xIAP_sample*5+13,x3*1+18,time_sample/10*11+19,xIAP_sample*7.35+30,x3+37.35,time_sample/10*5.85+38.35,xIAP_sample*4.94+44.2,x3*0.5+49.14,time_sample/10*5.16+49.64,xIAP_sample*0.67+54.8,x3*5+55.47,time_sample/10*1.66+60.47,xIAP_sample*3.99+62.13,x3*5.7+66.12,time_sample/10*1.18+71.82,xIAP_sample*2.45+73,x3*6+75.45,time_sample/10*10.43+81.45);
simuy = cat(2,IAP_sit,inc1,y31,IAP_stand,inc2,y32,IAP_stand,inc3,y33,IAP_stand,inc4,y34,IAP_stand,inc5,y35,IAP_stand,inc6,y36,IAP_stand);
yline_coor = [13,18,30,37.35,44.2,49.14,54.8,55.47,62.13,66.12,73,75.45];


actuation_target = zeros(1,size(simuy,2));
for i =1:size(simuy,2)
    if simuy(i)>40
        actuation_target(i)=1;
    end
    
end



figure(15)
subplot(3,1,1)
plot(test_time(1:end-9),test_bar,'Color',[0.4940 0.1840 0.5560],'LineWidth',1)
hold on 
yline(10.2,'r--','LineWidth',1)
area(x,FES*threshold,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.2)
text([5,13.5,20,31,45,53.5,62.5,73.5],ones(1,8)*75,{'sit','sit-to-stand','stand','squat','leg raise','Valsalva','two coughs','jump'},'FontSize',13,'Color',[150 150 150]./255)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,18,'amplitude threshold = 10.2','Color','r','FontSize',13)
ylim([0,100])
ylabel('Modulated EMG (% of MVC)')
xlabel('time (s)')
legend('sEMG amplitude','amplitude threshold','threshold triggers')
title('sEMG and threshold triggering')



subplot(3,1,2)

plot(simux,actuation_target*0.05,'c','LineWidth',1.5)
hold on
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
area(x,FES*0.1,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
ylim([0,0.5])
legend('Detection target (posture change/sound)','threshold triggers')
yticks([])
xlabel('time(s)')
title('Comparing target and treshold triggering episodes')


subplot(3,1,3)
plot(simux,simuy,'k')
hold on 

xlabel('time (s)')
ylabel('simulated IAP (cmH_2O)')
yline(70,'--','LineWidth',1,'Color',[0.8500 0.3250 0.0980])
yline(70-30,'--','LineWidth',1,'Color',[0.3010 0.7450 0.9330])
xline (yline_coor,'-.','Color',[150 150 150]./255)
area(simux,actuation_target*40,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
text(1,78,'input IAP threshold = 70','Color',[0.8500 0.3250 0.0980],'FontSize',13)
text(1,48,'actual IAP threshold = 33','Color',[0.3010 0.7450 0.9330],'FontSize',13)
xlim([0,90])
legend('Simulated IAP change','input IAP threshold','actual threshold (bias towards false positive)')
title('Simulated IAP and threshold values')





%% thresholding first derivative
test_sample_diff1 = abs(diff(test_sample));



rest_sampleI = cat(2,seg1,seg3, seg9,seg10);
mu = mean(diff(EMG_modulated(rest_sampleI)));
sigma = std(diff(EMG_modulated(rest_sampleI)));

% window length
W = 10;       % time increment between each sample is 20 ms, hence 3 samples cover the range of 60 ms
length = size(testI,2);
i =1;      % pointer

threshold1 = 0.37;
FES1 = zeros([1,length-(W-1)]);
H = zeros([1,length-(W-1)]);
test_bar1 = zeros([1,length-(W-1)]);

% apply the threshold
while((i+W-1)<length)
    test_bar1(i) = sum(test_sample_diff1(i:(i+W-1)))/W;
    
    if test_bar1(i)>threshold1
        FES1(i)=1;
    end
    i=i+1;

end


x1 = linspace(0,87.42,size(test_bar1,2));

figure(16)
subplot(2,1,1)
plot(test_time(1:end-1),test_sample_diff1,'Color',[150 150 150]./255)
hold on 
plot(test_time(1:end-(W-1)),test_bar1,'Color',[0.4660 0.6740 0.1880])
plot(test_time(1:end-(W-1)),H,'k')
xlabel('time (s)')
ylabel('Modulated EMG (% of MVC)')
legend('first derivative of pre-processed EMG','Further smoothed with a window of length 10','normalized EMG with respect to background noise')

subplot(2,1,2)
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
hold on
ylim([0,0.5])
area(x1,FES1*0.1,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.5)
legend('Detection target','Detected IAP spike by single threshold method')
yticks([])
xlabel('time')


%% single threshold algorithm demonstration 
test_sample_diff2 = abs(diff(test_sample_diff1));

rest_sampleI = cat(2,seg1,seg3, seg9,seg10);
mu = mean(diff(EMG_modulated(rest_sampleI)));
sigma = std(diff(EMG_modulated(rest_sampleI)));

% window length
W = 10;       % time increment between each sample is 20 ms, hence 3 samples cover the range of 60 ms
length = size(test_sample_diff2,2);
i =1;      % pointer

threshold2 = 0.1;
FES2 = zeros([1,length-(W-1)]);
H = zeros([1,length-(W-1)]);
test_bar2 = zeros([1,length-(W-1)]);

% apply the threshold
while((i+W-1)<length)
    test_bar2(i) = sum(test_sample_diff2(i:(i+W-1)))/W;
    
    if test_bar2(i)>threshold2
        FES2(i)=1;
    end
    i=i+1;

end

x2 = linspace(0,87.42,size(test_bar2,2));

figure(20)
subplot(3,2,1)
plot(test_time(1:end-9),test_bar,'Color',[0.4940 0.1840 0.5560],'LineWidth',1)
hold on 
yline(10.2,'r--','LineWidth',1)
area(x,FES*threshold,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.2)
text([5,13.5,23,31,45,53.5,62.5,73.5],[60,45,60,45,60,45,60,45],{'sit','sit-to-stand','stand','squat','leg raise','Valsalva','two coughs','jump'},'FontSize',13,'Color',[150 150 150]./255)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,18,'amplitude threshold = 10.2','Color','r','FontSize',13)
ylim([0,100])
xlim([52,57])
ylabel('Modulated EMG (% of MVC)')
xlabel('time (s)')
legend('sEMG absolute amplitude','amplitude threshold','threshold triggers')
title('sEMG absolute amplitude')
set(gca,'FontSize',15)


subplot(3,2,2)
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
hold on
area(x,FES*0.1,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
ylim([0,0.4])
xlim([52,57])
legend('Detection target (posture change/sound)','threshold triggers')
yticks([])
xlabel('time(s)')
title('Thresholding outcome using absolute amplitude')
set(gca,'FontSize',15)

subplot(3,2,3)
plot(test_time(1:end-9),test_bar1,'Color',[0 0.4470 0.7410],'LineWidth',1)
hold on 
yline(threshold1,'r--','LineWidth',1)
area(x1,FES1*threshold1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.2)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,1,'amplitude threshold = 0.37','Color','r','FontSize',13)
ylim([0,4])
xlim([52,57])
ylabel('Modulated EMG (% of MVC)')
xlabel('time (s)')
legend('first differeitation','amplitude threshold','threshold triggers')
title('sEMG first differentiation')
set(gca,'FontSize',15)

subplot(3,2,4)
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
hold on
area(x1,FES1*0.1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
ylim([0,0.4])
xlim([52,57])
legend('Detection target (posture change/sound)','threshold triggers')
yticks([])
xlabel('time(s)')
title('Thresholding outcome using the first differentiation')
set(gca,'FontSize',15)

subplot(3,2,5)
plot(test_time(1:end-11),test_bar2,'Color',[0.4660 0.6740 0.1880],'LineWidth',1)
hold on 
yline(threshold2,'r--','LineWidth',1)
area(x2,FES2*threshold2,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.2)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,0.5,'amplitude threshold = 0.08','Color','r','FontSize',13)
ylim([0,2])
xlim([52,57])
ylabel('Modulated EMG (% of MVC)')
xlabel('time (s)')
legend('second differentiation','amplitude threshold','threshold triggers')
title('sEMG second differentiation')
set(gca,'FontSize',15)

subplot(3,2,6)
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
hold on
area(x2,FES2*0.1,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
ylim([0,0.4])
xlim([52,57])
legend('Detection target (posture change/sound)','threshold triggers')
yticks([])
xlabel('time(s)')
title('Thresholding outcome using the second differentiation')
set(gca,'FontSize',15)

%% modified graph
on = zeros(size(FES));
on(FES==0)= NaN;

on1= zeros(size(FES1));
on1(FES1==0)= NaN;


on2= zeros(size(FES2));
on2(FES2==0)= NaN;


figure(20)
subplot(3,1,1)
rectangle('Position',[54.8,0,0.67,100],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[53,0,1.8,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,2.65,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
plot(x,on,'r','LineWidth',6)
plot(test_time(1:end-9),test_bar,'Color','k','LineWidth',1.5)
hold on 
yline(10.2,'m','LineWidth',1.5)
%stairs(x,FES*threshold,'LineWidth',2,'Color','r')
text([5,13.5,23,31,45,53.8,55,56,62.5,73.5],[60,45,60,45,60,60,60,60,60,45],{'sit','sit-to-stand','stand','squat','leg raise','Rest','Valsalva','Rest','two coughs','jump'},'FontSize',13,'Color',[150 150 150]./255,'FontSize',20)
xline (yline_coor,'-.','Color',[150 150 150]./255,'LineWidth',0.8)
text(53.1,18,'amplitude threshold = 10.2','Color','m','FontSize',20)
ylim([0,100])
xlim([53,57])
ylabel('Normalised EMG (% of MVC)')
xlabel('time (s)')
legend('sEMG absolute amplitude','Detection algo triggered','magnitude threshold')
title('sEMG absolute amplitude')
set(gca,'FontSize',18)



subplot(3,1,2)
rectangle('Position',[54.8,0,0.67,100],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[53,0,1.8,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,2.65,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
plot(test_time(1:end-9),test_bar1,'Color','k','LineWidth',1.5)
hold on 
plot(x1,on1,'r','LineWidth',6)
yline(threshold1,'m','LineWidth',1.5)
%area(x1,FES1*threshold1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.2)
xline (yline_coor,'-.','Color',[150 150 150]./255,'LineWidth',0.8)
text(53.1,1,'amplitude threshold = 0.37','Color','m','FontSize',20)
ylim([0,4])
xlim([53,57])
ylabel('EMG 1^{st} differentiation')
xlabel('time (s)')
legend('first differeitation','Detection algo triggered','magnitude threshold')
title('sEMG first differentiation')
set(gca,'FontSize',18)



subplot(3,1,3)
rectangle('Position',[54.8,0,0.67,100],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[53,0,1.8,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,2.65,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
plot(test_time(1:end-11),test_bar2,'Color','k','LineWidth',1.5)
hold on
plot(x2-0.03,on2,'r','LineWidth',6)
yline(threshold2,'m','LineWidth',1.5)
%area(x2,FES2*threshold2,'FaceColor',[0.4660 0.6740 0.1880],'FaceAlpha',0.2)
xline (yline_coor,'-.','Color',[150 150 150]./255,'LineWidth',0.8)
text(53.1,0.5,'amplitude threshold = 0.08','Color','m','FontSize',20)
ylim([0,1.5])
xlim([53,57])
ylabel('EMG 2^{nd} differentiation')
xlabel('time (s)')
legend('second differentiation','Detection algo triggered','magnitude threshold')
title('sEMG second differentiation')
set(gca,'FontSize',18)




%% latency demonstration

latency = xlsread('EMG recording.xlsx','M45:AP46');

xvalue = [1 1 1 2 2 2 3 3 3 3];
figure(21)

%scatter(latency(11:20,1), latency(11:20,2),250,'MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
%scatter(latency(21:30,1),latency(21:30,2),250,'MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
subplot(1,3,1)
plot(latency(1,1:3),latency(2,1:3),'.-','Color',[0.9290 0.6940 0.1250],'LineWidth',1,'MarkerSize',20)
hold on
plot(latency(1,4:6),latency(2,4:6),'.-','Color',[0.9290 0.6940 0.1250],'LineWidth',1,'MarkerSize',20)
plot(latency(1,7:9),latency(2,7:9),'.-','Color',[0.9290 0.6940 0.1250],'LineWidth',1,'MarkerSize',20)
xlim([0.5,3.5])
xticks([1,2,3])
xticklabels({'absolute amplitude','first differentiation','second differentiation'})
ylabel('Detection latency (s)')

set(gca,'FontSize',20)
title('Low Intensity Activities','FontSize',20)

subplot(1,3,2)
plot(latency(1,10:12),latency(2,10:12),'.-','Color',[0.8500 0.3250 0.0980],'LineWidth',1,'MarkerSize',20)
hold on
plot(latency(1,13:15),latency(2,13:15),'.-','Color',[0.8500 0.3250 0.0980],'LineWidth',1,'MarkerSize',20)
plot(latency(1,16:18),latency(2,16:18),'.-','Color',[0.8500 0.3250 0.0980],'LineWidth',1,'MarkerSize',20)
xlim([0.5,3.5])
xticks([1,2,3])
xticklabels({'absolute amplitude','first differentiation','second differentiation'})
ylabel('Detection latency (s)')
ylim([-0.3,0.3])
set(gca,'FontSize',20)
title('High Intensity Activities','FontSize',20)

subplot(1,3,3)
plot(latency(1,19:21),latency(2,19:21),'.-','Color',[0.6350 0.0780 0.1840],'LineWidth',1,'MarkerSize',20)
hold on
plot(latency(1,22:24),latency(2,22:24),'.-','Color',[0.6350 0.0780 0.1840],'LineWidth',1,'MarkerSize',20)
plot(latency(1,25:27),latency(2,25:27),'.-','Color',[0.6350 0.0780 0.1840],'LineWidth',1,'MarkerSize',20)
plot(latency(1,28:30),latency(2,28:30),'.-','Color',[0.6350 0.0780 0.1840],'LineWidth',1,'MarkerSize',20)
ylim([-0.3,0.3])
xlim([0.5,3.5])
xticks([1,2,3])
xticklabels({'absolute amplitude','first differentiation','second differentiation'})
ylabel('Detection latency (s)')
set(gca,'FontSize',20)
title('Extremely High Intensity Acitivities','FontSize',20)

%% double threshold
% if more than alpha consecutive samples lower than threshold -->
% results in a state change ON to OFF
m =55;      % consecutive 0.4s
j = m; % pointer
FES1double = zeros([1,size(FES1,2)]);
while(j<size(FES1,2))
    if sum(FES1((j-m+1):j))>0
        FES1double(j)=1;
    else
        FES1double(j)=0;
    end
    
    j = j+1;
    
end


% improve the plot presentation
% find ON/OFF time
switch_time=[];
switch_time_p = [];
for i = 1:size(FES,2)-1
    
    if FES(i+1)-FES(i)~=0
        switch_time_p(end+1)=(test_time(i));
    end
    
end

for i = 1:size(switch_time_p,2)/2
    switch_time(1,i) = switch_time_p(i*2-1);
    switch_time(2,i) = switch_time_p(i*2);
    
end

y2 = ones(2,size(switch_time_p,2)/2);
y2 = y2.*0.01;

switch_time_target=[];
switch_time_p_target = [];
for i = 1:size(target,2)-1
    
    if target(i+1)-target(i)~=0
        switch_time_p_target(end+1)=(test_time(i));
    end
    
end

%for i = 1:size(switch_time_p_diff1,2)/2
    %switch_time_target(1,i) = switch_time_p_target(i*2-1);
    %switch_time_target(2,i) = switch_time_p_target(i*2);
    
%end

y3 = ones(2,size(switch_time_p_target,2)/2);
y3 = y3.*0.05;

figure(17)
subplot(2,1,1)
plot(test_time(1:end-9),test_bar1,'Color',[0 0.4470 0.7410],'LineWidth',1)
hold on 
yline(threshold1,'r--','LineWidth',1)
area(x1,FES1*threshold1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.2)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,1,'amplitude threshold = 0.37','Color','r','FontSize',13)
ylim([0,4])
ylabel('Modulated EMG (% of MVC)')
xlabel('time (s)')
legend('first differeitation','amplitude threshold','threshold triggers')
title('sEMG first differentiation and threshold')



subplot(2,1,2)
area(x,FES1double*0.1,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5)
hold on
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
area(x1,FES1*0.05,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.5)
xlabel('time (s)')
ylim([0,0.3])
yticks([])
legend('double threshold detection outcome','Detection target (posture change/ sound)','single threshold detection outcome')
title('Detection Oucomes for adapted double threshold algorithm')

%% modified the diagram
ondouble = zeros(size(FES1double));
ondouble(FES1double==0)= NaN;


figure(17)
subplot(3,1,1)
rectangle('Position',[54.8,0,0.67,100],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[53,0,1.8,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,2.65,100],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
plot(test_time(1:end-9),test_bar1,'Color','k','LineWidth',2)
hold on
plot(x1,ondouble,'Color','r','LineWidth',6)
text([53.8,55,56],[3,3,3],{'Rest','Valsalva','Rest'},'FontSize',13,'Color',[150 150 150]./255,'FontSize',20)


xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,1,'amplitude threshold = 0.37','Color','r','FontSize',13)
ylim([0,4])
xlim([53,58])
ylabel('EMG 1^{st} differentiation')
xlabel('time (s)')
legend('first differeitation','Detection algo triggered')

set(gca,'FontSize',20)





%% latency and accuracy
laccur5 = xlsread('EMG recording.xlsx','X20:AB22');

plot(22)
yyaxis left

plot(laccur5(1,:),laccur5(3,:),'ro-','LineWidth',2,'MarkerFaceColor','r')
xlim([60,160])
ylim([60,102])
ylabel('Dectection Accuracy (%)','FontSize',20)
%subplot(2,1,1)
ylabel('Detection Accuracy (s)','FontSize',20)

yyaxis right
%subplot(2,1,2)
plot(laccur5(1,:),laccur5(2,:),'bs-','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','b')
xlim([60,160])
ylim([-0.15,0.05])
yline(0,'.--','Color','k')
ylabel('Mean Detection Latency (s)','FontSize',20)

set(gca,'FontSize',20)

legend('Detection accuracy','Mean detection latency','FontSize',20)
xlabel('Input IAP threshold value (cmH_2O)','FontSize',20)
title('Detection Accuracy and Mean Latency for Different Input IAP Thresholds')

set(gca,'FontSize',20)

%% actuation delay demonstration 
peak = 223;
rest = 32.5;

tao1 = 0.1;
tao2 = 1;

x1= linspace(0,2.525,100);
y1 = exp((x1-2)/tao1)+rest;
x2 = linspace(2.54,4.3,100);
y2 = peak* ones(1,size(x2,2));
x3 = linspace(4.3, 10, 500);
y3 = (peak-rest)*(exp(-1*(x3-4.3)/tao2))+rest;
x3 = (x3-4.3)/5.7;

xIAP_sample = cat(2,x1,x2);
xIAP_sample = xIAP_sample(1:124)/2.98;
yIAP_sample = cat(2,y1,y2);
yIAP_sample = yIAP_sample(77:200);
time_sample = linspace(0,10,100);
IAP_sit = ones(1,100)*32;
IAP_stand = ones(1,100)*33;


inc1= (yIAP_sample-33)/190*48.9+33;
y31 = (y3-33)/190*48.9+33;
inc2 = (yIAP_sample-33)/190*42+33;
y32 = (y3-33)/190*42+33;
inc3 = (yIAP_sample-33)/190*17.5+33;
y33 = (y3-33)/190*17.5+33;
inc4 = (yIAP_sample-33)/190*154.6+33;
y34 = (y3-33)/190*154.6+33;
inc5 = (yIAP_sample-33)/190*190+33;
y35 = (y3-33)/190*190+33;
inc6 =  (yIAP_sample-33)/190*242.4+33;
y36 = (y3-33)/190*242.4+33;
%simux = cat(2,time_sample/10*13,xIAP_sample(1:200)/10*6+13,time_sample/10*11+19,xIAP_sample/10*7.35+30,time_sample/10*6.85+37.35,xIAP_sample/10*4.94+44.2,time_sample/10*4.66+49.14,xIAP_sample/10*7+53.8,time_sample/10*0.33+60.8,xIAP_sample+61.13,time_sample/10*0+71.13,xIAP_sample/10*14+71.13,time_sample/10*10.43+85.13);
%simuy = cat(2,IAP_sit,inc1,IAP_stand,inc2,IAP_stand,inc3,IAP_stand,inc4,IAP_stand,inc5,IAP_stand,inc6,IAP_stand);

simux = cat(2,time_sample/10*13,xIAP_sample*5+13,x3*1+18,time_sample/10*11+19,xIAP_sample*7.35+30,x3+37.35,time_sample/10*5.85+38.35,xIAP_sample*4.94+44.2,x3*0.5+49.14,time_sample/10*5.16+49.64,xIAP_sample*0.67+54.8,x3*5+55.47,time_sample/10*1.66+60.47,xIAP_sample*3.99+62.13,x3*5.7+66.12,time_sample/10*1.18+71.82,xIAP_sample*2.45+73,x3*6+75.45,time_sample/10*10.43+81.45);
simuy = cat(2,IAP_sit,inc1,y31,IAP_stand,inc2,y32,IAP_stand,inc3,y33,IAP_stand,inc4,y34,IAP_stand,inc5,y35,IAP_stand,inc6,y36,IAP_stand);
yline_coor = [13,18,30,37.35,44.2,49.14,54.8,55.47,62.13,66.12,73,75.45];


actuation_target = zeros(1,size(simuy,2));
for i =1:size(simuy,2)
    if simuy(i)>36
        actuation_target(i)=1;
    end
    
end

m =250;      % consecutive 0.4s
j = m; % pointer
acuaDelay = zeros([1,size(FES1,2)]);
while(j<size(FES1,2))
    if sum(FES1((j-m+1):j))>0
        acuaDelay(j)=1;
    else
        acuaDelay(j)=0;
    end
    
    j = j+1;
    
end


figure(23)

subplot(3,1,1)
rectangle('Position',[0,0,13,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
hold on 
rectangle('Position',[13,0,5,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[18,0,12,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
rectangle('Position',[30,0,7.35,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,6.85,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
rectangle('Position',[44.2,0,4.94,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[49.14,0,5.66,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
rectangle('Position',[54.8,0,0.67,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,6.66,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
rectangle('Position',[62.13,0,3.99,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[66.12,0,6.88,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
rectangle('Position',[73,0,2.45,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[75.45,0,14.55,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
text([12,31,44,54.85,62.5,73],ones(1,6)*250,{'Sit-to-Stand','Squat','Leg raise','Valsalva','Coughs','Jump'},'FontSize',20,'Color','k')
plot(simux,simuy,'k','LineWidth',2)
hold on 
xlabel('time (s)')

%yline(70,'--','LineWidth',1,'Color',[0.8500 0.3250 0.0980])
yline(70-30,'--','LineWidth',1,'Color','r')
%xline (yline_coor,'-.','Color',[150 150 150]./255)
area(simux,actuation_target*40,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','c')
%text(1,85,'input IAP threshold = 70','Color',[0.8500 0.3250 0.0980],'FontSize',15)
%text(1,40,'IAP threshold','Color',[0.3010 0.7450 0.9330],'FontSize',20)
xlim([0,90])
%legend('Simulated IAP change','input IAP threshold','actual threshold (bias towards false positive)')
title('Intra-abdominal Pressure during Activities')
set(gca,'FontSize',20)
ylabel('simulated IAP (cmH_2O)','FontSize',15)


subplot(3,1,2)
area(x,FES1double*0.1,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5)
hold on
%plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
plot(simux,actuation_target*0.08,'c','LineWidth',1.5)
area(simux,actuation_target*0.08,'FaceColor','c','FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
%xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
xlim([0,90])
ylim([0,0.15])
yticks([0,90])
%legend('threshold detection outcome','Detection target (posture change/ sound)','Simulated incontience periods')
%title('Actuation periods without time delay (same as thresholding algorithm triggering periods)')
title('Detection Outome')
set(gca,'FontSize',20)

subplot(3,1,3)
area(x,acuaDelay*0.1,'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.5)
hold on
%plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
plot(simux,actuation_target*0.08,'c','LineWidth',1.5)
area(simux,actuation_target*0.08,'FaceColor','c','FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
%xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
xlim([0,90])
ylim([0,0.15])
yticks([0,90])
%legend('actuation periods','Detection target (posture change/ sound)','Simulated incontience periods')
%title('Actuation periods with introduced five-second time delay')
title('Stimulation Periods')
set(gca,'FontSize',20)

%% plot for ppt
figure (24)


subplot(3,1,1)
plot(test_time(1:end-9),test_bar1,'k','LineWidth',1)
hold on 
yline(threshold1,'r--','LineWidth',1)
text([5,12,23,31,45,53.5,62.5,75],[4,3,4,4,4,4,4,2],{'sit','sit-to-stand','stand','squat','leg raise','Valsalva','two coughs','jump'},'FontSize',20,'Color',[150 150 150]./255)
area(x1,FES1*threshold1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.4)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(25.2,0.5,'amplitude threshold = 0.37','Color','r','FontSize',18)
xlabel('time (s)')
legend('first differeitation','amplitude threshold','threshold triggers')
title('sEMG first differentiation and threshold')
set(gca,'FontSize',20)
ylabel('sEMG amplitude first derivative','FontSize',15)
ylim([0,3])


subplot(3,1,2)
area(x1,FES1*0.25,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.4)
hold on
plot(test_time(1:end-(W-1)),target*0.2*2.5,'r','LineWidth',1.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
ylim([0,0.3])
xlabel('time (s)')
yticks([])
legend('single threshold detection outcome','Detection target (posture change/ sound)')
title('Detection Oucomes for single threshold algorithm (only threshlding magnitude)')
set(gca,'FontSize',20)



subplot(3,1,3)
area(x,FES1double*0.25,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5)
hold on 
plot(test_time(1:end-(W-1)),target*0.2*2.5,'r','LineWidth',1.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
ylim([0,0.3])
yticks([])
legend('double threshold detection outcome','Detection target (posture change/ sound)')
title('Detection Oucomes for adapted double threshold algorithm (threshlding both magnitude and time)')
set(gca,'FontSize',20)



%% copy of plot for ppt

subplot(3,1,1)
rectangle('Position',[54.8,0,0.67,5],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[25,0,5,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,2.65,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
plot(test_time(1:end-9),test_bar1,'k','LineWidth',1)
%yline(threshold1,'r--','LineWidth',2)
%text([5,12,23,31,45,54.9,62.5,75],[4,3,4,4,4,4.2,4,2],{'sit','sit-to-stand','stand','squat','leg raise','Valsalva','two coughs','jump'},'FontSize',20,'Color','k')
%text([53,56.5],[4.2,4.2],{'Stand','Stand'},'FontSize',20,'Color','k')
%area(x1,FES1*threshold1,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.8)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%text(25.2,0.5,'Magnitude Threshold','Color','r','FontSize',18)
xlabel('time (s)')
%legend('first differeitation','amplitude threshold','threshold triggers')
set(gca,'FontSize',20)
%ylabel('sEMG amplitude first derivative','FontSize',15)
xlim([52,58])
ylim([0,5])
title('sEMG First Derivative Magnitude')
%text([26,32,38],ones(1,3)*1,{'Stand','Squat','Stand'},'FontSize',30,'Color','k')
set(gca,'XTick',52:6/7.5:58);
set(gca,'XTickLabel',0:1:7.5)

subplot(3,1,2)
rectangle('Position',[54.8,0,0.67,5],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[25,0,5,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,2.65,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
area(x1,FES1*0.25,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.8)
%plot(test_time(1:end-(W-1)),target*0.2*2.5,'r','LineWidth',1.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
ylim([0,0.4])
xlabel('time (s)')
yticks([])
%legend('single threshold detection outcome','Detection target (posture change/ sound)')
xlim([52,58])
title('Single Threshold Outcome')
set(gca,'FontSize',20)
set(gca,'XTick',52:6/7.5:58);
set(gca,'XTickLabel',0:1:7.5)


subplot(3,1,3)
rectangle('Position',[54.8,0,0.67,5],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[25,0,5,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,2.65,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on
area(x,FES1double*0.25,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.8)
%plot(test_time(1:end-(W-1)),target*0.2*2.5,'r','LineWidth',1.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
ylim([0,0.4])
yticks([])
%legend('double threshold detection outcome','Detection target (posture change/ sound)')
xlim([52,58])
title('Double Threshold Outcome')
set(gca,'FontSize',20)
set(gca,'XTick',52:6/7.5:58);
set(gca,'XTickLabel',0:1:7.5)
%%
subplot(4,1,4)
plot(simux(simux_index),simuy(simux_index),'k')
hold on 
yline(70,'--','LineWidth',1,'Color',[0.8500 0.3250 0.0980])
yline(70-30,'--','LineWidth',1,'Color',[0.3010 0.7450 0.9330])
xline (yline_coor(7:8),'-.','Color',[150 150 150]./255)
area(simux(simux_index),actuation_target(simux_index)*40,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
text(1,85,'input IAP threshold = 70','Color',[0.8500 0.3250 0.0980],'FontSize',15)
text(1,25,'actual IAP threshold = 36','Color',[0.3010 0.7450 0.9330],'FontSize',15)
xlim([25,40])
%legend('Simulated IAP change','input IAP threshold','actual threshold (bias towards false positive)')
title('IAP simulation, input IAP threshold value, and actual threshold')
set(gca,'FontSize',20)
ylabel('simulated IAP (cmH_2O)','FontSize',15)
xlabel('time (s)','FontSize',15)
ylim([32,90])
set(gca,'FontSize',20)

%% ppt picture 4

overflow = zeros(size(actuation_target));
overflow(actuation_target==0)= NaN;




figure(25)

subplot(3,1,1)
rectangle('Position',[0,0,13,300],'FaceColor',[0,1,0,0.1],'EdgeColor','none')

hold on 
rectangle('Position',[13,0,5,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[18,0,12,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[30,0,7.35,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,6.85,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[44.2,0,4.94,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[49.14,0,5.66,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[54.8,0,0.67,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,6.66,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[62.13,0,3.99,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[66.12,0,6.88,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[73,0,2.45,300],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[75.45,0,14.55,300],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')



hold on 
%yline(70,'--','LineWidth',1,'Color',[0.8500 0.3250 0.0980])
%yline(70-30,'LineWidth',1.5,'Color','m')
xline (yline_coor,'-.','Color',[150 150 150]./255)
%area(simux,actuation_target*40,'FaceColor','c','FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
IAPplot = plot(simux,simuy,'k','LineWidth',1,'DisplayName','Simulated IAP change');
overflowplot = plot(simux,overflow,'Color','c','LineWidth',6,'DisplayName','Estimated incontinece episode');
%text(1,85,'input IAP threshold = 70','Color',[0.8500 0.3250 0.0980],'FontSize',15)
text(1,65,'Target IAP Threshold','Color','r','FontSize',18)
xlim([52,61])
ylim([0,200])
text([53,55,58],[150,150,150],{'Rest','Valsalva','Rest'},'Color',[150,150,150]./255,'FontSize',20)
%legend('Simulated IAP change','input IAP threshold','actual threshold (bias towards false positive)')
%title('IAP simulation, input IAP threshold value, and actual threshold')
title('IAP Simulation')
ylabel('simulated IAP (cmH_2O)','FontSize',15)
xlabel('time (s)','FontSize',15)
legend([IAPplot,overflowplot])
set(gca,'XTick',52:1:61);
set(gca,'XTickLabel',0:1:9)
set(gca,'FontSize',20)



subplot(3,1,2)
rectangle('Position',[0,0,13,8],'FaceColor',[0,1,0,0.1],'EdgeColor','none')
hold on 
rectangle('Position',[13,0,5,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[18,0,12,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[30,0,7.35,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,6.85,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[44.2,0,4.94,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[49.14,0,5.66,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[54.8,0,0.67,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,6.66,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[62.13,0,3.99,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[66.12,0,6.88,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[73,0,2.45,8],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[75.45,0,14.55,8],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
%area(x1,FES1double*threshold1,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.8)
detectionplot = plot(x1,ondouble+0.175,'Color','r','LineWidth',6,'DisplayName','Detection algo triggered episode');
overflowplot = plot(simux,overflow,'Color','c','LineWidth',6,'DisplayName','Estimated incontinece episode');
sEMGplot = plot(test_time(1:end-9),test_bar1,'k','LineWidth',1,'DisplayName','sEMG signal');
%yline(threshold1,'m','LineWidth',2)
%text([5,12,23,31,44,54.85,62.5,73],[4,6,4,6,6,4,6,6],{'Sit','Sit-to-Stand','Stand','Squat','Leg raise','Valsalva','Coughs','Jump'},'FontSize',20,'Color','k')
%text([53,58],[4,4],{'Stand','Stand'},'FontSize',20,'Color','k')
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,1,'Magnitude Threhsold','Color','r','FontSize',18)
xlim([52,61])
ylim([0,5])
legend([sEMGplot,detectionplot,overflowplot])
%legend('first differeitation','amplitude threshold')
%title('sEMG first differentiation and threshold')
title('Thresholding Detection Outcome')
%ylabel('sEMG amplitude first derivative','FontSize',15)
xlabel('time (s)','FontSize',15)
set(gca,'XTick',52:1:61);
set(gca,'XTickLabel',0:1:9)
set(gca,'FontSize',20)

%%
subplot(3,1,3)
rectangle('Position',[0,0,13,0.15],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
hold on 
rectangle('Position',[13,0,5,0.15],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[18,0,12,0.15],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[30,0,7.35,0.15],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[37.35,0,6.85,0.15],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[44.2,0,4.94,1],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[49.14,0,5.66,1],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[54.8,0,0.67,1],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[55.47,0,6.66,1],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[62.13,0,3.99,1],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[66.12,0,6.88,2],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
rectangle('Position',[73,0,2.45,2],'FaceColor',[1,0,0,0.1],'EdgeColor','none')
rectangle('Position',[75.45,0,14.55,0.15],'FaceColor',[0.4660 0.6740 0.1880,0.1],'EdgeColor','none')
xline (yline_coor,'-.','Color',[150 150 150]./255)
%area(x,FES1double*0.1,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',1)
hold on
plot(x1,ondouble+0.175,'Color','r','LineWidth',6)
plot(simux,overflow,'Color','c','LineWidth',6)
%plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
%plot(simux,actuation_target*0.08,'c','LineWidth',2)
%area(simux,actuation_target*0.08,'FaceColor','c','FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])

%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
xlim([52,61])
ylim([0,0.5])
yticks([0,90])
%legend('threshold detection outcome','Detection target (posture change/ sound)','Simulated incontience periods')
%title('Actuation periods without time delay (same as thresholding algorithm triggering periods)')
title('Detection Outcome & Incontinence Episodes')
set(gca,'FontSize',20)
set(gca,'XTick',52:1:61);
set(gca,'XTickLabel',0:1:9)

%% ppt picutre 5
time_index = find (test_time >26 & test_time <40);
simux_index = find (simux >26 & simux <40);
x_index = find (x >26 & x <40);

figure(26)

subplot(3,1,2)
plot(simux(simux_index),simuy(simux_index),'k')
hold on 
yline(70,'--','LineWidth',1,'Color',[0.8500 0.3250 0.0980])
yline(70-30,'--','LineWidth',1,'Color',[0.3010 0.7450 0.9330])
xline (yline_coor(7:8),'-.','Color',[150 150 150]./255)
area(simux(simux_index),actuation_target(simux_index)*40,'FaceColor',[0.3010 0.7450 0.9330],'FaceAlpha',0.3,'EdgeColor',[0.3010 0.7450 0.9330])
text(1,85,'input IAP threshold = 70','Color',[0.8500 0.3250 0.0980],'FontSize',15)
text(1,25,'actual IAP threshold = 36','Color',[0.3010 0.7450 0.9330],'FontSize',15)
xlim([26,40])
legend('Simulated IAP change','input IAP threshold','actual threshold (bias towards false positive)')
title('IAP simulation, input IAP threshold value, and actual threshold')
set(gca,'FontSize',20)
ylabel('simulated IAP (cmH_2O)','FontSize',15)
xlabel('time (s)','FontSize',15)


subplot(3,1,1)
plot(test_time(1:end-9),test_bar1,'k','LineWidth',1)
hold on 
yline(threshold1,'r--','LineWidth',1)
text([5,12,23,31,45,53.5,62.5,75],[4,3,4,4,4,4,4,2],{'sit','sit-to-stand','stand','squat','leg raise','Valsalva','two coughs','jump'},'FontSize',20,'Color',[150 150 150]./255)
xline (yline_coor,'-.','Color',[150 150 150]./255)
text(1,1,'amplitude threshold = 0.37','Color','r','FontSize',13)
xlim([26,40])
legend('first differeitation','amplitude threshold')
title('sEMG first differentiation and threshold')
set(gca,'FontSize',20)
ylabel('sEMG amplitude first derivative','FontSize',15)
xlabel('time (s)','FontSize',15)

subplot(3,1,3)
area(x,FES1double*0.1,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.5)
hold on
plot(test_time(1:end-(W-1)),target*0.2,'r','LineWidth',1.5)
plot(simux,actuation_target*0.05,'c','LineWidth',1.5)
xline (yline_coor,'-.','Color',[150 150 150]./255)
%plot(switch_time, y2,'Color',[0.4660 0.6740 0.1880],'LineWidth',8)
xlabel('time (s)')
xlim([26,40])
ylim([0,0.3])
yticks([0,90])
legend('threshold detection outcome','Detection target (posture change/ sound)','Simulated incontience periods')
title('Actuation periods without time delay (same as thresholding algorithm triggering periods)')
set(gca,'FontSize',20)
