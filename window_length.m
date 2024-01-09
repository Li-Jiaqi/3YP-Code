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

f_low= 300; % attenuate low-frequency noise such as breathing
f_high = 40; 
norder = 4;
d_bp=design(fdesign.bandpass('N,F3dB1,F3dB2',norder,f_high,f_low,fs),'butter');
[b_bp,a_bp] = tf(d_bp);

RA_filtered = filtfilt(b_bp,a_bp,RA);
RA_filtered = abs(RA_filtered);
Oblique_filtered = filtfilt(b_bp,a_bp,Oblique);
Oblique_filtered = abs(Oblique_filtered);

segment1_unsmoothed = Oblique_filtered(173403:237008);
r = size( segment1_unsmoothed,1);

window_list = [200 400 600 800 1000 1200 1400 1600 1800 2000];
overlap_list = [100 300 500 700 900 1100 1300 1500 1700 1900];
SNR = zeros(1,size(window_list,2));
k =1;
time = zeros (1, size(window_list,2));
while k<=size(window_list,2)
    window_size = window_list(k);
    overlap_size = overlap_list(k);
    window_step = window_size - overlap_size;
    segment1_window = zeros(1,floor(r/(window_size -overlap_size)));
    RA_window = zeros(1,floor(r/(window_size -overlap_size)));
    i = 1;
    n_start = 1;
   
    tic
    while (n_start+window_size<=r)
        segment1_window(i) =sqrt(mean(segment1_unsmoothed(n_start:n_start+window_size-1).^2));
        n_start = n_start + window_size - overlap_size;
        i = i+1;
    end
    time(k)=toc;
    window_mean = mean(segment1_window);
    window_std= std(segment1_window);
    SNR(k) = (window_mean/window_std)^2;
    k = k+1;
end
figure(18)
plot(window_list/10, SNR,'k','LineWidth',1)
xlabel('window length (ms)')
ylabel('SNR = s^2/\sigma^2')
title('SNR against window length, for T_{inc} = 10ms')
xline(80,'--','Color',[150 150 150]./255)
text(85,6.5,'chosen window length = 80 ms','Color',[150 150 150]./255,'FontSize',15)
set(gca,'FontSize',15)