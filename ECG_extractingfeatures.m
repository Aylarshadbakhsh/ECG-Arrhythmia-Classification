clc
clear
close all


load 106
load 106_atr

fs = 360;
n = 2;
T=10;
Fea = [];
Y = [];

for i = 1:numel(ANNOT)-1 %RR intervals 
    
    ix1 = round(ATRTIME(i)*fs);
    ix2 = round(ATRTIME(i+1)*fs-1);
    ecg = val(1,ix1:ix2);
    t = 1/fs:1/fs:(ix2-ix1+1)/fs;

    if (ANNOT(i) == 1)||(ANNOT(i) == 5)

        smoth_ecg = sgolayfilt(ecg,5,21);


        C = polyfit(t,smoth_ecg,7);
        BaseLine = polyval(C,t);


        denoised = smoth_ecg - BaseLine;
        x = denoised;

        % Feature Extraction 
        % Statistical    
        avg = mean(x);
        stdev = std(x);
        skew = skewness(x);
        kurt = kurtosis(x);

        Stats = [avg stdev skew kurt];

        % Freq. Features
        A = 2*abs(fft(x,3600))/3600;
        f = fs/2*linspace(0,1,3600/2+1);
        A = A(1:3600/2+1);
        [amp, dominantFreq] = findpeaks(A,f,'NPeaks',3,'SortStr','descend');

        Freq = [amp dominantFreq];

        % Wavelet
        wave_name = 'db10';
        [C,L] = wavedec(x,4,wave_name);
        App4 = appcoef(C,L,wave_name,4);
        Det4 = detcoef(C,L,4);
        Det3 = detcoef(C,L,3);
        Det2 = detcoef(C,L,2);
        Det1 = detcoef(C,L,1);

        WaveletFeatures = [sum(App4.^2)/numel(App4) sum(Det4.^2)/numel(Det4) ...
                           sum(Det3.^2)/numel(Det3) sum(Det2.^2)/numel(Det2) ...
                           sum(Det1.^2)/numel(Det1)];

        % Nonlinear
        ShannonEntropy = wentropy(x,'shannon');

        F = [Stats, Freq, WaveletFeatures, ShannonEntropy];

        Fea = [Fea; F];
        
        Y = [Y; (ANNOT(i)-1)/4];
    end

end

save('ExtractedFeatures','Fea','Y')
