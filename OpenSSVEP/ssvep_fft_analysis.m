function out = ssvep_fft_analysis(~)

OpenBMI('C:\Users\cvpr\Desktop\OpenBMI') % Edit the variable BMI if necessary
global BMI;
BMI.EEG_DIR=['C:\Users\cvpr\Desktop\Data\StarLab Data\StarlabDB\subject3 Á¤ÁöÈÆ\Session2'];
file=fullfile(BMI.EEG_DIR, '\161031_jhjeong_s2_ssvep2');
% BMI.EEG_DIR=['C:\Users\cvpr\Desktop\¹Þ¾Æ¶ù\data'];
% file=fullfile(BMI.EEG_DIR, '\hjpark_offline1');

fRange = [1 40];
fSample = 100;

marker={'1','up', '7.5';'2','left', '10';'3','center', '12';'4','right', '15';'5','down', '20'};
% marker={'1','up', '9';'2','left', '11';'3','center', '13';'4','right', '15';'5','down', '17'};
[EEG.data, EEG.marker, EEG.info]=Load_EEG(file,{'device','brainVision';'marker', marker});
% [EEG.data, EEG.marker, EEG.info]=Load_EEG(file,{'device','brainVision'});

field={'x','t','fs','y_dec','y_logic','y_class','class', 'chan'};
CNT=opt_eegStruct({EEG.data, EEG.marker, EEG.info}, field);

T = 5;
window=T*CNT.fs-1;
tInterval = [0 window];

% CNT = prep_selectChannels(CNT, {'Name',{'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2'}});
CNT = prep_selectChannels(CNT, {'Name',{'Oz'}});
% CNT = prep_filter(CNT, {'frequency', fRange; 'fs',fSample});
SMT=prep_segmentation(CNT, {'interval',tInterval});

%%
%Trigger Á¤º¸: up: 7.5, left: 10, center: 12, right: 15, down: 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%my CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[seg, nTrials, ~] = size(SMT.x);
[nClasses, ~] = size(SMT.class);
selFreq = str2double(SMT.class(:,3));
amp = zeros(nClasses,1);
% myClass = zeros(5, 10);
% myAvg = zeros(500,5);
figure;
% % % % % % for i = 1:nClasses
% % % % % %     myClass(i,:)=find(SMT.y_logic(i,:));
% % % % % %
% % % % % %     for j = 1:500
% % % % % %         myAvg(j,i) = mean(SMT.x(j,myClass(i,:)));
% % % % % %     end
% % % % % %     [YfreqDomain,frequencyRange] = positiveFFT(myAvg(:,i),SMT.fs);
% % % % % %
% % % % % %     try
% % % % % %         q = find(frequencyRange == selFreq(i));
% % % % % %     catch
% % % % % %         q = 0;
% % % % % %     end
% % % % % %
% % % % % %     if(q)
% % % % % %          amp(i) = YfreqDomain(q);
% % % % % %     else
% % % % % %         a = selFreq(i)-floor(selFreq(i));
% % % % % %         w = find(frequencyRange>7.5, 1);
% % % % % %         amp(i) = (YfreqDomain(w) + YfreqDomain(w-1)) * a;
% % % % % %     end
% % % % % %
% % % % % %     bar(1:nClasses,abs(amp))
% % % % % %
% % % % % % %     fft poower range
% % % % % % %     power - 50 trials
% % % % % % %     average_class1 = ;
% % % % % % %     2 = ;
% % % % % % %
% % % % % % %     plot();
% % % % % %
% % % % % %
% % % % % % %     plot(frequencyRange,abs(YfreqDomain)); hold on;
% % % % % % end
myAvg1 = zeros(seg,1);
cor = 0;
q=0;
for i = 1:nTrials
    myAvg1 = SMT.x(:,i);
    
    [YfreqDomain,frequencyRange] = positiveFFT(myAvg1,SMT.fs);
    
    for j = 1:nClasses
        try
            q = find(frequencyRange == selFreq(j));
        catch
            q = 0;
        end
        
        if(q)
            amp(j) = YfreqDomain(q);
        else
            w = find(frequencyRange>selFreq(j), 1);
            amp(j) = interp1(frequencyRange, YfreqDomain, selFreq(j));
        end
    end
    [~, ind] = max(abs(amp));
    SMT.class(ind,2)
    cor = strcmp(SMT.class(ind,2), SMT.y_class(i)) + cor;
    
    suptitle(sprintf('true lable: %s, predict label: %s',SMT.y_class{i}, SMT.class{ind,2}));
    
    subplot(2,1,1); bar(1:nClasses,abs(amp));
    ax = gca;
    set(ax, 'xticklabels', SMT.class(:,2));
    labels = arrayfun(@(value) num2str(value,'%2.1f'),abs(amp),'UniformOutput',false);
    text(1:nClasses,abs(amp),labels,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    ylim([0 4]);drawnow;
    subplot(2,1,2); plot(frequencyRange,abs(YfreqDomain));
    xlim([7 23]);
end
loss = cor / i *100
visual_fft(SMT,{'channel','Oz';'xlim',[5 23]})


function [X,freq] = positiveFFT(x,Fs)
N=length(x); %get the number of points
k=0:N-1;     %create a vector from 0 to N-1
T=N/Fs;      %get the frequency interval
freq=k/T;    %create the frequency range
X=fft(x)/N*2; % normalize the data

%only want the first half of the FFT, since it is redundant
cutOff = ceil(N/2);

%take only the first half of the spectrum
X = X(1:cutOff);
freq = freq(1:cutOff);