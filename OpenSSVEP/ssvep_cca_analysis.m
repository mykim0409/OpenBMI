% function output = ssvep_cca_analysis(~)
OpenBMI('C:\Users\cvpr\Desktop\OpenBMI') % Edit the variable BMI if necessary
global BMI;
% BMI.EEG_DIR=['C:\Users\cvpr\Desktop\Data\StarLab Data\StarlabDB\subject3 정지훈\Session2'];
% file=fullfile(BMI.EEG_DIR, '\161031_jhjeong_s2_ssvep2');
BMI.EEG_DIR=['C:\Users\cvpr\Desktop\받아랏\data'];
file=fullfile(BMI.EEG_DIR, '\hjpark_offline1');

fRange = [1 40];
fSample = 100;

% marker={'1','up', '7.5';'2','left', '10';'3','center', '12';'4','right', '15';'5','down', '20'};
marker={'1','up', '9';'2','left', '11';'3','center', '13';'4','right', '15';'5','down', '17'};
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

[seg, nTrials, ~] = size(SMT.x);
[nClasses, ~] = size(SMT.class);
selFreq = str2double(SMT.class(:,3));

t= 0.001:0.001:T;

Y = zeros(length(t),4,nClasses);
for i = 1:nClasses
    ref = 2*pi*selFreq(i)*t;
    Y(:,:,i) = [sin(ref)',cos(ref)',sin(ref*2)', cos(ref*2)'];
end

r =zeros(nTrials, nClasses, nClasses-1);

for i = 1:nTrials
    for j = 1:nClasses
        [~,~,r(i,j,:)] = canoncorr(SMT.x(:,i),Y(:,:,j));
    end
end

corrCount = zeros(nClasses,1);

% % % % for i = 1:nTrials
% % % %     for j = 1:nClasses
% % % %         if(max(r(i,:,:))) == max(r(i,j,:))
% % % %             corrCount(j) = corrCount(j) + 1;
% % % %             disp(sprintf('클래스 %d r값: %f',j, max(r(i,j,:))));
% % % %         end
% % % %     end
% % % % end

for i = 1:nTrials
    for j = 1:nClasses
        if(max(r(i,:,:))) == max(r(i,j,:))
            corrCount(j) = corrCount(j) + 1;
            disp(sprintf('클래스 %d r값: %f',j, max(r(i,j,:))));
        end
    end
end

Accuracy = sum(corrCount)/length(CNT.t)*100;