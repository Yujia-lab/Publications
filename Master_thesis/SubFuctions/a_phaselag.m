function [GS_phase,ROISignals_phase] = a_phaselag(GS,ROISignals)
GS = squeeze(GS);
ROISignals = permute(ROISignals,[2,1,3]);
GS_phase = angle(hilbert(GS));
ROISignals_phase = angle(hilbert(ROISignals));

% phase_lag = zeros(1200,246,size(ROISignals,3));
% for i = 1:size(ROISignals,3)
%     phase_lag(:,:,i) = 1 - abs(sin(bsxfun(@minus,ROISignals_phase(:,:,i),GS_phase(:,i))));
% end

