function R_correct = a_multicorrect(R,P,mode,thre)
% three method for multiple comparison correction
if mode == 'FDR'
    P_correct = mafdr(P,'BHFDR','true');

end
if mode == 'FWE'
    P_correct = P .* length(P);
end

if mode == 'NOC'
    P_correct = P;
end

FDR_pos = find(P_correct > thre);
R(FDR_pos) = 0;
R_correct = R;