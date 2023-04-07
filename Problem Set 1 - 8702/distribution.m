%Allocate array for the stationary distribution
Pi = zeros(NB, NSS);

BP_vec = unique(bpSIM);
bp_vec = repmat(BP_vec,1,NSS);

%Compute the stationary distribution 
for i_b=1:NB
    for i_y=1:NSS
        logical = (S_index==i_y) .* (bpSIM==BP_vec(i_b));
        Pi(i_b, i_y) = sum(logical(:)) / T;
    end
end


