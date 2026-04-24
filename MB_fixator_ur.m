function [muRel, M] = MB_fixator_ur(B)
% Express the relationship between H and relative permeability of the fixator
% B: magnetic flux density vector [Bx, By, Bz]
% muRel: relative permeability
% M: magnetization vector

mu0 = 4 * pi * 1e-7;

H_data = [0    270       400       800       1600      8000      16000     30000];
B_data = [0    0.3000    0.6000    1.0000    1.2000    1.3000    1.3101    1.3276];

H_output=zeros(3,1);

for ind=1:3
    H_output(ind) = interp1(B_data, H_data, abs(B(ind)), 'linear');

    % Symmetric extension for negative B
    if B(ind)<0
        H_output(ind) = -H_output(ind);
    end
end

muRel = norm(B) ./ (norm(H_output) * mu0);

M = (muRel - 1).*H_output;

end