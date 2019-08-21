function [Qout] = updateQ(U, F)
%%%%%%%%%%%%%%
%% solve alpha ||U - FQ||^2_F, s.t. QQ^T = I
%%%%%%%%%%%%%%

[Omega1, ~, Omega2] = svd(F' * U);
Qout =  Omega1 * Omega2';
end