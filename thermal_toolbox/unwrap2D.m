function phi_uw = unwrap2D(phi_w,idx)

[M N] = size(phi_w);

x_wrap = 0;

% Unwrap a row of phi
merit = unwrap(phi_w(fix(M/2)+1-x_wrap,:));
% Unwrap the columns of phi
for j = 1:N
    col(:,j) = unwrap(phi_w(:,j));
    delta(j) = merit(j)- col(fix(M/2)-x_wrap,j);
    phi_uw(:,j) = col(:,j) + delta(j)*ones(M,1);
end

phi_uw = phi_uw;% - mean(phi_uw(idx==1))*ones(N,N);
end