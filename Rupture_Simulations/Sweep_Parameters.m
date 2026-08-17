% generate space of k0 and delat parameters
K0space = logspace(-1,1,10);
d_space = logspace(-0.7,0.5,9);

% Output RMSE for three proteins
R1 = zeros(length(K0space),length(d_space));
R2 = zeros(length(K0space),length(d_space));
R3 = zeros(length(K0space),length(d_space));

DD = [];
KK = [];

for i1=1:length(K0space)
    KK(i1) = mean(K0space(i1)*[2.7e-3 2.9e-3]);
    for j1=1:length(d_space)
        K0c = K0space(i1);
        dc = d_space(j1);
        DD(j1) = mean(dc*[1.23 1.63]);
        Model_Distribution;
        R1(i1,j1) = RMSEcond;
        R2(i1,j1) = RMSEcoh;
        R3(i1,j1) = RMSE_h_c;
    end
end


% [i,j]=min(R2,[],"all");
% [i1,i2]=ind2sub(size(R1),j);
% 
% K0space(i1)*[2.7e-3 2.9e-3]
% d_space(i2)*[1.23 1.63];