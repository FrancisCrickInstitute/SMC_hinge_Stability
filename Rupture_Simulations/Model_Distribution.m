
F_Pombe = [25	16	10	36	12	35	10	10	20	20	45	14	14	27	18	17	5.9	17	6	15	24	24	14	22	27	34	16	17	7.6	4.2	3.6	16	20	8	6.1	6.4	11	16	12	20	16];
F_Human = [29    18    18    37    20    19    30    20    38    33    18    20    28    10    16    40    70    20    63    35    45    20    21];

F_pombe_condensin = [20 8.6 17 10 15 12 15 5.5 14 23 16 24 9.7 19 18 23 16 16 9.5 13 9 21 20 22 16 12 16];

global L10 L20 D x0

v = 0.16;   % loading rate (microns/second)
L20 = 16.32;  % first DNA contour length (microns), persistence length is assumed 50 nm
L10 = 8;   % second DNA contour length (microns), persistence length is assumed 50 nm
Ntrials = round(250/2);
T = zeros(Ntrials,1);
Fd = zeros(Ntrials,1);
Kd = zeros(Ntrials,1);
N = 1350;
dt = 0.05;           % time step (seconds)

% Start Parameter Set
k0 = K0c*[2.7e-3 2.9e-3];
delta = dc*[1.23 1.63];


Nmol = 1; % number of cohesin molecules
NmolActive = Nmol; % number of molecules not broken

F = zeros(N,1);
L = zeros(N,1);
ypos = zeros(N,1);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');  

% de = (randi(2,Ntrials,1)-1)*0.3+1.3;


tic
hw = waitbar(0,'Please wait... ');

ee = rand(Ntrials,1)*1.6+8;   % this is from experimental distribution measured by Martina
ind = randi(2,Ntrials,1);
for j=1:Ntrials
    r = rand(N,Nmol);
    x = 6;
%     ee = rand(N,1)*1.6+8;   % this is from experimental distribution measured by Martina
    ddx = zeros(N,1);      % displacement with respect to the center of DNA
%     delta = de(j);
    
    params = zeros(N,4);
    startt = [0 ee(j) 0 0];
    NmolActive = Nmol; % number of molecules not broken
    for i=1:N
        if ddx(i)>ee(j)/3
            ddx(i)=ee(j)/3;
        end
        x = x + dt*v;       % new position
        ypos(i) = x;

        
        fun = @findFandL;
        D = ee(j);
        x0 = x;
        t0 = startt;
        t = fsolve(fun,t0,options);
        params(i,:) = t;
        startt = t;
        F(i) = t(3);
        
        % nf = 1:NmolActive;
        % fdr = exp(-nf);
        % fdr = fdr*F(i)/sum(fdr);
        for q = 1:NmolActive
            CF = (F(i)/NmolActive);
            % CF = fdr(q);
            k = k0(ind(j))*exp(CF*delta(ind(j))/4.14);
            p = -log(1-r(i,q))/k;
            if p<= dt
                NmolActive = NmolActive-1;
                % L20 = L20 + 0.025;
                % L10 = L10 + 0.02;
                % x=x-0.005;
            end
        end
        if NmolActive<1
            break
        end

    end
%     plot(ypos(1:5000),F(1:5000))
    T(j) = i*dt;
    Fd(j) = F(i);
    % Kd(j) = k0(ind(i));
        waitbar(j/Ntrials,hw);
        
        % plot(T,Fd,'.')
        % pause
end
toc
close(hw)

% [p,h,stats]=ranksum(F_Pombe,F_Human);

% figure
% h = histogram(Fd,25);
% emax = round(ceil(max(h.BinEdges)),1,'significant');
% emax = max(emax,50);
% title(['mean force: ' num2str(mean(Fd)) ' pN; d=' num2str(delta) '; k0 = ' num2str(k0)])

% e = 0:4:emax;
% n = histc(Fd,e);
% n = n/max(n);
% nPb = histc(F_Pombe,e);
% nPb = nPb/max(nPb);
% nHu = histc(F_Human,e);
% nHu = nHu/max(nHu);
% nCo = histc(F_pombe_condensin,e);
% nCo = nCo/max(nCo);
% figure
% hold on
% plot(e,n)
% plot(e,nPb)
% plot(e,nHu)
% plot(e,nCo)
% legend('simulated','Pombe','Human','Condensin')
% title(['k0=' num2str(k0) '; d=' num2str(delta) '; Nmol = ' num2str(Nmol)])

e = 0:1:emax;
n = histc(Fd,e);
e = e';
e = e(4:end);
n = n(4:end);
n = n/sum(n);
m = sum(e.*n);
s = sqrt(sum(n.*((e-m).^2)));

% disp(['Pombe: ' num2str(mean(F_Pombe)) ' +/- ' num2str(std(F_Pombe)) '; Human: ' num2str(mean(F_Human)) ' +/- ' num2str(std(F_Human)) '; Current: ' num2str(m) ' +/- ' num2str(s)])


RMSEcond = sqrt(((m-mean(F_pombe_condensin))^2+(s-std(F_pombe_condensin))^2)/2);
RMSEcoh = sqrt(((m-mean(F_Pombe))^2+(s-std(F_Pombe))^2)/2);

RMSE_h_c = sqrt(((m-mean(F_Human))^2+(s-std(F_Human))^2)/2);

