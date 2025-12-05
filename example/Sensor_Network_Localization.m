%% Generate the sensor network localization problem
rng(1);
n = 10; % number of sensors
loc = rand(2, n);
Aset = [];
a = [0.25 0.75 0.3 0.8; 0.75 0.25 0.8 0.3]; % anchors
for i = 1:n-1
    for j = i+1:n
        if norm(loc(:,i) - loc(:,j))^2 <= 0.5
            Aset = [Aset [i;j]];
        end
    end
end
Bset = [];
for i = n
    for j = 1:size(a,2)
        if norm(loc(:,i) - a(:,j))^2 <= 0.5
            Bset = [Bset [i;j]];
        end
    end
end

x = msspoly('x', 2*n);
f = 0;
for i = 1:size(Aset,2)
    f = f + ((x(Aset(1,i))-x(Aset(2,i)))^2 + (x(Aset(1,i)+n)-x(Aset(2,i)+n))^2 - (loc(1,Aset(1,i))-loc(1,Aset(2,i)))^2 - (loc(2,Aset(1,i))-loc(2,Aset(2,i)))^2)^2;
end
for i = 1:size(Bset,2)
    f = f + ((x(Bset(1,i))-a(1,Bset(2,i)))^2 + (x(Bset(1,i)+n)-a(2,Bset(2,i)))^2 - (loc(1,Bset(1,i))-a(1,Bset(2,i)))^2 - (loc(2,Bset(1,i))-a(2,Bset(2,i)))^2)^2;
end

%% Generate moment relaxations
cliques{1} = 1:2*n;
[At, b, c, K] = snl_mom_sparse(f, x, cliques);

%% Solve with ManiSDP
rng(0);
clear options;
options.tol = 1e-4;
options.sigma0 = 1;
options.sigma_min = 1e1;
options.theta = 1e-3;
options.TR_maxiter = 8;
options.line_search = 0;
options.alpha = 0.01;
maxc = max(abs(c));
tic
[~, fval, data] = ManiSDP(At, b, c/maxc, K, options);
fval = fval*maxc;
emani = max([data.gap, data.pinf, data.dinf]);
tmani = toc;


%% Generate SOS relaxations
% [blk, At, C, b] = sosprogram(-1, [f, -1], [], [], x, 2);
% [At, b, c, K] = SDPT3data_SEDUMIdata(blk, At, C, b);

%% Solve with ManiDSDP
% rng(0);
% clear options;
% options.tol = 1e-4;
% options.sigma0 = 1;
% options.theta = 1e-3;
% options.TR_maxiter = 6;
% options.tau1 = 1e-2;
% maxb = full(max(abs(b)));
% tic
% [~, dfval, data] = ManiDSDP(At', b/maxb, -c, K, options);
% dfval = dfval*maxb;
% emanid = max([data.gap, data.pinf, data.dinf]);
% tmanid = toc;

%% Solve with MOSEK
% prob       = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
% prob       = convert_sedumi2mosek(At, b, c, K);
% tic
% [~,res]    = mosekopt('minimize echo(3)',prob);
% [X,y,S,mobj] = recover_mosek_sol_blk(res, blk);
% by = b'*y;
% gap = abs(mobj(1)-by)/(abs(by)+abs(mobj(1))+1);
% x = zeros(sum(K.s.^2)+K.f, 1);
% x(1:K.f) = X{end};
% mS = zeros(length(K.s), 1);
% ind = K.f+1;
% for i = 1:length(K.s)
%     x(ind:ind+K.s(i)^2-1) = X{i}(:);
%     ind = ind + K.s(i)^2;
%     [~, dS] = eig(S{i}, 'vector');
%     mS(i) = max(0, -dS(1))/(1+abs(dS(end)));
% end
% eta = norm(At'*x - b)/(1+norm(b));
% emosek = max([eta, gap, max(mS)]);
% tmosek = toc;


% fprintf('Mosek: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', mobj(1), emosek, tmosek);
% fprintf('SDPLR: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', vlr, elr, tlr);
% fprintf('SDPNAL: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', objnal(1), enal, tnal);
fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
% fprintf('ManiDSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', dfval, emanid, tmanid);