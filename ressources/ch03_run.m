function R = ch03_run()
%CH03_RUN  Numerical study of Chapter 3: active-estimation observer
%          (Spica & Robuffo Giordano, CDC 2013) for the depth of a point.
%
%   A. Reproduction of the reference script (Omega built on the estimate,
%      H = (2 sqrt(k)/sigma) I) and comparison with the observer of the paper.
%   B. Depth transient for four damping choices: reference script, script H
%      with measured s, critical damping of the paper, constant damping;
%      predicted vs measured decay rates (mass-spring-damper model).
%   C. Active estimation (paper, Sect. V-A): constant ||v||, sigma^2
%      maximised (case I) or not (case II), omega used to keep the point at
%      the image centre; camera pose integrated on SE(3).
%
%   Outputs: ../fig/ch03_*.csv, ../fig/ch03_numbers.tex,
%   code/results/ch03_results.mat. ode45 only (no Simulink, no toolbox).

here   = fileparts(mfilename('fullpath'));
figDir = fullfile(here, '..', 'fig');
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% A/B: trajectory of the reference script -------------------------------
u  = @(t) [0.7*cos(2*pi*t); 0.5*cos(pi*t/2); cos(pi*t); 0.1*cos(pi*t/2); 0; 0.1];
k  = 10;
q0 = [0; 0; 0.25];
TA = 20;  tA = (0:0.01:TA).';

prRef   = struct('k', k, 'mode', 'reference');
prPaper = struct('k', k, 'mode', 'paper');
A1 = simObs(prRef,   u, q0, [1; 2; 0.25], tA, opts);
A2 = simObs(prPaper, u, q0, [1; 2; 0.25], tA, opts);
fprintf('\n=== Chapter 3 ===\nA: reference |e(T)| = %.3e, paper |e(T)| = %.3e\n', ...
        norm(A1.q(end, :) - A1.qh(end, :)), norm(A2.q(end, :) - A2.qh(end, :)));

% B: depth transient (image initialised exactly, depth wrong)
qh0 = [0; 0; 0.5];
TB  = 10;  tB = (0:0.01:TB).';
sig = zeros(numel(tB), 1);
base = simObs(prPaper, u, q0, qh0, tB, opts);
for i = 1:numel(tB)
    ui = u(tB(i));
    sig(i) = norm([base.q(i, 1)*ui(3) - ui(1), base.q(i, 2)*ui(3) - ui(2)]);
end
ellMean = mean(2*sqrt(k)*sig);
modes = {'reference', struct('k', k, 'mode', 'reference'),    'Ref'; ...
         'refH',      struct('k', k, 'mode', 'refH'),         'RefH'; ...
         'paper',     struct('k', k, 'mode', 'paper'),        'Paper'; ...
         'const',     struct('k', k, 'mode', 'const', 'ell', ellMean), 'Const'};
nm = size(modes, 1);
eB = zeros(numel(tB), nm);  rMeas = nan(1, nm);  tSet = nan(1, nm);
for j = 1:nm
    r = simObs(modes{j, 2}, u, q0, qh0, tB, opts);
    eB(:, j) = abs(r.q(:, 3) - r.qh(:, 3));
    rMeas(j) = envelopeRate(tB, eB(:, j));
    tSet(j)  = settling(tB, eB(:, j), 0.01*eB(1, j));
end
% predicted (frozen-time mass-spring-damper) rates averaged along the run
ellRef = 2*sqrt(k)./sig;  ellPap = 2*sqrt(k)*sig;
rRef   = msdRate(ellRef, k*sig.^2);
rPap   = msdRate(ellPap, k*sig.^2);
rCst   = msdRate(ellMean*ones(size(sig)), k*sig.^2);
rPred  = [mean(rRef), mean(rRef), mean(rPap), mean(rCst)];
for j = 1:nm
    fprintf('B: %-9s rate pred %.3f  meas %.3f  t1%% %.2f s\n', modes{j, 1}, rPred(j), rMeas(j), tSet(j));
end

% B2: gain sweep -- the averaged model assumes a slowly rotating Omega;
% its validity is governed by the ratio damping / rotation rate of Omega.
kList = [10 100 1000];  sw = struct();
th = zeros(numel(tB), 1);
for i = 1:numel(tB)
    ui = u(tB(i));
    th(i) = atan2(base.q(i, 2)*ui(3) - ui(2), base.q(i, 1)*ui(3) - ui(1));
end
thd = [0; diff(unwrap(th))]/(tB(2) - tB(1));
for j = 1:numel(kList)
    kk = kList(j);
    rp = simObs(struct('k', kk, 'mode', 'paper'),     u, q0, qh0, tB, opts);
    rr = simObs(struct('k', kk, 'mode', 'reference'), u, q0, qh0, tB, opts);
    sw.measP(j) = envelopeRate(tB, abs(rp.q(:, 3) - rp.qh(:, 3)));
    sw.measR(j) = envelopeRate(tB, abs(rr.q(:, 3) - rr.qh(:, 3)));
    sw.predP(j) = mean(msdRate(2*sqrt(kk)*sig,  kk*sig.^2));
    sw.predR(j) = mean(msdRate(2*sqrt(kk)./sig, kk*sig.^2));
    sw.hOverTh(j) = median(2*sqrt(kk)*sig./max(abs(thd), 1e-6));
    fprintf('B2: k = %5d  paper meas %.2f pred %.2f | reference meas %.2f pred %.2f | median h/|thdot| %.2f\n', ...
            kk, sw.measP(j), sw.predP(j), sw.measR(j), sw.predR(j), sw.hOverTh(j));
end

%% C: active estimation --------------------------------------------------
kA  = 2000;  v0 = [0.03; 0; 0.04];  vn = norm(v0);
act.k1 = 10;  act.kdes = 0.5*vn^2;  act.lc = 1;  act.k = kA;
qA0  = [0; 0; 1];  qhA0 = [0; 0; 0.2];
TC = 10;  tC = (0:0.01:TC).';
z0 = [zeros(3, 1); reshape(eye(3), 9, 1); qA0; qhA0; v0];
act.k2 = 10;  [~, ZI]  = ode45(@(t, z) activeRhs(t, z, act), tC, z0, opts);
act.k2 = 0;   [~, ZII] = ode45(@(t, z) activeRhs(t, z, act), tC, z0, opts);
CI = activePost(ZI);  CII = activePost(ZII);
e0 = qA0(3) - qhA0(3);
sdI  = sqrt(kA)*vn;                      % case I: sigma -> ||v||
sdII = sqrt(kA)*mean(CII.sig(end-100:end));
etaI  = e0*(1 + sdI*tC).*exp(-sdI*tC);
etaII = e0*(1 + sdII*tC).*exp(-sdII*tC);
fprintf('C: sigma^2 final  I %.3e  II %.3e   (||v||^2 = %.3e)\n', CI.sig(end)^2, CII.sig(end)^2, vn^2);
fprintf('C: |e_chi(2 s)|  I %.3e  II %.3e\n', abs(CI.echi(201)), abs(CII.echi(201)));

%% Export ------------------------------------------------------------------
fl = 1e-12;
writeCsv(fullfile(figDir, 'ch03_ref.csv'), ...
    {'t','x','y','chi','xh','yh','chih','xhp','yhp','chihp'}, [tA, A1.q, A1.qh, A2.qh]);
writeCsv(fullfile(figDir, 'ch03_modes.csv'), ...
    {'t','eref','erefh','epaper','econst','sigma','href','hpaper'}, ...
    [tB, max(eB, fl), sig, ellRef, ellPap]);
writeCsv(fullfile(figDir, 'ch03_rot.csv'), {'t','thdot','h10','h1000'}, ...
    [tB, abs(thd), 2*sqrt(10)*sig, 2*sqrt(1000)*sig]);
writeCsv(fullfile(figDir, 'ch03_active.csv'), ...
    {'t','s2I','s2II','eI','eII','etaI','etaII','cxI','cyI','czI','cxII','cyII','czII', ...
     'xI','yI','vzI','vzII'}, ...
    [tC, CI.sig.^2, CII.sig.^2, CI.echi, CII.echi, etaI, etaII, CI.c, CII.c, ...
     CI.s, CI.v(:, 3), CII.v(:, 3)]);

fid = fopen(fullfile(figDir, 'ch03_numbers.tex'), 'w');
fprintf(fid, '%% Generated by code/ch03_run.m -- do not edit by hand\n');
fprintf(fid, '%% MATLAB %s, %s\n', version('-release'), datestr(now, 'yyyy-mm-dd HH:MM'));
mac(fid, 'MatlabRelease', version('-release'));
mac(fid, 'AErrRef',   sci(norm(A1.q(end, :) - A1.qh(end, :))));
mac(fid, 'AErrPaper', sci(norm(A2.q(end, :) - A2.qh(end, :))));
mac(fid, 'AChiMaxRef',   sprintf('%.3f', max(abs(A1.q(:, 3) - A1.qh(:, 3)))));
mac(fid, 'AChiMaxPaper', sprintf('%.3f', max(abs(A2.q(:, 3) - A2.qh(:, 3)))));
mac(fid, 'ASetRef',   tsetStr(settling(tA, vecnorm(A1.q - A1.qh, 2, 2), 0.01*norm(A1.q(1, :) - A1.qh(1, :)))));
mac(fid, 'ASetPaper', tsetStr(settling(tA, vecnorm(A2.q - A2.qh, 2, 2), 0.01*norm(A2.q(1, :) - A2.qh(1, :)))));
mac(fid, 'SigMin', sprintf('%.3f', min(sig)));
mac(fid, 'SigMax', sprintf('%.3f', max(sig)));
mac(fid, 'SigMean', sprintf('%.3f', mean(sig)));
mac(fid, 'HrefMin', sprintf('%.2f', min(ellRef)));  mac(fid, 'HrefMax', sprintf('%.2f', max(ellRef)));
mac(fid, 'HpapMin', sprintf('%.2f', min(ellPap)));  mac(fid, 'HpapMax', sprintf('%.2f', max(ellPap)));
mac(fid, 'EllMean', sprintf('%.2f', ellMean));
for j = 1:nm
    mac(fid, ['B' modes{j, 3} 'Pred'], sprintf('%.2f', rPred(j)));
    mac(fid, ['B' modes{j, 3} 'Meas'], sprintf('%.2f', rMeas(j)));
    mac(fid, ['B' modes{j, 3} 'Set'],  tsetStr(tSet(j)));
    mac(fid, ['B' modes{j, 3} 'End'],  sci(eB(end, j)));
end
names = {'A', 'B', 'C'};
for j = 1:numel(kList)
    mac(fid, ['Sw' names{j} 'k'], sprintf('%d', kList(j)));
    mac(fid, ['Sw' names{j} 'MeasP'], sprintf('%.2f', sw.measP(j)));
    mac(fid, ['Sw' names{j} 'PredP'], sprintf('%.2f', sw.predP(j)));
    mac(fid, ['Sw' names{j} 'MeasR'], sprintf('%.2f', sw.measR(j)));
    mac(fid, ['Sw' names{j} 'PredR'], sprintf('%.2f', sw.predR(j)));
    mac(fid, ['Sw' names{j} 'Ratio'], sprintf('%.1f', sw.hOverTh(j)));
end
mac(fid, 'ThdMedian', sprintf('%.2f', median(abs(thd))));
mac(fid, 'ThdMax', sprintf('%.1f', max(abs(thd))));
mac(fid, 'CkA', sprintf('%d', kA));
mac(fid, 'CvNorm', sprintf('%.2f', vn));
mac(fid, 'CsigI',  sci(CI.sig(end)^2));
mac(fid, 'CsigII', sci(CII.sig(end)^2));
mac(fid, 'CvnormSq', sci(vn^2));
mac(fid, 'CsdI',  sprintf('%.2f', sdI));
mac(fid, 'CsdII', sprintf('%.2f', sdII));
mac(fid, 'CeTwoI',  sci(abs(CI.echi(201))));
mac(fid, 'CeTwoII', sci(abs(CII.echi(201))));
mac(fid, 'CsetI',  tsetStr(settling(tC, abs(CI.echi), 0.01*abs(e0))));
mac(fid, 'CsetII', tsetStr(settling(tC, abs(CII.echi), 0.01*abs(e0))));
mac(fid, 'CvnDev', sci(max(abs(vecnorm(CI.v, 2, 2) - vn))));
mac(fid, 'CsMax', sci(max(vecnorm(CI.s, 2, 2))));
mac(fid, 'CidealDevI', sci(max(abs(CI.echi - etaI))));
mac(fid, 'CdistDev', sci(max(abs(vecnorm(CI.c - [0 0 1], 2, 2) - 1))));
mac(fid, 'CriseI', sprintf('%.2f', tC(find(CI.sig.^2 >= 0.99*vn^2, 1, 'first'))));
fclose(fid);

save(fullfile(resDir, 'ch03_results.mat'), 'A1', 'A2', 'eB', 'modes', 'rPred', 'rMeas', 'sw', 'CI', 'CII');
R = struct('A1', A1, 'A2', A2, 'eB', eB, 'CI', CI, 'CII', CII);
fprintf('Results written to %s\n', figDir);
end

%% ==========================================================================
function r = simObs(prm, u, q0, qh0, tg, opts)
f = @(t, z) pairRhs(t, z, u, prm);
[~, z] = ode45(f, tg, [q0; qh0], opts);
r.q = z(:, 1:3);  r.qh = z(:, 4:6);
end

function dz = pairRhs(t, z, u, prm)
ut = u(t);  v = ut(1:3);  w = ut(4:6);
dz = [cam_sys(z(1:3), v, w); obs_spica(z(4:6), v, w, z(1:2), prm)];
end

function dz = activeRhs(~, z, a)
Rm = reshape(z(4:12), 3, 3);
q  = z(13:15);  qh = z(16:18);  v = z(19:21);
s  = q(1:2);  x = s(1);  y = s(2);
Om = [x*v(3) - v(1), y*v(3) - v(2)];
Lw = [x*y, -(1 + x^2), y;  1 + y^2, -x*y, -x];
% angular velocity: keep the point at the image centre (estimated depth)
w  = -pinv(Lw)*(a.lc*s + Om.'*qh(3));
% linear acceleration, eq. (31) of the paper
Jv = 2*[-Om(1), -Om(2), Om(1)*x + Om(2)*y];
vn2 = v.'*v;
dv = v/vn2*a.k1*(a.kdes - 0.5*vn2) + a.k2*(eye(3) - (v*v.')/vn2)*Jv.';
W  = [0, -w(3), w(2); w(3), 0, -w(1); -w(2), w(1), 0];
dc = Rm*v;
dR = Rm*W;
dq = cam_sys(q, v, w);
dqh = obs_spica(qh, v, w, s, struct('k', a.k, 'mode', 'paper'));
dz = [dc; dR(:); dq; dqh; dv];
end

function C = activePost(Z)
n = size(Z, 1);
C.c = Z(:, 1:3);
C.s = Z(:, 13:14);
C.v = Z(:, 19:21);
C.echi = Z(:, 15) - Z(:, 18);
C.sig = zeros(n, 1);
for i = 1:n
    C.sig(i) = norm([Z(i, 13)*Z(i, 21) - Z(i, 19), Z(i, 14)*Z(i, 21) - Z(i, 20)]);
end
end

function r = msdRate(ell, stiff)
disc = ell.^2/4 - stiff;
r = ell/2;
od = disc > 0;
r(od) = ell(od)/2 - sqrt(disc(od));
end

function r = envelopeRate(t, e)
env = flipud(cummax(flipud(e)));  e0 = e(1);
sel = env <= 0.3*e0 & env >= 1e-9*e0;
r = NaN;
if nnz(sel) > 10
    cf = polyfit(t(sel), log(env(sel)), 1);  r = -cf(1);
end
end

function ts = settling(t, e, thr)
ts = NaN;
idx = find(e > thr, 1, 'last');
if isempty(idx), ts = 0; elseif idx < numel(t), ts = t(idx + 1); end
end

function s = tsetStr(x)
if isnan(x), s = '--'; else, s = sprintf('%.2f', x); end
end

function mac(fid, name, val)
fprintf(fid, '\\newcommand{\\chThree%s}{%s}\n', name, val);
end

function s = sci(x)
if x == 0, s = '0'; return, end
e = floor(log10(abs(x)));  m = x/10^e;
if abs(round(m*10)/10) >= 10, m = m/10; e = e + 1; end
s = sprintf('\\num{%.1fe%d}', m, e);
end

function writeCsv(file, header, M)
fid = fopen(file, 'w');
fprintf(fid, '%s', header{1});
fprintf(fid, ',%s', header{2:end});
fprintf(fid, '\n');
fmt = [repmat('%.10g,', 1, size(M, 2) - 1), '%.10g\n'];
fprintf(fid, fmt, M.');
fclose(fid);
end
