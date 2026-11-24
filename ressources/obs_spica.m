function dqh = obs_spica(qh, v, omega, y, prm)
%OBS_SPICA  Active-estimation observer of Spica & Robuffo Giordano (CDC 2013)
%           for the inverse depth of a point feature.
%
%   dqh = obs_spica(qh, v, omega, y, prm)
%
%   Observer (paper, eq. (2)) with x_m = s = [x; y], x_u = chi:
%     d sh/dt   = f_m(s, omega) + Omega' chih + H xi
%     d chih/dt = f_u(s, chih, v, omega) + k * Omega * xi ,   xi = y - sh
%   with Omega = [x vz - vx, y vz - vy] (1-by-2), k = alpha*beta (gain on the
%   depth row) and the damping matrix H built from the SVD of Omega:
%     V = [Omega' , J Omega']/sigma  (J = 90 deg rotation),  sigma = ||Omega||
%     H = V diag(D1, D2) V' ,  D1 = 2 sqrt(k) sigma  (critical damping).
%
%   qh   estimate [xh; yh; chih]
%   y    measured image point [x; y]
%   prm  structure
%        .k      depth gain alpha*beta               (default 10)
%        .mode   'paper'     : Omega, f_m, f_u evaluated at the measured s,
%                              H with the orthonormal V above (default)
%                'reference' : reproduces the reference script: Omega, f_m,
%                              f_u evaluated at the estimate sh and
%                              V = [Omega', J Omega']/sigma^2, i.e.
%                              H = (2 sqrt(k)/sigma) I
%                'refH'      : measured s, but the H of 'reference'
%                'const'     : measured s, H = prm.ell * I
%        .D2     second damping value for 'paper' (default D1)
%        .ell    constant damping for 'const'
%        .sfloor floor on sigma (default 1e-8)

if nargin < 5, prm = struct(); end
if ~isfield(prm, 'k'),      prm.k = 10;          end
if ~isfield(prm, 'mode'),   prm.mode = 'paper';  end
if ~isfield(prm, 'sfloor'), prm.sfloor = 1e-8;   end

sh = qh(1:2);  chih = qh(3);
if strcmp(prm.mode, 'reference'), s = sh; else, s = y(:); end
x = s(1);  ys = s(2);
vx = v(1); vy = v(2); vz = v(3);
wx = omega(1); wy = omega(2);

Om    = [x*vz - vx, ys*vz - vy];
sigma = max(norm(Om), prm.sfloor);
Lw    = [x*ys, -(1 + x^2), ys;  1 + ys^2, -x*ys, -x];
fm    = Lw*omega(:);
fu    = vz*chih^2 + (ys*wx - x*wy)*chih;

switch prm.mode
    case 'paper'
        D1 = 2*sqrt(prm.k)*sigma;
        if isfield(prm, 'D2'), D2 = prm.D2; else, D2 = D1; end
        V  = [Om(1), -Om(2); Om(2), Om(1)]/sigma;          % orthonormal
        H  = V*diag([D1, D2])*V.';
    case {'reference', 'refH'}
        D1 = 2*sqrt(prm.k)*sigma;
        V  = [Om(1), -Om(2); Om(2), Om(1)]/sigma^2;        % as in the script
        H  = V*D1*V.';                                      % = (2 sqrt(k)/sigma) I
    case 'const'
        H  = prm.ell*eye(2);
    otherwise
        error('obs_spica:mode', 'Unknown mode %s.', prm.mode);
end

xi   = y(:) - sh;
dsh  = fm + Om.'*chih + H*xi;
dchi = fu + prm.k*Om*xi;
dqh  = [dsh; dchi];
end
