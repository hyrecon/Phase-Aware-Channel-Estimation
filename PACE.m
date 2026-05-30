function [h_hat, phi1_hat, theta_hat, diag_out] = phase_aware_lifted_CE(yp, xp, n_idx, sigma_w2, mode, prior)
% PHASE_AWARE_LIFTED_CE  Algorithm 1: Phase-Aware Lifted Channel Estimation
%
%   Inputs:
%     yp       : (|Pk| x 1) complex, debiased received pilot samples  y_tilde_p[n]
%     xp       : (|Pk| x 1) complex, known pilot symbols  x_p[n]
%     n_idx    : (|Pk| x 1) real, pilot sample indices n (USE CENTERED indices
%                n = -(N-1)/2 ... (N-1)/2 so that s1 = 0; see Eq.(47))
%     sigma_w2 : scalar, noise variance  sigma_w^2
%     mode     : 'LS' or 'LMMSE'
%     prior    : (LMMSE only) struct with fields
%                  mu    = [mu_theta0; mu_theta1]   (complex 2x1)
%                  Lam   = [sigma_theta0_2; sigma_theta1_2] (real 2x1, diag of Lambda_theta)
%
%   Outputs:
%     h_hat    : estimated effective channel  h_hat  (= theta0)
%     phi1_hat : residual phase-slope estimate  phi1_hat  [rad/sample], Eq.(29)
%     theta_hat: [theta0_hat; theta1_hat]
%     diag_out : struct (Delta, cond number, etc.) for diagnostics

    yp = yp(:); xp = xp(:); n = n_idx(:);

    % Sufficient statistics
    w  = abs(xp).^2;                 % |x_p[n]|^2
    s0 = sum(w);                     % pilot energy
    s1 = sum(n .* w);                % time-centering of energy
    s2 = sum(n.^2 .* w);             % temporal spread
    t0 = sum(conj(xp) .* yp);        % matched filter   (complex)
    t1 = sum(n .* conj(xp) .* yp);   % index-weighted MF (complex)

    switch upper(mode)
        case 'LS'                                  % lines 7-11
            Delta = s0*s2 - s1^2;                  % Eq.(26)
            if Delta <= 0
                error('Gram matrix not positive definite: Delta = %g <= 0.', Delta);
            end
            theta0 = (s2*t0 - s1*t1) / Delta;      % Eq.(27)/(28)
            theta1 = (s0*t1 - s1*t0) / Delta;      % Eq.(27)
            h_hat  = theta0;
            diag_out.Delta = Delta;

        case 'LMMSE'                               % lines 12-19
            if nargin < 6 || ~isfield(prior,'mu') || ~isfield(prior,'Lam')
                error('LMMSE mode requires prior.mu and prior.Lam.');
            end
            mu  = prior.mu(:);  Lam = prior.Lam(:);
            % Eq.(36)-(37) coefficients
            b00 = 1/Lam(1) + s0/sigma_w2;
            b11 = 1/Lam(2) + s2/sigma_w2;
            b01 = s1/sigma_w2;
            a0  = t0/sigma_w2 + mu(1)/Lam(1);
            a1  = t1/sigma_w2 + mu(2)/Lam(2);
            DeltaB = b00*b11 - b01^2;              % Eq.(36)
            if DeltaB <= 0
                error('LMMSE normal matrix singular: DeltaB = %g.', DeltaB);
            end
            theta0 = (b11*a0 - b01*a1) / DeltaB;   % Eq.(37)/(38)
            theta1 = (b00*a1 - b01*a0) / DeltaB;   % Eq.(37)
            h_hat  = theta0;
            diag_out.DeltaB = DeltaB;

        otherwise
            error('mode must be ''LS'' or ''LMMSE''.');
    end

    % Residual phase-slope diagnostic
    phi1_hat  = imag(theta1 / theta0);
    theta_hat = [theta0; theta1];
end
