function [b_hat, est_log] = end_to_end_backscatter(y0, xp, Pk, Dk_list, Ts, W, n0, sigma_w2, mode, prior, fL, fs)
% END_TO_END_BACKSCATTER  Algorithm 2: End-to-End Backscatter Baseband Processing
%
%
%   Inputs:
%     y0       : (Ntot x 1) raw complex baseband samples y0[n]      (Eq.58)
%     xp       : (|Pk| x 1) per-frame known pilot symbols x_p[n]
%     Pk       : (Nframes x 1) cell; Pk{k} = pilot sample indices of frame k
%                (1-based indices into y0; require min index >= 2, see line 3)
%     Dk_list  : (Nframes x 1) cell; Dk_list{k} = payload sample indices of frame k
%     Ts       : sampling period [s]
%     W        : CFO estimation window length [samples]   (Table III: 240)
%     n0       : reference start index for CFO derotation  (1-based, n0 >= 2)
%     sigma_w2 : noise variance sigma_w^2
%     mode     : 'LS' or 'LMMSE'  (passed to Algorithm 1)
%     prior    : (LMMSE only) struct with fields
%                  mu  = [mu_theta0; mu_theta1]        (complex 2x1)
%                  Lam = [sigma_theta0_2; sigma_theta1_2] (real 2x1)
%                pass [] for LS mode
%     fL       : FFT-LPF cutoff [Hz]  (Table III: 200 kHz)
%     fs       : sampling rate [Hz] = 1/Ts
%
%   Outputs:
%     b_hat    : decoded bits over all payload samples (column vector)
%     est_log  : (Nframes x 1) struct array with fields h_hat, phi1_hat


    % CFO estimation via lag-1 autocorrelation
    idxW  = (n0 : n0 + W - 1).';
    u_bar = mean( y0(idxW) .* conj(y0(idxW - 1)) );
    df_hat = angle(u_bar) / (2*pi*Ts);                       % Delta f_hat

    Phi0_hat = angle(y0(n0));                              
    y = y0 .* exp(-1j * (2*pi*df_hat*(n - n0)*Ts + Phi0_hat));

    % FFT-domain low-pass filter (CE branch)
    y = fft_lpf(y, fL, fs);

    % dataset-wide DC offset removal
    allP   = vertcat(Pk{:});
    d_DC   = mean( y(allP) );                                
    y_tilde = y - d_DC;                                      

    Nframes = numel(Pk);
    b_hat   = nan(N, 1);
    est_log = repmat(struct('h_hat', [], 'phi1_hat', []), Nframes, 1);

    for k = 1:Nframes
        pk = Pk{k}(:);
        % CENTER pilot indices so that s1 = 0 (Eq.47) -> bias-free LS
        n_idx = pk - mean(pk);

        [h_hat, phi1_hat] = phase_aware_lifted_CE( ...
            y_tilde(pk), xp, n_idx, sigma_w2, mode, prior);

        est_log(k).h_hat    = h_hat;
        est_log(k).phi1_hat = phi1_hat;

        dk = Dk_list{k}(:);
        zd = y_tilde(dk) / h_hat;
        b_hat(dk) = double( real(zd) >= 0 );
    end

    b_hat = b_hat(~isnan(b_hat));
end


% Algorithm 1: Phase-Aware Lifted Channel Estimation
function [h_hat, phi1_hat, theta_hat, diag_out] = phase_aware_lifted_CE(yp, xp, n_idx, sigma_w2, mode, prior)

    yp = yp(:);  xp = xp(:);  n = n_idx(:);

    % sufficient statistics
    w  = abs(xp).^2;                 % |x_p[n]|^2
    s0 = sum(w);                     % pilot energy
    s1 = sum(n .* w);                % time-centering of energy
    s2 = sum(n.^2 .* w);             % temporal spread
    t0 = sum(conj(xp) .* yp);        % matched filter        (complex)
    t1 = sum(n .* conj(xp) .* yp);   % index-weighted MF     (complex)

    diag_out = struct();

    switch upper(mode)
        case 'LS'
            Delta = s0*s2 - s1^2;                    % Eq.26
            if Delta <= 0
                error('Gram matrix not positive definite: Delta = %g <= 0.', Delta);
            end
            theta0 = (s2*t0 - s1*t1) / Delta;        % Eq.27 / Eq.28
            theta1 = (s0*t1 - s1*t0) / Delta;        % Eq.27
            h_hat  = theta0;
            diag_out.Delta = Delta;

        case 'LMMSE'
            if isempty(prior.mu) || isempty(prior.Lam)
                error('LMMSE mode requires prior.mu and prior.Lam.');
            end
            mu  = prior.mu(:);  Lam = prior.Lam(:);
            b00 = 1/Lam(1) + s0/sigma_w2;            % Eq.36-37
            b11 = 1/Lam(2) + s2/sigma_w2;
            b01 = s1/sigma_w2;
            a0  = t0/sigma_w2 + mu(1)/Lam(1);
            a1  = t1/sigma_w2 + mu(2)/Lam(2);
            DeltaB = b00*b11 - b01^2;                % Eq.36
            if DeltaB <= 0
                error('LMMSE normal matrix singular: DeltaB = %g.', DeltaB);
            end
            theta0 = (b11*a0 - b01*a1) / DeltaB;     % Eq.37 / Eq.38
            theta1 = (b00*a1 - b01*a0) / DeltaB;     % Eq.37
            h_hat  = theta0;
            diag_out.DeltaB = DeltaB;

        otherwise
            error('mode must be ''LS'' or ''LMMSE''.');
    end

    % ---- residual phase-slope diagnostic (Eq.29) ----
    phi1_hat  = imag(theta1 / theta0);
    theta_hat = [theta0; theta1];
end


% FFT-domain low-pass filter for the CE branch
function y_lp = fft_lpf(y, fL, fs)
    y = y(:);  N = numel(y);
    Y = fft(y);
    f = (0:N-1).' * (fs/N);
    f(f > fs/2) = f(f > fs/2) - fs;     % map to [-fs/2, fs/2)
    Y(abs(f) > fL) = 0;
    y_lp = ifft(Y);
end
