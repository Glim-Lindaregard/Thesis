function visualize_PWM(freq, pwm_resolution, thrust, nBlocks)
% visualize_PWM  Faithful visualization of ftc_node.py create_pwm()
%
% visualize_PWM(freq, pwm_resolution, thrust, max_force, nBlocks)
%
% freq           : PWM frequency [Hz]
% pwm_resolution : slots per PWM period (integer)
% thrust         : desired thrust (scalar)
% max_force      : maximum thrust (maps to 100% duty)
% nBlocks        : number of PWM periods to show (default = 3)
%
% Output plotted as 0/1 valve command, with alternating shaded periods.
    
    max_force = 0.7;
    if nargin < 4 || isempty(nBlocks), nBlocks = 3; end

    res = pwm_resolution;

    % ----- EXACT create_pwm() logic -----
    duty = thrust / max_force;
    duty = max(0.0, min(1.0, duty));      % clamp [0,1]
    number_of_pulses = floor(duty * res); % Python int() for positive -> floor

    pwm_block = [ones(1, number_of_pulses), ...
                 zeros(1, res - number_of_pulses)];
    % -----------------------------------

    % ----- Timing -----
    T_pwm = 1 / freq;     % one period [s]
    dt    = T_pwm / res;  % one slot / pulse width [s]

    signal = repmat(pwm_block, 1, nBlocks);
    t = (0:numel(signal)-1) * dt;

    % ----- Plot -----
    figure; clf; hold on; box on;

    yl = [-0.2, 1.2];

    % Alternating shaded background per PWM period
    for k = 0:nBlocks-1
        x0 = k * T_pwm;
        x1 = (k+1) * T_pwm;
        if mod(k,2) == 0
            patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], [0 0 0], ...
                  'FaceAlpha', 0.05, 'EdgeColor', 'none');
        end
    end

    % PWM waveform
    stairs(t, signal, 'LineWidth', 2);

    % Period boundaries
    for k = 0:nBlocks
        xline(k * T_pwm, '--', 'LineWidth', 1);
    end

    ylim(yl);
    xlim([0, nBlocks*T_pwm]);
    yticks([0 1]);
    yticklabels({'OFF','ON'});
    xlabel('Time [s]');
    ylabel('Valve command (0/1)');

    % Helpful title info
    title(sprintf(['PWM: f=%.4g Hz, res=%d | slot width dt=%.6g s (%.3f ms)\n' ...
                   'duty=%.4f, pulses=%d/%d'], ...
        freq, res, dt, 1e3*dt, duty, number_of_pulses, res));

    grid on;
end
