function animateTrajectory(t, X_cl, U_cl, cfg, failureTime)

    if nargin < 5 || isempty(failureTime)
        failureTime = inf(1, cfg.N_thrusters);   % default: never fails
    end

    tol_u          = 1e-6;
    simSpeed       = 1;      % 1 = real time, >1 faster, <1 slower
    broken_len     = 0.6;
    refArrowLength = 0.5;
    ANIM_SUB       = 4;      % animation upsampling per sim step

    % ---- meta
    ref   = cfg.xRef;
    t     = t(:)';                  % 1×Nt
    Nt    = numel(t);               % number of time steps
    N_thr = cfg.N_thrusters;        % number of thrusters

    % ---- table bounds from cfg (with wall buffer) ----
    xMin = cfg.xMin - cfg.wallBuffer;
    xMax = cfg.xMax + cfg.wallBuffer;
    yMin = cfg.yMin - cfg.wallBuffer;
    yMax = cfg.yMax + cfg.wallBuffer;

    % X_cl is 6×Nt
    Xraw = X_cl.';                  % Nt×6
    x    = Xraw(:,1);
    y    = Xraw(:,2);
    th   = Xraw(:,3);

    % U_cl is N_thr×(Nt-1) or N_thr×Nt (or Nt×N_thr); normalize to N_thr×Nt
    U = normalizeU(U_cl, N_thr);
    if size(U,2) == Nt-1
        U = [U, U(:,end)];
    elseif size(U,2) ~= Nt
        error('U_cl must have Nt or Nt-1 columns to match t.');
    end

    % ---- geometry / scales
    body_r = cfg.a;
    pos    = cfg.pos;        % N_thr×2
    beta   = cfg.beta(:);    % N_thr×1 (body-frame directions)

    umax_guess  = max(max(abs(U),[],'all'), 1e-3);
    arrow_scale = 1.5*body_r/umax_guess;

    % ---- build smoother animation grid (upsampling) ----
    if Nt > 1
        dt_sim = mean(diff(t));
    else
        dt_sim = cfg.Ts;   % fallback if single sample
    end
    dt_anim  = dt_sim / ANIM_SUB;
    t_anim   = t(1):dt_anim:t(end);
    Nt_anim  = numel(t_anim);

    % interpolate x, y, th for animation
    x_anim  = interp1(t, x,  t_anim, 'pchip');
    y_anim  = interp1(t, y,  t_anim, 'pchip');
    th_anim = interp1(t, th, t_anim, 'pchip');

    % zero-order hold for U over animation grid
    U_anim = zeros(N_thr, Nt_anim);
    for i_thr = 1:N_thr
        U_anim(i_thr,:) = interp1(t, U(i_thr,:), t_anim, 'previous', 'extrap');
    end

    % ---- figure
    fig = figure('Color',0.1*[1 1 1]); hold on; axis equal;
    set(fig,'Name','Real Time Slider Animation');
    title('Real Time Slider Animation');
    xlabel('Table X');
    ylabel('Table Y');

    % table rectangle (static)
    tableX = [xMin xMax xMax xMin xMin];
    tableY = [yMin yMin yMax yMax yMin];
    plot(tableX, tableY, '-', 'Color', 0.7*[0.7 1 1], 'LineWidth', 5);

    % trajectory preview (original samples)
    plot(x, y, ':', 'LineWidth', 0.6, 'Color', 0.6*[1 1 1]);
    trace = plot(NaN,NaN,'-','LineWidth',1.4,'Color',[0.9 0.9 0.95]);

    % axis limits: include trajectory, ref, AND full table
    [xmin_traj,xmax_traj] = bounds(x);
    [ymin_traj,ymax_traj] = bounds(y);

    xmin = min([xmin_traj, ref(1), xMin]);
    xmax = max([xmax_traj, ref(1), xMax]);
    ymin = min([ymin_traj, ref(2), yMin]);
    ymax = max([ymax_traj, ref(2), yMax]);

    pad = 1; % margin around
    axis([xmin-pad xmax+pad ymin-pad ymax+pad]);

    % ---- timer text (top-left, just above square border) ----
    dx = 0.02*(xMax - xMin + eps);
    dy = 0.03*(yMax - yMin + eps);
    timerX = xMin + dx;       % slightly in from left border
    timerY = yMax + dy;       % just above top border

    hTimer = text(timerX, timerY, 't = 0.00 s', ...
                  'Color',[1 1 1], ...
                  'FontSize',12, ...
                  'HorizontalAlignment','left', ...
                  'VerticalAlignment','bottom');

    % reference orientation arrow
    quiver(ref(1),ref(2),refArrowLength*cos(ref(3)),refArrowLength*sin(ref(3)), ...
           'AutoScale','off','Color',[1 0.3 0.3], ...
           'LineWidth',2,'MaxHeadSize',0.8,'HandleVisibility','off');

    % graphics objects
    [hBody, hPts, hQuiv] = createGraphics(body_r, N_thr);   % green thruster arrows

    headLen = 1.0*body_r;
    % body axes: x_b^+ (red) and y_b^+ (blue)
    hX = quiver(NaN,NaN,NaN,NaN,0,'LineWidth',1.5,'MaxHeadSize',2.5,'Color','r'); % x_b^+
    hY = quiver(NaN,NaN,NaN,NaN,0,'LineWidth',1.5,'MaxHeadSize',2.5,'Color','b'); % y_b^+

    % --- multi-failure red markers (one per thruster, initially hidden)
    hBroken = gobjects(N_thr,1);
    for i = 1:N_thr
        hBroken(i) = quiver(NaN,NaN,NaN,NaN,0, ...
            'LineWidth',1.8,'MaxHeadSize',2.5,'Color',[0.9 0.2 0.2], ...
            'HandleVisibility','off','Visible','off');

        set(hBroken(i), 'AutoScale', 'off', 'Clipping','off');
        uistack(hBroken(i), 'top');
    end

    legend([trace, hQuiv(1)], {'trajectory','active thruster'}, ...
           'TextColor',[0.95 0.95 0.95], 'Location','best');

    % ---- per-thruster failure times (same ordering as cfg.pos / rows of U)
    ft = prepareFailureTimes_onlyTimes(failureTime, N_thr);

    realStart = tic;       % wall-clock reference
    t0        = t_anim(1); % simulation start time

    makeVideo = true;                     % set false to disable
    videoFile = 'sliderAnimation.avi';
    if makeVideo
        v = VideoWriter(videoFile, 'Motion JPEG AVI');
        v.FrameRate = 30;                 % playback FPS (independent of real-time sim)
        open(v);
    end

    % ---- animation loop
    for k = 1:Nt_anim

        % use interpolated trajectory for drawing
        set(trace,'XData',x_anim(1:k),'YData',y_anim(1:k));
        xc  = x_anim(k);
        yc  = y_anim(k);
        thk = th_anim(k);
        R   = [cos(thk) -sin(thk); sin(thk) cos(thk)];

        % body axes
        set(hX,'XData',xc,'YData',yc, ...
               'UData', headLen*cos(thk), 'VData', headLen*sin(thk));      % x_b^+
        set(hY,'XData',xc,'YData',yc, ...
               'UData',-headLen*sin(thk), 'VData', headLen*cos(thk));      % y_b^+

        % world positions of thrusters
        thr_w = (R*pos.').'+[xc yc];   % N_thr×2

        % body patch + points
        updateCircle(hBody, [xc yc], body_r);
        set(hPts,'XData',thr_w(:,1),'YData',thr_w(:,2));

        % green thrust arrows (hide once that thruster has failed)
        uf     = U_anim(:,k);                   % N_thr×1 (same ordering as cfg)
        active = abs(uf) > tol_u;

        dir_l = [cos(beta) sin(beta)];     % N_thr×2, body-frame dirs
        dir_w = (R*dir_l.').';             % N_thr×2, world-frame dirs

        L   = arrow_scale * uf;            % signed lengths
        Ux  = dir_w(:,1).*L;
        Uy  = dir_w(:,2).*L;

        tk = t_anim(k);

        % update timer text (simulation elapsed time)
        sim_t = t_anim(k) - t_anim(1);
        set(hTimer, 'String', sprintf('t = %.2f s', sim_t));

        for i = 1:N_thr
            isFailedNow = (tk >= ft(i));
            set(hQuiv(i), ...
                'XData', thr_w(i,1), ...
                'YData', thr_w(i,2), ...
                'UData', Ux(i), ...
                'VData', Uy(i), ...
                'Visible', tern(active(i) && ~isFailedNow, 'on','off'));
        end

        % red broken markers
        Lb = broken_len * body_r;
        for i = 1:N_thr
            if tk >= ft(i)
                bdir = dir_w(i,:);
                set(hBroken(i), ...
                    'XData', thr_w(i,1), ...
                    'YData', thr_w(i,2), ...
                    'UData', Lb*bdir(1), ...
                    'VData', Lb*bdir(2), ...
                    'Visible','on');
            else
                set(hBroken(i),'Visible','off');
            end
        end

        drawnow;

        if makeVideo
            frame = getframe(fig);
            writeVideo(v, frame);
        end

        % desired wall-clock time since start (scaled by simSpeed)
        simElapsed   = (t_anim(k) - t0) / simSpeed;
        wallElapsed  = toc(realStart);
        remaining    = simElapsed - wallElapsed;

        if k == 1
            % optional initial pause *in addition* to real-time sync
            pause(2.5);
            % reset reference so main run is still real-time
            realStart = tic;
        else
            if remaining > 0
                pause(remaining);
            end
        end
    end

    if makeVideo
        close(v);
        fprintf('Video saved: %s (and MP4)\n', videoFile);
    end


end

% ----------------------- helpers (same file) -----------------------------
function U = normalizeU(D, N_thr)
    D = squeeze(D);

    % We *expect* N_thr×Nt
    if ~ismatrix(D)
        error('U_cl must be a 2-D array, got size %s.', mat2str(size(D)));
    end

    if size(D,1) == N_thr
        % Already N_thr×Nt
        U = D;
    elseif size(D,2) == N_thr
        % Time × N_thr, transpose once
        U = D.';   % N_thr×Nt
    else
        error('U_cl must have %d thruster channels, got size %s.', ...
              N_thr, mat2str(size(D)));
    end
end


function [hBody, hPts, hQuiv] = createGraphics(r, N_thr)
    th    = linspace(0,2*pi,100);
    hBody = patch(r*cos(th), r*sin(th), 0.92*[1 1 1], ...
                  'EdgeColor', 0.3*[1 1 1], 'LineWidth', 1.0);
    hPts  = plot(NaN,NaN,'k.','MarkerSize',12);
    hQuiv = gobjects(N_thr,1);
    for i = 1:N_thr
        hQuiv(i) = quiver(NaN,NaN,NaN,NaN,0, ...
                          'LineWidth',1.6, 'MaxHeadSize',2.5, ...
                          'Color',[0.2 0.8 0.2], 'Visible','off');
    end
end

function updateCircle(h, c, r)
    th = linspace(0,2*pi,100);
    set(h, 'XData', c(1) + r*cos(th), ...
           'YData', c(2) + r*sin(th));
end

function s = tern(tf, a, b)
    if tf, s = a; else, s = b; end
end

function ft = prepareFailureTimes_onlyTimes(failureTime, N_thr)
    % failureTime: scalar or 1×N_thr (seconds). Inf => never fails. <0 => fail at t=0.
    if nargin < 1 || isempty(failureTime)
        ft = inf(1, N_thr);
        return;
    end
    if isscalar(failureTime)
        ft = repmat(double(failureTime), 1, N_thr);
    else
        if numel(failureTime) ~= N_thr
            error('failureTime must be scalar or length-%d vector.', N_thr);
        end
        ft = double(failureTime(:)).';
    end
    ft(~isfinite(ft)) = inf;   % NaN -> never
    ft(ft < 0)        = 0;     % negative -> fail at t=0
end
