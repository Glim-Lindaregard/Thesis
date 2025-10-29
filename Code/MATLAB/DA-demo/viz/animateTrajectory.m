function animateTrajectory(simout, ref, opts,broken)
% AI made this but it seems to work
    if nargin < 4 || isempty(broken), broken = 0; end
    if ~isfield(opts,'broken_len'), opts.broken_len = 1.0; end 

    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'tol_u'),        opts.tol_u      = 1e-6;      end
    if ~isfield(opts,'fps'),          opts.fps        = 60;        end
    if ~isfield(opts,'saveVideo'),    opts.saveVideo  = false;     end
    if ~isfield(opts,'videoName'),    opts.videoName  = 'slider_traj.mp4'; end
    if ~isfield(opts,'simSpeed'),    opts.simSpeed  = 1; end



    % --- setup / logs ---
    cfg    = config();
    logs   = simout.logsout;

    % Pose
    state_ts = logs.get('State').Values;
    t_state  = state_ts.Time;
    Xraw     = state_ts.Data;

    x = Xraw(:,1); y = Xraw(:,2); th = Xraw(:,3);

    t = seconds(t_state - t_state(1)); t = t(:);

    % Thrusters
    u_ts = logs.get('u').Values;
    tu   = seconds(u_ts.Time - u_ts.Time(1)); tu = tu(:);
    U    = normalizeU(u_ts.Data);      % 8×N_u

    % Resample controls to pose timeline (step/zero-order hold)
    if ~isequal(size(U,2), numel(t)) || any(t ~= tu)
        U = interp1(tu, U.', t, 'previous','extrap').';  % 8×N
    end

    % --- geometry / scales ---
    if isfield(cfg,'a'), body_r = cfg.a; else, body_r = 0.2; end
    pos  = cfg.pos;   beta = cfg.beta;
    assert(isequal(size(pos),[8,2]),   'cfg.pos must be 8x2.');
    assert(numel(beta)==8,             'cfg.beta must have 8 entries.');
    umax_guess = max(max(abs(U),[],'all'), 1e-3);
    if ~isfield(opts,'arrow_scale')
        opts.arrow_scale = 0.75*body_r/umax_guess;
    end

    % --- figure ---
    fig = figure('Color',[0.18 0.18 0.18]); hold on; axis equal;
    set(fig,'Name','Slider Animation');
    plot(x, y, ':', 'LineWidth', 0.6, 'Color', 0.6*[1 1 1]);   % full path (faint)
    trace = plot(NaN,NaN,'-','LineWidth',1.4,'Color',[0.9 0.9 0.95]);

    % Axis box with margin incl. ref
    [xmin,xmax] = bounds(x); [ymin,ymax] = bounds(y);
    if ~isempty(ref)
        xmin = min(xmin, ref(1)); xmax = max(xmax, ref(1));
        ymin = min(ymin, ref(2)); ymax = max(ymax, ref(2));
    end
    pad  = 3*body_r + 0.1*max(1, max(xmax-xmin, ymax-ymin));
    axis([xmin-pad xmax+pad ymin-pad ymax+pad]);

    % Reference arrow (once)
    if ~isempty(ref)
        quiver(ref(1),ref(2),cos(ref(3)),sin(ref(3)), ...
            'AutoScale','off','Color',[0.9 0.3 0.3], ...
            'LineWidth',1.5,'MaxHeadSize',0.8,'HandleVisibility','off');
    end

    % Body + thrusters (reusable handles)
    [hBody, hPts, hQuiv] = createGraphics(body_r);
    headLen  = 1.0*body_r;
    hHeading = quiver(NaN,NaN,NaN,NaN,0,'LineWidth',1.5,'MaxHeadSize',2.5,'Color',[0.3 0.6 1.0]);

    hBroken = quiver(NaN,NaN,NaN,NaN,0, ...
    'LineWidth',1.8,'MaxHeadSize',2.5,'Color',[0.9 0.2 0.2], ...
    'HandleVisibility','off','Visible','off');


    legend([trace, hQuiv(1)], {'trajectory','active thruster'}, 'TextColor',[0.95 0.95 0.95], 'Location','best');

    % Video?
    if opts.saveVideo
        vw = VideoWriter(opts.videoName,'MPEG-4');
        vw.FrameRate = max(1, min(120, round(opts.fps)));
        open(vw);
    end

    % --- animation ---
    N = numel(t);
    dt_target = 1/max(1,opts.fps);

    for k = 1:N
        % Update trace
        set(trace,'XData',x(1:k),'YData',y(1:k));

        % Pose / rotation
        xc=x(k); yc=y(k); thk=th(k);
        R = [cos(thk) -sin(thk); sin(thk) cos(thk)];

        % Heading vector
        set(hHeading,'XData',xc,'YData',yc, ...
                     'UData',headLen*cos(thk),'VData',headLen*sin(thk));

        % Thruster world positions
        thr_w = (R*pos.').'+[xc yc];

        % Body + points
        updateCircle(hBody, [xc yc], body_r);
        set(hPts,'XData',thr_w(:,1),'YData',thr_w(:,2));

        % Thruster arrows
        uf      = U(:,k);
        active  = abs(uf) > opts.tol_u;
        dir_l   = [cos(beta) sin(beta)];           % local dirs
        dir_w   = (R*dir_l.').';                   % world dirs
        L       = opts.arrow_scale * uf;           % signed length
        Ux      = dir_w(:,1).*L;  Uy = dir_w(:,2).*L;

        if broken > 0 && broken <= 8
            % constant length (scaled by body radius)
            Lb  = opts.broken_len * body_r;
            bx  = thr_w(broken,1);  by = thr_w(broken,2);
            bUx = dir_w(broken,1) * Lb;
            bUy = dir_w(broken,2) * Lb;
        
            set(hBroken,'XData',bx,'YData',by,'UData',bUx,'VData',bUy,'Visible','on');
        else
            set(hBroken,'Visible','off');
        end


        for i = 1:8
            set(hQuiv(i), 'XData',thr_w(i,1), 'YData',thr_w(i,2), ...
                          'UData',Ux(i),      'VData',Uy(i), ...
                          'Visible', tern(active(i),'on','off'));
        end

        drawnow limitrate;

        if opts.saveVideo, writeVideo(vw, getframe(fig)); end

        % pacing
        if k < N
            dt = seconds(t(min(k+1,N)) - t(k));
            pause(opts.simSpeed*max(0, min(dt_target, dt)));    % keep it smooth but not sluggish
        end
    end

    if opts.saveVideo
        close(vw);
        fprintf('Saved video: %s\n', opts.videoName);
    end
end

% ----------------------- helpers (same file) -----------------------------

function U = normalizeU(D)
% Normalize arbitrary Simulink logging shapes to 8×N numeric
    D = squeeze(D);
    if      ismatrix(D) && size(D,2)==8        % N×8
        U = D.';                                % -> 8×N
    elseif  ismatrix(D) && size(D,1)==8        % 8×N
        U = D;
    elseif  ndims(D)==3 && size(D,3)==8        % N×1×8, etc.
        U = permute(D,[3 1 2]);                % 8×N×1/?
        U = reshape(U, 8, []);                  % 8×N
    else
        error('Logged ''u'' must be 8-wide. Got size %s after squeeze.', mat2str(size(D)));
    end
end

function [hBody, hPts, hQuiv] = createGraphics(r)
    th = linspace(0,2*pi,100);
    hBody = patch(r*cos(th), r*sin(th), 0.92*[1 1 1], ...
                  'EdgeColor', 0.3*[1 1 1], 'LineWidth', 1.0);
    hPts  = plot(NaN,NaN,'k.','MarkerSize',12);
    hQuiv = gobjects(8,1);

    for i = 1:8
        hQuiv(i) = quiver(NaN,NaN,NaN,NaN,0, ...
                          'LineWidth',1.6, 'MaxHeadSize',2.5, ...
                          'Color',[0.2 0.8 0.2], 'Visible','off');
    end


end

function updateCircle(h, c, r)
    th = linspace(0,2*pi,100);
    set(h,'XData', c(1) + r*cos(th), 'YData', c(2) + r*sin(th));
end

function s = tern(tf, a, b), if tf, s=a; else, s=b; end
end
