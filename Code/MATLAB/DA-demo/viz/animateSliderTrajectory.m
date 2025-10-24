function animateSliderTrajectory(simout,ref)
    close all;
    cfg = config();
    logs = simout.logsout;
    
    % Pose (must be a 3-wide vector per sample, logged as 'state')
    state_ts = logs.get('state').Values;
    t_state  = state_ts.Time;
    X        = state_ts.Data;             % could be N×3 or 3×N depending on block
    if size(X,2)==3
        x = X(:,1); y = X(:,2); th = X(:,3);
    elseif size(X,1)==3
        x = X(1,:).'; y = X(2,:).'; th = X(3,:).';
    else
        error('Logged ''state'' must have width 3. Got size %s.', mat2str(size(X)));
    end
    xytheta = timetable(seconds(t_state - t_state(1)), x, y, th, ...
                        'VariableNames', {'x','y','theta'});
    
    % Thrusters (must be 8-wide vector per sample, logged as 'u')
    u_ts    = logs.get('u').Values;
    t_u     = u_ts.Time;
    D       = u_ts.Data;    % can be N×8, 8×N, N×1×8, etc.
    
    % --- normalize to 8×N numeric U ---
    D = squeeze(D);
    if      size(D,2)==8        % N×8
        U = D.';                % -> 8×N
    elseif  size(D,1)==8        % 8×N
        U = D;                  % already 8×N
    elseif  ndims(D)==3 && size(D,3)==8
        U = permute(D,[3 2 1]); % 8×?×?  (rare path)
        U = reshape(U,8,[]);    % 8×N
    else
        error('Logged ''u'' must be 8-wide. Got size %s after squeeze.', mat2str(size(D)));
    end
    
    % Build timetable with exactly 8 variable names
    uNames = "u" + (1:8);
    ut = array2timetable(U.', ...
        'RowTimes', seconds(t_u - t_u(1)), ...
        'VariableNames', cellstr(uNames));

    opts = struct('body_radius',0.2,'fps',120);
    
    % Animate
    animateSliderTrajectory2(xytheta, ut, cfg,opts,ref);
end


function animateSliderTrajectory2(xytheta, ulog, cfg, opts,ref)
% ANIMATESLIDERTRAJECTORY  Animate slider pose and thruster firing.
%
% Inputs
%   xytheta : pose log. Any ONE of:
%       - timetable with variables x, y, theta
%       - struct with fields x, y, theta, t (vectors)
%       - Nx3 numeric [x y theta], in which case t = (0:N-1)'
%   ulog    : thruster log. Any ONE of:
%       - timetable with variables u1..u8 (or a single 8x1 vector entry)
%       - struct with fields U (8xN) and t
%       - 8xN numeric
%   cfg     : struct with geometry (same style you've used before)
%       .a     (characteristic half-size, used for body radius if present)
%       .pos   (8x2) local thruster positions, rows = [xi yi]
%       .beta  (8x1) local thruster pointing angles (rad)
%   opts    : optional
%       .tol_u        (default 1e-6)  threshold to consider "active"
%       .body_radius  (default = cfg.a or 0.15)
%       .arrow_scale  (default = 0.25*body_radius/umax_guess)
%       .fps          (default = min(60, 1/median(dt)))
%       .saveVideo    (false)  set true to write MP4
%       .videoName    ('slider_traj.mp4')
%
% Example (after Simulink run):
%   % suppose you have xyt and ut as timetables from postProcess(simout)
%   animateSliderTrajectory(xyt, ut, cfg);

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'tol_u'),       opts.tol_u       = 1e-6; end
    if ~isfield(opts,'body_radius')
        if isfield(cfg,'a'), opts.body_radius = cfg.a;
        else,                opts.body_radius = 0.15; end
    end
    if ~isfield(opts,'saveVideo'),   opts.saveVideo   = false; end
    if ~isfield(opts,'videoName'),   opts.videoName   = 'slider_traj.mp4'; end

    % --- normalize inputs to arrays ---
    [t, X, Y, TH] = normalizePose(xytheta);
    [tu, U]       = normalizeThrusters(ulog);

    refx = ref(1); refy = ref(2); reft  = ref(3);


    % sync lengths by resampling/interp if needed
    if ~isequal(tu, t)
        % linear interp to pose timestamps
        U = interp1(tu, U.', t, 'previous', 'extrap').';  % step-like control
    end

    N = numel(t);
    assert(size(U,1)==8, 'ulog must have 8 thruster channels.');

    % --- geometry ---
    pos  = cfg.pos;   % 8x2
    beta = cfg.beta;  % 8x1
    assert(size(pos,1)==8 && size(pos,2)==2, 'cfg.pos must be 8x2.');
    assert(numel(beta)==8, 'cfg.beta must have 8 entries.');

    % rough umax guess for arrow scaling
    uabsmax = max( max(abs(U), [], 'all'), 1e-3 );
    if ~isfield(opts,'arrow_scale')
        opts.arrow_scale = 0.75*opts.body_radius/uabsmax;
    end

    % --- plotting scaffold ---
    figure('Color',[0.2,0.2,0.2]); hold on; axis equal;
    % set a stable axis that fits the whole trajectory with margin

    % --- axis limits including reference if present ---
    xmin = min(X); xmax = max(X);
    ymin = min(Y); ymax = max(Y);
    
    if exist('refx','var') && ~isempty(refx)
        xmin = min([xmin; refx(:)]);
        xmax = max([xmax; refx(:)]);
    end
    if exist('refy','var') && ~isempty(refy)
        ymin = min([ymin; refy(:)]);
        ymax = max([ymax; refy(:)]);
    end
    
    pad = 3*opts.body_radius + 0.1*max(1, max([xmax-xmin, ymax-ymin]));
    axis([xmin-pad xmax+pad ymin-pad ymax+pad]);



    % plot full trajectory (faint), and a "so-far" trace that we update
    plot(X, Y, ':', 'LineWidth', 0.5, 'Color', 0.7*[1 1 1]);
    trace = plot(NaN, NaN, '-', 'LineWidth', 1.5);

    % create body + thruster graphics (reused each frame)
    [hBody, hThrPts, hThrQuiv] = createSliderGraphics(opts.body_radius);

    % heading arrow fixed to the body (shows slider orientation theta)
    headLen = 1.0 * opts.body_radius;   % tune length if you like
    hHeading = quiver(NaN,NaN,NaN,NaN, 0, ...
    'LineWidth', 1.6, 'MaxHeadSize', 2.5, 'Color', [0 0 0.8]);  % dark blue


    % legend
    lg = legend([trace, hThrQuiv(1)], {'trajectory', 'active thruster'}, ...
                'Location','best'); %#ok<NASGU>

    % optional video
    if opts.saveVideo
        vw = VideoWriter(opts.videoName, 'MPEG-4');
        vw.FrameRate = min(60, max(1, round(1/median(max(diff(t),eps)))));
        open(vw);
    end

    % --- animation loop (update handles only) ---
    for k = 1:N
        % update trace
        set(trace, 'XData', X(1:k), 'YData', Y(1:k));
        pause(0.2);   % 0.02 s pause → 50 FPS


        quiver(refx,refy,cos(reft),sin(reft), ...
        'AutoScale','off','Color','r','LineWidth',1.5,'MaxHeadSize',...
        0.8,'HandleVisibility','off');

        % pose
        x = X(k); y = Y(k); th = TH(k);
        R = [cos(th) -sin(th); sin(th) cos(th)];


        % heading vector in world frame (unit * headLen)
        hx = headLen*cos(th);
        hy = headLen*sin(th);
        set(hHeading, 'XData', x, 'YData', y, 'UData', hx, 'VData', hy);



        % transform thruster points to world
        thr_world = (R*pos.').'+[x y];

        % update body circle
        updateCircle(hBody, [x y], opts.body_radius);

        % update thruster marker dots
        set(hThrPts, 'XData', thr_world(:,1), 'YData', thr_world(:,2));

        % update thruster arrows
        uf = U(:,k);
        active = abs(uf) > opts.tol_u;

        % Start points = thruster world positions
        X0 = thr_world(:,1);  Y0 = thr_world(:,2);

        % Direction vectors from beta rotated by body yaw
        dir_local = [cos(beta), sin(beta)];
        dir_world = (R*dir_local.').';  % 8x2

        % Length proportional to |u|
        L = opts.arrow_scale * uf;
        Ux = dir_world(:,1) .* L;
        Uy = dir_world(:,2) .* L;

        for i = 1:8
            set(hThrQuiv(i), 'XData', X0(i), 'YData', Y0(i), ...
                             'UData', Ux(i), 'VData', Uy(i), ...
                             'Visible', onOff(active(i)));
        end

        drawnow;

        if opts.saveVideo
            writeVideo(vw, getframe(gcf));
        end
    end

    if opts.saveVideo
        close(vw);
        fprintf('Saved video: %s\n', opts.videoName);
    end
end

% ---------- helpers ----------

function [t, X, Y, TH] = normalizePose(xytheta)
    if istimetable(xytheta)
        % expect variables named (case-insensitive) x, y, theta
        vn = lower(string(xytheta.Properties.VariableNames));
        X  = xytheta{:, vn=="x"};
        Y  = xytheta{:, vn=="y"};
        TH = xytheta{:, vn=="theta"};
        t  = seconds(xytheta.Properties.RowTimes - xytheta.Properties.RowTimes(1));
        t  = t(:);
    elseif isstruct(xytheta) && isfield(xytheta,'x')
        X = xytheta.x(:); Y = xytheta.y(:); TH = xytheta.theta(:);
        if isfield(xytheta,'t'), t = xytheta.t(:);
        else, t = (0:numel(X)-1)'; end
    else
        A = xytheta;
        assert(size(A,2)==3, 'xytheta numeric must be Nx3 [x y theta].');
        X = A(:,1); Y = A(:,2); TH = A(:,3); t = (0:size(A,1)-1)';
    end
end

function [t, U] = normalizeThrusters(ulog)
    if istimetable(ulog)
        % Accept either variables u1..u8 or single matrix column with 8 rows per sample
        V = ulog.Variables;
        if ismatrix(V) && size(V,2)==8
            U = V.';                       % 8xN
        else
            % try to pick vars named u1..u8
            vn = string(ulog.Properties.VariableNames);
            idx = zeros(1,8);
            for k = 1:8
                cand = strcmpi(vn, sprintf('u%d',k));
                if any(cand), idx(k) = find(cand,1); end
            end
            assert(all(idx>0), 'ulog timetable must have u1..u8 or a single 8-wide matrix variable.');
            U = ulog{:, idx}.';
        end
        t = seconds(ulog.Properties.RowTimes - ulog.Properties.RowTimes(1));
        t = t(:);
    elseif isstruct(ulog) && isfield(ulog,'U')
        U = ulog.U;  assert(size(U,1)==8, 'ulog.U must be 8xN');  t = ulog.t(:);
    else
        U = ulog;     assert(size(U,1)==8, 'ulog numeric must be 8xN'); t = (0:size(U,2)-1)';
    end
end

function [hBody, hThrPts, hThrQuiv] = createSliderGraphics(r)
    % body circle (as patch)
    th = linspace(0,2*pi,100);
    hBody = patch(r*cos(th), r*sin(th), 0.92*[1 1 1], ...
                  'EdgeColor', 0.2*[1 1 1], 'LineWidth', 1.0);

    % thruster dots (8)
    hThrPts = plot(NaN, NaN, 'k.', 'MarkerSize', 12);

    % 8 quiver arrows (created invisible, shown when active)
    hold on;
    hThrQuiv = gobjects(8,1);
    for i = 1:8
        hThrQuiv(i) = quiver(NaN,NaN,NaN,NaN, 0, ...
                             'LineWidth', 1.6, ...
                             'MaxHeadSize', 2.5, ...
                             'Color', [0.1 0.6 0.1], ...
                             'Visible','off');
    end
end

function updateCircle(h, c, r)
    th = linspace(0,2*pi,100);
    set(h, 'XData', c(1) + r*cos(th), 'YData', c(2) + r*sin(th));
end

function s = onOff(tf)
    if tf, s='on'; else, s='off'; end
end
