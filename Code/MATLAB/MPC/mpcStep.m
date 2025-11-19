function u = mpcStep(xk, xref, MPC)
% xk   : current state (nx x 1)
% xref : reference state (nx x 1)
% MPC  : struct with A,B,Q,R,Np,Nc,u_min,u_max

    A  = MPC.Ad;
    B  = MPC.Bd;
    Q  = MPC.Q;
    R  = MPC.R;
    Np = MPC.Np;
    Nc = MPC.Nc;

    [nx,nu] = size(B);

    % --- Prediction matrices Sx, Su ---
    Apow = repmat(eye(nx), 1, 1, Np);
    for i = 2:Np
        Apow(:,:,i) = Apow(:,:,i-1)*A;
    end

    Sx = zeros(nx*Np, nx);
    Su = zeros(nx*Np, nu*Nc);

    for i = 1:Np
        Sx_i = Apow(:,:,i);
        Sx((i-1)*nx+1:i*nx,:) = Sx_i;

        for j = 1:Nc
            if j <= i
                Aij = Apow(:,:,i-j+1);
                row = (i-1)*nx+1:i*nx;
                col = (j-1)*nu+1:j*nu;
                Su(row,col) = Aij*B;
            end
        end
    end

    % --- Cost matrices ---
    Qbar = kron(eye(Np), Q);
    Rbar = kron(eye(Nc), R);

    Xref = repmat(xref, Np, 1);
    Sx_xk = Sx*xk;

    H = 2*(Su'*Qbar*Su + Rbar);
    f = 2*Su'*Qbar*(Sx_xk - Xref);

    % --- Input bounds ---
    umin = MPC.u_min(:);
    umax = MPC.u_max(:);
    lb = repmat(umin, Nc, 1);
    ub = repmat(umax, Nc, 1);

    % --- Solve QP ---
    if isfield(MPC,'qp_opts')
        opts = MPC.qp_opts;
    else
        opts = optimoptions('quadprog','Display','off');
    end

    U = quadprog(H, f, [], [], [], [], lb, ub, [], opts);

    if isempty(U)
        u = zeros(nu,1);
    else
        u = U(1:nu);
    end
end
