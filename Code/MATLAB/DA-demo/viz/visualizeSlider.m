function visualizeSlider(cfg, ud)
a = cfg.a; pos = cfg.pos; beta = cfg.beta(:);
ud = ud(:);

A = [cos(beta)'; ...
     sin(beta)'; ...
     (pos(:,1)'.*sin(beta)' - pos(:,2)'.*cos(beta)')];

y = A*ud; f = y(1:2); tau = y(3);

cla; hold on; axis equal
xlim([-1.2*a 1.2*a]); ylim([-1.2*a 1.2*a])
rectangle('Position',[-a -a 2*a 2*a])
plot(pos(:,1),pos(:,2),'ko','MarkerFaceColor',[0.8 0.8 0.8])

um = max(ud); if um<=0, um=1; end
Lmax = 0.5*a;
for i = 1:length(ud)
    if ud(i) > 0
        d = [cos(beta(i)) sin(beta(i))];
        L = Lmax*(ud(i)/um);
        quiver(pos(i,1), pos(i,2), L*d(1), L*d(2), 0, 'MaxHeadSize',0.8, 'Color',[0.3 0.8 0.3])
    end
end

sf = (0.6*a)/max(norm(f), eps);
quiver(0,0, sf*f(1), sf*f(2), 0, 'LineWidth',2, 'Color','r', 'MaxHeadSize',1.0)

R = 0.7*a;
tau_scale = a*sum(abs(ud));
ang = min(2*pi, abs(tau)/max(tau_scale,eps)*2*pi);
if ang > 0
    sgn = sign(tau); if sgn==0, sgn=1; end
    t = linspace(0, sgn*ang, 100);
    x = R*cos(t); y = R*sin(t);
    plot(x, y, 'm-', 'LineWidth',2)
    te = sgn*ang;
    tv = R*[-sin(te); cos(te)]; tv = tv / max(norm(tv),eps);
    quiver(R*cos(te), R*sin(te), 0.1*a*tv(1), 0.1*a*tv(2), 0, 'Color','m', 'MaxHeadSize',1.2)
end

hold off

fprintf('Resultant vector [Fx, Fy, Tau] = [%.3f, %.3f, %.3f]\n', f(1), f(2), tau);

end
