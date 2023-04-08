VidMatrix

X = data;
dt = 1/vidObj.FrameRate;
t = linspace(0,vidObj.duration,677);
X1 = X(:,1:end-1);
X2 = X(:,2:end);

[U, Sigma, V] = svd(X1,'econ');
S = U'*X2*V*diag(1./diag(Sigma));

A1 = S; % Shape: 675x675 double

[eigvec,eigval] = eig(S);
mu = diag(eigval);
omega = log(mu)/dt;
[mineig, argmineig] = min(abs(omega));
mineigg = omega(argmineig);

Phi = U*eigvec(:,argmineig);
y0 = Phi\X1(:,1); % pseudoinverse to get initial conditions
u_modes = zeros(length(y0),length(t));
for iter = 1:length(t)
   u_modes(:,iter) = y0.*exp(mineigg*t(iter)); 
end
u_dmd = Phi*u_modes;

A2 = mineig; % scalar
A3 = u_dmd(:,339); % 100368x1 double

u_dmd = u_dmd(:,1:size(u_dmd,2)-1);
dataprocess = data-u_dmd;

A4 = abs(u_dmd); % 100368x1 double
A5 = abs(dataprocess); % 100368x1 double
MakeVid