% a script to compare implicit and explicit schemes for doing FDTD
% simulations with a convolutional PML. 
% A little history:
% Began with a 1D code, have expanded to a 2D TM (Hz, Ex, Ey, modeled)
% Branched from Rob Marshall's CPML script, which was an explicit,
% matrix-free implementation of 2D FDTD for a goofy lens.
% 
% Implicit PML can be realized in a very simple manner. The only
% equations to modify are the Psi updates.
% in the Explicit FDTD method, they have the general form:
% Psi+ = beta*Psi + alpha*Nabla*E
% For the implicit formula, we just write out that:
% dPsi/dt = (beta-1)*Psi + alpha*Nabla*E
% which can then be solved however you want, Trapesoidal or higher
% order RK methods

% in this script, we compare explicit and implicit methods, we use
% the trapesoidal rule for the implicit scheme.


epso = 8.854e-12;
muo = 4*pi*1e-7;

% setup basic parameters
vp = 1/sqrt(epso*muo);
f = 15e6;   % wavelength = 20 m
dx = 1;  % 1 m
dy = 1;
CFL = dx/(vp*sqrt(2));
% Explicit solver blows up if dt > CFL.
% Implicit solver remains stable if dt > CFL.
dt = 2*CFL; 

% grid size for the "normal grid" 
% half-grid will have nx+1 units.
nx = 100;

% number of time steps to run 
nts = 100;

% crE is the curl that acts on E
% maps full grid to half grid
crE = spdiags(repmat([-1 1]/dx,nx+1,1),-1:0, nx+1,nx);


%crH is the curl that acts on H
% maps half grid to full grid
crH = spdiags(repmat([-1 1]/dx,nx,1), 0:1, nx,nx+1); % or d2

hz2ex = kron(crH,speye(nx+1));
hz2ey = kron(speye(nx+1), crH);

ex2hz = kron(crE,speye(nx+1));
ey2hz = kron(speye(nx+1),crE);

n = nx*(nx+1); % full vector dimension of Ex, Ey
N = (nx+1)^2;  % full vector dimension of Hz

% +++++++ Specify PML parameters +++++++++++

sxmax = -(4+1)*log(1e-6)/2/377/(10*dx);

sx = zeros(nx,1); sxm = zeros(nx+1,1);
sy = zeros(nx,1); sym = zeros(nx+1,1);


for mm = 1:10
    sy(mm) = sxmax*((11-mm)/10)^4;
    sym(mm) = sxmax*((11-mm+0.5)/10)^4;
    sy(nx+1-mm) = sxmax*((11-mm)/10)^4;
    sym(nx+2-mm) = sxmax*((11-mm+0.5)/10)^4;
    sx(mm) = sxmax*((11-mm)/10)^4;
    sxm(mm) = sxmax*((11-mm+0.5)/10)^4;
    sx(nx+1-mm) = sxmax*((11-mm)/10)^4;
    sxm(nx+2-mm) = sxmax*((11-mm+0.5)/10)^4;
end

%% a and b values for CPML ------------------------------

aex = exp(-(sx*dt/epso)) - 1; bex = aex + 1;
aex = repmat(aex,1,nx+1); bex = repmat(bex,1,nx+1);
aey = exp(-(sy*dt/epso)) - 1; bey = aey + 1;
aey = repmat(aey,1,nx+1)'; bey = repmat(bey,1,nx+1)';
ahx = exp(-(sxm*dt/epso)) - 1; bhx = ahx + 1;
ahx = repmat(ahx,1,nx+1); bhx = repmat(bhx,1,nx+1);
ahy = exp(-(sym*dt/epso)) - 1; bhy = ahy + 1;
ahy = repmat(ahy,1,nx+1)'; bhy = repmat(bhy,1,nx+1)';

% build the matrix for Maxwell dynamics
normalA = dt*[sparse(n,n) sparse(n,n) (1/epso)*hz2ex;
        sparse(n,n) sparse(n,n) -(1/epso)*hz2ey;
        (1/muo)*ex2hz (-1/muo)*ey2hz sparse(N,N)];

% build the matrix to add in the PML components to the standard
% field components
% Pex Pey Phy Phx 
PsiAdd  = [sparse(n,n) (dt/epso)*speye(n) sparse(n,N) sparse(n,N);
           (-dt/epso)*speye(n) sparse(n,n) sparse(n,N) sparse(n,N);
           sparse(N,n) sparse(N,n) (dt/muo)*speye(N) (-dt/muo)* ...
           speye(N)];

% updates from the fields to the PSI components
PsiFields = [sparse(n,n) sparse(n,n) spdiags(aex(:),0,n,n)*hz2ey;
             sparse(n,n) sparse(n,n) spdiags(aey(:),0,n,n)*hz2ex;
             spdiags(ahy(:),0,N,N)*ex2hz sparse(N,n) sparse(N,N);
             sparse(N,n) spdiags(ahx(:),0,N,N)*ey2hz sparse(N,N)];

% self updates for the Psi components
PsiSelf = [spdiags(bex(:)-1,0,n,n) sparse(n,n) sparse(n,N) ...
           sparse(n,N);
           sparse(n,n) spdiags(bey(:)-1,0,n,n) sparse(n,N) ...
           sparse(n,N);
           sparse(N,n) sparse(N,n) spdiags(bhy(:)-1,0,N,N) sparse(N, ...
                                                  N);
           sparse(N,n) sparse(N,n) sparse(N,N) spdiags(bhx(:)-1,0, ...
                                                  N,N)];
% smush it into one big matrix
A = [normalA PsiAdd;
     PsiFields PsiSelf];
           
% create trapezoidal rule matrices             
L = speye(n+n+N+n+n+N+N) - (1/2)*A;
R = speye(n+n+N+n+n+N+N) + (1/2)*A;

% factor and cache the factorization once
factorTime = tic;
fsl = factorizedLU(L);
disp(['Matrix factor time for ' num2str(size(L)) ' matrix = ' num2str(toc(factorTime))])

%putting a point source right in the middle
sourceIndx = n+n+4999;

% start with u = zero at time 1
u = zeros(n+n+N+n+n+N+N,1);

% create holders for the explicit scheme
Phy = zeros(N,1);
Phx = zeros(N,1);
Pex = zeros(n,1);
Pey = zeros(n,1);

ex = zeros(n,1);
ey = zeros(n,1);
hz = zeros(N,1);

% hold on to all history if you want
superU = zeros(n+n+N+n+n+N+N,nts);
fullLoop = tic;
for itr = 1:nts
    % create source terms
    dkn = zeros(n+n+N+n+n+N+N,1);
    dkn(sourceIndx) = -sin(2*pi*f*dt*itr);
    ssk(itr) = -sin(2*pi*f*dt*itr);
    dkh = zeros(n+n+N+n+n+N+N,1);
    dkh(sourceIndx) = -sin(2*pi*f*dt*(itr-1));
    
    % solve the implicit system
    SR = (1/2)*(dkn+dkh);
    u = fsl(R*u + SR);
    
    % uncomment to save history from destruction (its all up to you)
    % superU(:,itr) = u;
    
    % uncomment the following to have it show the fields after
    % every step (wicked slow, but great for debugging)
    %    Ex = u(n+n+(1:N));
    %    figure(33);
    %    imagesc(reshape(Ex,nx+1,nx+1))
    %    colorbar
    %    caxis([-0.2 0.2])
    
    % ========== explicit stepping section ========
    % ========== uncomment to run =================
 % tic   
 %    Phy = spdiags(bhy(:),0,N,N)*Phy + spdiags(ahy(:),0,N,N)*ex2hz* ...
 %          ex;
 %    Phx = spdiags(bhx(:),0,N,N)*Phx + spdiags(ahx(:),0,N,N)*ey2hz* ...
 %          ey;
    
 %    hz = hz + (dt/muo)*ex2hz*ex - (dt/muo)*ey2hz*ey + (dt/muo)* ...
 %         (Phy-Phx);
    
 %    hz(4999) = hz(4999) - sin(2*pi*f*dt*itr);
 %    Pey = spdiags(bey(:),0,n,n)*Pey + spdiags(aey(:),0,n,n)*hz2ex*hz;
 %    Pex = spdiags(bex(:),0,n,n)*Pex + spdiags(aex(:),0,n,n)*hz2ey*hz;
    
 %    ex = ex +(dt/epso)*hz2ex*hz + (dt/epso)*Pey;
 %    ey = ey -(dt/epso)*hz2ey*hz - (dt/epso)*Pex;
 %   toc 
 
 % plot the explicit step
 %    figure(34);
 %    imagesc(reshape(hz,nx+1,nx+1))
 %    colorbar
 %    caxis([-0.2 0.2]) 
    
end
disp(['Full loop time ' num2str(toc(fullLoop))])

% parse the fields
Ex = u(1:n);
Ey = u(n+(1:n));
Hz = u(n+n+(1:N));

% plot them at the end
figure(33);
imagesc(reshape(Hz,nx+1,nx+1)')
colorbar
xlabel('x')
ylabel('y')
title('Hz')
caxis([-0.1 0.1])

figure(31)
imagesc(reshape(Ex, nx+1,nx)')
colorbar
xlabel('x')
ylabel('y')
title('Ex')

figure(32)
imagesc(reshape(Ey, nx,nx+1)')
xlabel('x')
ylabel('y')
title('Ey')
colorbar

% figure(83); 
% imagesc(sU)
% colorbar
% caxis([-0.2 0.2])
