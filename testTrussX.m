%% create a truss

[x,y] = meshgrid(1:2);
x = [x(:), y(:)];

N = length(x);

myBody = rigidBody(2,0.1,eye(2),[0 0],0.1,[0,0.1]);

bodies(1:N) = myBody;

[a,b] = meshgrid(1:12,13:25);
edges = [1 2 1 1
         2 3 1 1
         3 4 1 1
         4 1 1 1
         1 3 1 1
         2 4 1 1];

k = 1;

mytruss = trussx(x,edges,k,bodies);

%% plot
plot2d(mytruss,[],1)

%% assemble
mytruss = mytruss.assemble();

%% static solve test
[u,f] = mytruss.staticsSolve([1 0   0   0; ...
                              2 nan nan 0; ...
                              3 nan 0   0; ...
                              4 nan nan 0], ...
                             [4 1   0   0]);

%% harmonic solve test
[u,f] = mytruss.harmonicSolve([1 0   0   0; ...
                              2 nan nan 0; ...
                              3 nan 0   0; ...
                              4 nan nan 0], ...
                             [4 1   0   0], 10);

%% ev solve test
[lambda,u,f] = mytruss.evSolve([], 12, 0);

%% triangular truss
% dimensions
n_rows = 30;
n_cols = 30;

% side length and height
dx = 1;
dy = sqrt(3)/2*dx;

% bodies
myBody = rigidBody(2,0.3,eye(2),[0 0],0.2,[0.4,0]);
myBody.N = 6;
myBody.x = myBody.R*[cos(2*pi/6*(0:myBody.N-1))' sin(2*pi/6*(0:myBody.N-1))'];

% build truss
tt = trussx.tritruss(n_rows, n_cols, dx, dy, myBody, 1);

% constraints
constraints = [(1:n_rows)'*n_cols (0:n_rows-1)'*n_cols+1 ones(n_rows,1)];
tt = tt.assemble();
tt = tt.constrain(constraints);

%% ev solve
[lambda,u,f] = tt.evSolve([], 5, 0.01);

%% plot
plot2d(tt,real(50*u(:,:,2)),0)
