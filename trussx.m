% Current version by Hussein Nassar (nassarh@missouri.edu)
% Earlier versions by HN and Phillip Brucks
% you are free to use and modify for research and education purposes with
% proper citation and attribution.
% for commercial use, bugs and questions, contact HN.


classdef trussx
    % trussPlus class for "trusses" where nodes are rigid bodies with a
    % finite size
    
    properties
        D       % dimension 2 or 3
        N       % number of bodies
        E       % number of edges
        bodies  % N x 1 or 1 x 1 array of rigidBody objects
        x       % N x D array of body positions
        edges   % E x 4 array of connectivity by indices
        k       % E x 1 array of spring constants
        K       % 3(D-1)N x 3(D-1)N stiffness matrix
        dir     % E x D array of unit directors
        M       % 3(D-1)N x 3(D-1)N mass matrix
        cstr    % list of constraints
    end
    
    methods
%         function x = port3d(obj,p,u)
%             % returns the port position for port p(2) of object p(1)
%             % displaced through u in 3D
%             x = obj.x(p(1)) + u(p(1),1:obj.D) + cross(u(p(1),obj.D+1:end),obj.bodies(p(1)).x(p(2)));
%         end
%         
%         function x = port2d(obj,p,u)
%             % returns the port position for port p(2) of object p(1)
%             % displaced through u in 2D
%             x = obj.x(p(1)) + u(p(1),1:obj.D) + u(p(1),obj.D+1)*circshift(obj.bodies(p(1)).x(p(2)),1).*[-1,1];
%         end
        
        function obj = trussx(x,edges,k,bodies)
            % constructor
            % in: x, edges, k
            %     body = rigidBody object
            % out: trussx object
            
            obj.N = size(x,1);
            obj.D = size(x,2);
            obj.E = size(edges,1);
            obj.bodies = bodies;
            obj.x = x;
            obj.edges = edges;
            
            if size(k) ~= obj.E
                obj.k = ones(obj.E,1)*k(1);
            else
                obj.k = k;
            end
            
            obj.cstr = [];
        end      
        
        function obj = assemble(obj)
            if obj.D == 2
                obj = assemble2d(obj);
            else
                obj = assemble3d(obj);
            end
        end
        
        function obj = assemble2d(obj)
            % assembles the stiffness/mass matrix for D=2
            
            % stats
            N = obj.N;
            D = 2;
            E = obj.E;
            
            % stiffness matrix
            % initialize sparse stiffness matrix K
            PQ = zeros((2*3*(D-1))^2*E,3);
                        
            % initialize edge directors
            dir = zeros(E,D);
            
            % assemble stiffness matrix
            for i = 1:E
                % edge i = [p q]
                p = obj.edges(i,[1,3]);
                q = obj.edges(i,[2,4]);
                
                % ports
                a = obj.bodies(p(1)).x(p(2),:);
                b = obj.bodies(q(1)).x(q(2),:);
                p = p(1);
                q = q(1);
                
                % director
                d = obj.x(q,:) + b - obj.x(p,:) - a;
                d = d/norm(d);
                dir(i,:) = d;
                
                % cross product matrices
                A = [-a(2),a(1)]*d';
                B = [-b(2),b(1)]*d';
                
                % element matrix
                dd = d'*d;
                Ad = A*d;
                Bd = B*d;
                
                Ke = [ dd     Ad'     -dd     -Bd'
                       Ad     A'*A    -Ad     -A'*B
                      -dd    -Ad'      dd      Bd'
                      -Bd    -B'*A     Bd      B'*B];
                
                % node p has 2+ displacements [ux, uy, ...]
                % corresponding rows/columns in K are [p, p+N, ...]
                [pq, qp] = meshgrid([p+(0:3*(D-1)-1)*N q+(0:3*(D-1)-1)*N]);
                PQ(1+(2*3*(D-1))^2*(i-1):(2*3*(D-1))^2*i,:) = [pq(:),qp(:),obj.k(i)*Ke(:)];
            end
            
            obj.dir = dir;
            obj.K = sparse(PQ(:,1),PQ(:,2),PQ(:,3),3*(D-1)*N,3*(D-1)*N);

            % mass matrix            
            % initialize sparse mass matrix M
            IJ = zeros((3*(D-1))^2*N,3);
            
            % assemble mass matrix
            for i = 1:N        
                m  = obj.bodies(i).M*eye(D);
                j  = obj.bodies(i).J;
                cm = obj.bodies(i).cm;
                cm = m*[-cm(2); cm(1)];
                
                mj = [m    cm
                      cm'  j];
                
                [ij,ji] = meshgrid(i+(0:3*(D-1)-1)*N);
                IJ(1+(3*(D-1))^2*(i-1):(3*(D-1))^2*i,:) = [ij(:),ji(:),mj(:)];
            end
        
            obj.M = sparse(IJ(:,1),IJ(:,2),IJ(:,3),3*(D-1)*N,3*(D-1)*N,length(IJ));
        end
        
        function obj = assemble3d(obj)
            % assembles the stiffness matrix for D=3
            
            % stats
            N = obj.N; %#ok<*PROP>
            D = 3;
            E = obj.E;
            
            % initialize sparse stiffness matrix K
            pqk = zeros((2*3*(D-1))^2*E,3);
            
            % initialize sparse mass matrix M
            
            % initialize edge directors
            dir = zeros(E,D);
            
            % assemble stiffness matrix
            for i = 1:E
                % edge i = [p q]
                p = obj.edges(i,[1,3]);
                q = obj.edges(i,[2,4]);
                
                % ports
                a = obj.bodies(p(1)).x(p(2),:);
                b = obj.bodies(q(1)).x(q(2),:);
                p = p(1);
                q = q(1);
                
                % director
                d = obj.x(q,:)+b - obj.x(p,:)-a;
                d = d/norm(d);
                dir(i,:) = d;
                
                % cross product matrices
                A = cross(a,d);
                B = cross(b,d);
                
                % element matrix
                dd = d'*d;
                Ad = A'*d;
                Bd = B'*d;
                
                Ke = [ dd     Ad'     -dd     -Bd'
                       Ad     A'*A    -Ad     -A'*B
                      -dd    -Ad'      dd      Bd'
                      -Bd    -B'*A     Bd      B'*B];
                
                % node p has 2+ displacements [ux, uy, ...]
                % corresponding rows/columns in K are [p, p+N, ...]
                [pq, qp] = meshgrid([p+(0:3*(D-1)-1)*N q+(0:3*(D-1)-1)*N]);
                pqk(1+(2*3*(D-1))^2*(i-1):(2*3*(D-1))^2*i,:) = [pq(:),qp(:),obj.k(i)*Ke(:)];
            end
            
            obj.dir = dir;
            obj.K = sparse(pqk(:,1),pqk(:,2),pqk(:,3),3*(D-1)*N,3*(D-1)*N);

            % mass matrix            
            % initialize sparse mass matrix M
            IJ = zeros((3*(D-1))^2*N,3);
            
            % assemble mass matrix
            for i = 1:N        
                m  = obj.bodies(i).M*eye(D);
                j  = obj.bodies(i).J*eye(D);
                cm = (obj.bodies(i).cm)';
                cm = m*cross(eye(D),[cm,cm,cm]); % ??
                
                mj = [m   cm
                      cm'  j];
                
                [ij,ji] = meshgrid(i+(0:3*(D-1)-1)*N);
                IJ(1+(3*(D-1))^2*(i-1):(3*(D-1))^2*i,:) = [ij(:),ji(:),mj(:)];
            end
        
            obj.M = sparse(IJ(:,1),IJ(:,2),IJ(:,3),3*(D-1)*N,3*(D-1)*N,length(IJ));
        end
        
        function [u,f, freeDofs] = prepare(obj,u0,f0)
            % prepares the user defined
            % boundary conditions u0, and
            % applied "reactions" f0
            % for the different solvers
            
            % just in case ...
            if isempty(f0)
                f0 = [1,zeros(1,3*(obj.D-1))];
            end
            
            if isempty(u0)
                u0 = [1,nan*zeros(1,3*(obj.D-1))];
            end
            
            % initialize forces            
            f0nb = size(f0,1)*3*(obj.D-1);
            f0dofs = reshape(f0(:,1)+(0:3*(obj.D-1)-1)*obj.N,f0nb,1);
            f0loads = reshape(f0(:,2:end),f0nb,1);
            % f = sparse(f0dofs,ones(f0nb,1),f0loads,3*(obj.D-1)*obj.N,1,3*(obj.D-1)*obj.N);
            f = zeros(3*(obj.D-1)*obj.N,1);
            f(f0dofs) = f0loads;
            
            % initialize dispalcements
            u0nb = size(u0,1)*3*(obj.D-1);
            u0dofs = reshape(u0(:,1)+(0:3*(obj.D-1)-1)*obj.N,u0nb,1);
            u0loads = reshape(u0(:,2:end),u0nb,1);
            
            % excepts where not assigned
            u0dofs(isnan(u0loads)) = [];
            u0loads(isnan(u0loads)) = [];
            % u0nb = length(u0loads);
            % u = sparse(u0dofs,ones(u0nb,1),u0loads,3*(obj.D-1)*obj.N,1,3*(obj.D-1)*obj.N);
            u = zeros(3*(obj.D-1)*obj.N,1);
            u(u0dofs) = u0loads;
            
            % find the free degrees of freedom
            freeDofs = setdiff(1:3*(obj.D-1)*obj.N,u0dofs);
            
            if ~isempty(obj.cstr)
                cstrDofs = obj.cstr(:,1)+(0:3*(obj.D-1))*obj.N;
                freeDofs = setdiff(freeDofs,cstrDofs(:));
            end
        end
            
        function [u,f] = staticsSolve(obj,u0,f0)
            % Statics solver
            % in:   truss obj
            %       f0 = _ x (3(D-1)+1) array of nodal loads
            %       u0 = _ x (3(D-1)+1) array of given nodal displacements
            % out:  u  = N x 3(D-1) array of nodal displacements
            %       f  = N x 3(D-1) array of reactions
            
            % prepare for the right-hand side = f0 - Ku0
            [u,f,freeDofs] = obj.prepare(u0,f0);
            
            % the right hand side is
            rhs = f - obj.K*u;
            
            % solve for free degrees of freedom
            u(freeDofs) = obj.K(freeDofs,freeDofs)\rhs(freeDofs);
            
            % find reactions
            f = obj.K*u;
            
            % reshape into x-like arrays
            u = reshape(u,obj.N,3*(obj.D-1));
            f = reshape(f,obj.N,3*(obj.D-1));
            
            % adjust for constraints
            if ~isempty(obj.cstr)
                u = follow(obj, u);
            end
        end
        
        function [u,f] = harmonicSolve(obj,u0,f0,omega)
            % Harmonic solver
            % in:   truss obj
            %       f0 = _ x (3(D-1)+1) array of nodal loads
            %       u0 = _ x (3(D-1)+1) array of given nodal displacements
            %       w  = angular frequency of excitation
            % out:  u  = N x 3(D-1) array of nodal displacements
            %       f  = N x 3(D-1) array of reactions
            
             % prepare for the right-hand side = f0 - Ku0
            [u,f,freeDofs] = obj.prepare(u0,f0);
            
            % the right hand side is
            rhs = f - obj.K*u + omega^2*obj.M*u;
            
            % solve for free degrees of freedom
            u(freeDofs) = (obj.K(freeDofs,freeDofs)-omega^2*obj.M(freeDofs,freeDofs))\rhs(freeDofs);
            
            % find reactions
            f = obj.K*u-omega^2*obj.M*u;
            
            % reshape into x-like arrays
            u = reshape(u,obj.N,3*(obj.D-1));
            f = reshape(f,obj.N,3*(obj.D-1));
            
            % adjust for constraints
            if ~isempty(obj.cstr)
                u = follow(obj, u);
            end
        end
        
        function [lambda,u,f] = evSolve(obj,u0,n,SIGMA)
            % Eigenvalues/vectors solver
            % in:   truss obj
            %       u0 = _ x (D+1) array of given nodal displacements
            %       n  = number of (smallest) eigenvalues/vectors to be returned
            %       SIGMA = eigs option
            % out:  u  = N x 3(D-1) x n array of nodal displacements
            %       f  = N x 3(D-1) x n array of reactions
            %       lambda = n x 1 array of eigenvalues
            
             % find the free degrees of freedom
            [~,~,freeDofs] = obj.prepare(u0,[]);

            % solve the eigenvalue problem
            u = zeros(obj.N*3*(obj.D-1),n);
            % [u(freeDofs,:),lambda] = eigs(obj.K(freeDofs,freeDofs),obj.M(freeDofs,freeDofs),n,SIGMA,'IsSymmetricDefinite',true);
            [u(freeDofs,:),lambda] = eigs(obj.K(freeDofs,freeDofs),obj.M(freeDofs,freeDofs),n,SIGMA);
            lambda = diag(lambda);
            
            % find reactions
            f = zeros(obj.N*3*(obj.D-1),n);
            reactionsDofs = setdiff(1:obj.N*3*(obj.D-1),freeDofs);
            f(reactionsDofs,:) = obj.K(reactionsDofs,:)*u-lambda'.*(obj.M(reactionsDofs,:)*u);
            
            % reshape u and f
            u = reshape(u,obj.N,3*(obj.D-1),n);
            f = reshape(f,obj.N,3*(obj.D-1),n);
            
            % adjust for constraints
            if ~isempty(obj.cstr)
                u = follow(obj, u);
            end
            
            % return frequencies
            lambda = sqrt(lambda);
        end
        
        function [lambda,u,f] = evFull(obj,u0)
            % Eigenvalues/vectors solver
            % in:   truss obj
            %       u0 = _ x (D+1) array of given nodal displacements
            %       n  = number of (smallest) eigenvalues/vectors to be returned
            %       SIGMA = eigs option
            % out:  u  = N x 3(D-1) x n array of nodal displacements
            %       f  = N x 3(D-1) x n array of reactions
            %       lambda = n x 1 array of eigenvalues
            
            n = obj.N*3*(obj.D-1);
            
             % find the free degrees of freedom
            [~,~,freeDofs] = obj.prepare(u0,[]);

            % solve the eigenvalue problem
            u = zeros(obj.N*3*(obj.D-1),n);
            % [u(freeDofs,:),lambda] = eigs(obj.K(freeDofs,freeDofs),obj.M(freeDofs,freeDofs),n,SIGMA,'IsSymmetricDefinite',true);
            [u(freeDofs,:),lambda] = eig(full(obj.K(freeDofs,freeDofs)),full(obj.M(freeDofs,freeDofs)));
            lambda = diag(lambda);
            
            % find reactions
            f = zeros(obj.N*3*(obj.D-1),n);
            reactionsDofs = setdiff(1:obj.N*3*(obj.D-1),freeDofs);
            f(reactionsDofs,:) = obj.K(reactionsDofs,:)*u-lambda'.*(obj.M(reactionsDofs,:)*u);
            
            % reshape u and f
            u = reshape(u,obj.N,3*(obj.D-1),n);
            f = reshape(f,obj.N,3*(obj.D-1),n);
            
            % adjust for constraints
            if ~isempty(obj.cstr)
                u = follow(obj, u);
            end
            
            % return frequencies
            lambda = sqrt(lambda);
        end

        function [time,u,el,t,f] = transientSolve(obj,f0,u0,x0,v0,T)
            % Transient solver
            % in:   truss obj
            %       f0 = @(t) _ x (D+1) array of nodal loads
            %       u0 = @(t) _ x (D+1) array of given nodal displacements
            %       x0, v0 = initial positions and velocities
            %       T = final time
            % out:  u  = ND x steps array of nodal displacements
            %       el = E x steps array of edge elongations
            %       time = 1 x steps array of times
        end
        
        function obj = constrain(obj, constraints)
            % Returns a truss with a reduced stiffness matrix
            % that enforces the constraints
            % in: constraints = n x (2+m) where m = 1 OR (3*(D-1))^2
            %     constraints(i) = [pupil guru matrix(:)'] where pupil dofs = matrix * guru dofs 
            
            n = size(constraints,1);
            D = obj.D;
            N = obj.N;
            
            for i = 1:n
                % pupil
                p = constraints(i,1);
                
                % guru
                q = constraints(i,2);
                
                % matrix
                m = constraints(i,3:end);
                m = sparse(reshape(m,sqrt(length(m)),sqrt(length(m))));
                
                % relevant dofs    
                obj.K(:,q+(0:3*(D-1)-1)*N) = obj.K(:,q+(0:3*(D-1)-1)*N) + obj.K(:,p+(0:3*(D-1)-1)*N)*m;
                obj.K(q+(0:3*(D-1)-1)*N,:) = obj.K(q+(0:3*(D-1)-1)*N,:) + m'*obj.K(p+(0:3*(D-1)-1)*N,:);
            end
            
            % update list of constraints
            obj.cstr = [obj.cstr; constraints];
        end
        
        function u = follow(obj, u)
            % construct u so that the constrained dofs follow the lead of v
            
            n = size(obj.cstr,1);
            un = size(u,3);
            
            for j = 1:un
                for i = 1:n
                    % pupil
                    p = obj.cstr(i,1);

                    % guru
                    q = obj.cstr(i,2);

                    % matrix
                    m = obj.cstr(i,3:end);
                    m = sparse(reshape(m,sqrt(length(m)),sqrt(length(m))));

                    % relevant dofs    
                    u(p,:,j) = u(q,:,j)*m.';    
                end
            end
        end
        
        function plot3d(obj,u,annotate)
            % plots the truss in obj deformed through u
            % in current figure in 3D
            
            % deformed truss port-to-port array        
            p2p = zeros(obj.E,2*obj.D);
            
            % go port to port
            for i = 1:obj.E
                p = obj.edges(i,[1,3]);
                q = obj.edges(i,[2,4]);
                p2p(i,:) = [obj.x(p(1),:)+obj.bodies(p(1)).x(p(2),:)+u(p(1),1:obj.D)+cross(u(p(1),obj.D+1:end),obj.bodies(p(1)).x(p(2),:)), ... 
                            obj.x(q(1),:)+obj.bodies(q(1)).x(q(2),:)+u(q(1),1:obj.D)+cross(u(q(1),obj.D+1),obj.bodies(q(1)).x(q(2),:))];
            end
            
            % trace edges
            plot3([p2p(:,1) p2p(:,obj.D+1)]',...
                  [p2p(:,2) p2p(:,obj.D+2)]',...
                  'k','LineWidth',2)
            hold on
            axis equal off
            
            % scatter nodes
            nodes = obj.x + u(:,1:obj.D);
            scatter3(nodes(:,1),nodes(:,2),nodes(:,3),200,[0. 0.5 0.7],'filled')
            
            % annotate nodes
            if annotate
                text(nodes(:,1),nodes(:,2),string(1:obj.N),...
                    "FontSize",15,"HorizontalAlignment","center", ...
                    "Color",'black')
            end
            
            hold off
        end
        
        function plot2d(obj,u,annotate)
            % plots the truss in obj deformed through u
            % in current figure in 2D
            
            % default
            if isempty(u)
                u = zeros(obj.N,3*(obj.D-1));
            end
            
            % deformed truss port-to-port array        
            p2p = zeros(obj.E,2*obj.D);
            
            % go port to port
            for i = 1:obj.E
                p = obj.edges(i,[1,3]);
                q = obj.edges(i,[2,4]);
                p2p(i,:) = [obj.x(p(1),:)+obj.bodies(p(1)).x(p(2),:)+u(p(1),1:obj.D)+circshift(obj.bodies(p(1)).x(p(2),:),1).*[-1,1]*u(p(1),obj.D+1), ... 
                            obj.x(q(1),:)+obj.bodies(q(1)).x(q(2),:)+u(q(1),1:obj.D)+circshift(obj.bodies(q(1)).x(q(2),:),1).*[-1,1]*u(q(1),obj.D+1)];
            end
            
            % trace edges
            plot([p2p(:,1) p2p(:,obj.D+1)]',...
                  [p2p(:,2) p2p(:,obj.D+2)]',...
                  'k','LineWidth',2)
            hold on
            axis equal off
            
            % scatter nodes
            nodes = obj.x + u(:,1:obj.D);
            viscircles(nodes,[obj.bodies.R],'Color',[0. 0.5 0.7],'linewidth',4);
            
            % annotate nodes
            if annotate
                text(nodes(:,1),nodes(:,2),string(1:obj.N),...
                    "FontSize",15,"HorizontalAlignment","center", ...
                    "Color",'black')
            end
            
            hold off
        end
    end
    
    methods(Static)
        function tt = tritruss(n_rows, n_cols, dx, dy, myBody, k)
            % returns a triangular truss
            
            % number of nodes
            N = n_rows*n_cols;

            % nodes
            [X,Y] = meshgrid(1:n_cols,1:n_rows);

            % modify into triangular lattice
            X(1:2:n_rows,:)=X(1:2:n_rows,:)+1/2;
            X = (X'-floor(n_cols/2))*dx;
            Y = Y'*dy;
            Y = Y - min(Y,[],'all') + 2*eps;

            x = [X(:) Y(:)];

            % edges
            E = 3*(n_rows-1)*(n_cols-1)+n_cols-1;
            edges = zeros(E,4);

            counter = 0;

            for i = 0:n_rows-2
                for j = 0:n_cols-2
                    node = i*n_cols + j + 1;
                    % for every node, assign three edges
                    edges(counter+1:counter+3,:) = ...
                            [node            node+1                   1     4
                             node            node+n_cols              3-mod(i,2)     6-mod(i,2)
                             node + mod(i,2) node+n_cols + mod(i+1,2) 2+mod(i,2)     5+mod(i,2)];

                    counter = counter + 3;
                end    
            end

            % take care of "edge cases"
            for i = N-n_cols+1:N-1
                edges(counter+1,:) = [i      i+1    1   4];
                counter = counter + 1;
            end
            
            % bodies
            bodies(1:N) = myBody;
            
            % build truss
            tt = trussx(x,edges,k,bodies);
        end

        function tt = fiber(N, dX, myBody, k)
            % returns a fiber
            % connecting ports 1 and 2 in myBody
            
            % nodes
            X = (0:N-1)*dX;
            Y = zeros(N,1);
            
            x = [X(:) Y(:)];

            % edges
            E = N-1;
            edges = [(1:N-1)', (2:N)', ones(E,1), 2*ones(E,1)];
            
            % bodies
            bodies(1:N) = myBody;
            
            % build truss
            tt = trussx(x,edges,k,bodies);
        end
    end        
end