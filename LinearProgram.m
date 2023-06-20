classdef LinearProgram < handle
    
% primal program:
% min c'*x: x >= 0, Ax = b
% primal vector space: col(Aperp')
% dual program:
% min Ainvb'*s - bAc: s >= 0, Aperp*s = cperp
% dual vector space: col(A')

    properties
        
        % original data
        Aoriginal double 
        boriginal double 
        coriginal double 
        N double 
        
        % current state of LP
        n double % number of variables
        m double % number of linear equalities
        A double % coefficient matrix
        b double % right-hand side
        c double % cost function
        c0 double % cost function offset
        undetermined logical % boolean vector indicating undetermined variables
        primal_zero logical % boolean vector indicating x_i = 0
        dual_zero logical % boolean vector indicating s_i = 0
        
        test_x logical % which i to test for x_i = 0
        test_s logical % which i to test for s_i = 0
        sumx double % how many columns in inv_vectors will be occupied by quantities used for testing lower bounds
        sums double % how many columns in inv_vectors will be occupied by quantities used for testing lower bounds
        
        % auxiliary quantities
        w double % scaling point
        iw double % inverse scaling point
        w2 double % scaling point squared
        iw2 double % inverse scaling point squared
        w2c double % W^2*c
        Aw2c double % A*w2c
        Aw double % A*w
        AwIAw double % (A*w)'*inv(A*W^2*A')*Aw
        nAwIAw double % n - (A*w)'*inv(A*W^2*A')*Aw
        Awb double % (A*w)'*inv(A*W^2*A')*b
        Pc double % (I - A'*inv(A*W^2*A')*A*W^2)*c
        wPc double % w'*Pc
        inv_vectors double % inv(A*W^2*A')*[A*w, b, A*w2c, w_ia_i], i only for variables to be tested
        cAWAc double % w2c'*Pc
        addon_c double % wPc^2/cAWAc
        bAWAb double % obj.b'*inv_vectors(:,2)
        addon_b double % Awb^2/bAWAb

        % bounds on primal and dual objective
        % outer bounds on interval have to intersect the tau > 0 (mu > 0)
        % sector
        transition_outer logical
        max_outer double
        min_outer double
        inner_max_positive logical
        max_inner double
        min_inner double
        
        % more simple bounds on optimal value
        lower_bound = -Inf % from unbounded
        upper_bound = +Inf % to infeasible
        
        Aperp double % complement of A
        cperp double % Aperp*c
        Ainv double % (A*A')\A
        Ainvb double % Ainv'*b
        bAc double % c'*Ainvb
        M double
        Mperp double % feasible directions in x,s,t space
        
        % global parameter
        tol = 10^(-10) % tolerance
        ga double % neighbourhood of central path
        % ga = 10^(-3); % neighbourhood of central path
        test_threshold = 0.001 % for w_i smaller we test whether x_i = 0
        
    end
    
    methods
        
        function obj = LinearProgram(n,m,A,b,c,c0)
            % constructor
            if nargin == 2
                A = [randn(m-1,n); ones(1,n)]; % feasible set bounded
                b = A*abs(randn(n,1)); % random positive vector is feasible
                c = randn(n,1);
            end
            obj.Aoriginal = A;
            obj.boriginal = b;
            obj.coriginal = c;
            obj.A = A;
            obj.b = b;
            obj.c = c;
            if nargin == 6
                obj.c0 = c0;
            else
                obj.c0 = 0;
            end
            obj.N = n;
            obj.n = n;
            obj.m = m;
            obj.ga = (n-2*sqrt(n-1))/(2*(n-1));
            obj.undetermined = true(n,1);
            obj.primal_zero = false(n,1);
            obj.dual_zero = false(n,1);
            obj.transition_outer = [];
            obj.inner_max_positive = [];
        end
      
        
        function [solved,solution] = check_special_cases(obj)
            solved = true;
            if obj.m == 0
                % no linear constraints
                if min(obj.c) >= 0
                    x = zeros(obj.n,1);
                    s = obj.c;
                    value = obj.c0;
                    small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),true(obj.n,1),false(obj.n,1));
                    fprintf("No primal linear constraints left, solution at origin\n")
                else
                    small_solution = LinearSolution(true,false);
                    fprintf("No primal linear constraints left, problem unbounded\n")
                end
            elseif max(abs(obj.b)) <= obj.tol
                Ap = null(obj.A)';
                Ag = [Ap, Ap*(obj.c - ones(obj.n,1))];
                bg = Ap*obj.c;
                cg = [zeros(obj.n,1); 1];
                fprintf("Primal linear constraints homogeneous, solving auxiliary problem\n")
                auxLP = LinearProgram(size(Ag,2),size(Ag,1),Ag,bg,cg);
                solLP = solve(auxLP);
                value = solLP.value;
                s = solLP.x;
                delete(auxLP);
                delete(solLP);
                if value <= obj.tol
                    % x = 0 optimal
                    x = zeros(obj.n,1);
                    value = obj.c0;
                    small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),true(obj.n,1),false(obj.n,1));
                else
                    % problem unbounded
                    small_solution = LinearSolution(true,false);
                end
            elseif obj.m == obj.n
                % linear constraints determine point
                x = obj.A\obj.b;
                if min(x) >= 0
                    value = obj.c'*x + obj.c0;
                    s = zeros(obj.n,1);
                    small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),false(obj.n,1),true(obj.n,1));
                    fprintf("No dual linear constraints left, primal solution feasible\n")
                else
                    small_solution = LinearSolution(false);
                    fprintf("No dual linear constraints left, primal problem infeasible\n")
                end
            elseif obj.m == 1
                % dual problem is a line search
                % s = c - y*A'
                lby = -Inf;
                uby = +Inf;
                for k = 1:obj.n
                    if obj.A(k) == 0
                        if obj.c(k) < 0
                            % dual infeasible
                            lby = +Inf;
                            uby = -Inf;
                        end
                    elseif obj.A(k) > 0
                        uby = min(uby,obj.c(k)/obj.A(k));
                    else
                        lby = max(lby,obj.c(k)/obj.A(k));
                    end
                end
                if lby > uby
                    % dual infeasible
                    small_solution = LinearSolution(true,false);
                else
                    if obj.b > 0
                        y = uby;
                    else
                        y = lby;
                    end
                    if isinf(y)
                        % dual unbounded
                        small_solution = LinearSolution(false);
                    else
                        s = obj.c - y*obj.A';
                        value = obj.b*y + obj.c0;
                        [~,minind] = min(s); % element of s which vanishes
                        pz = true(obj.n,1);
                        pz(minind) = false;
                        dz = ~pz;
                        x = zeros(obj.n,1);
                        x(minind) = obj.b/obj.A(minind);
                        small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),pz,dz);
                    end
                end
            elseif obj.n == obj.m + 1
                % primal problem is a line search
                xinh = obj.A\obj.b;
                x0 = null(obj.A);
                % x = xinh - t*x0
                % compute bounds on t
                lbt = -Inf;
                ubt = +Inf;
                for k = 1:obj.n
                    if x0(k) == 0
                        if xinh(k) < 0
                            % primal infeasible
                            lbt = +Inf;
                            ubt = -Inf;
                        end
                    elseif x0(k) > 0
                        ubt = min(ubt,xinh(k)/x0(k));
                    else
                        lbt = max(lbt,xinh(k)/x0(k));
                    end
                end
                if lbt > ubt
                    % primal infeasible
                    small_solution = LinearSolution(false);
                else
                    cost_coef = obj.c'*x0;
                    if cost_coef > 0
                        t = ubt;
                    else
                        t = lbt;
                    end
                    if isinf(t) && (cost_coef ~= 0)
                        % primal unbounded
                        small_solution = LinearSolution(true,false);
                    else
                        x = xinh - t*x0;
                        value = obj.c'*x + obj.c0;
                        [~,minind] = min(x); % element of x which vanishes
                        dz = true(obj.n,1);
                        dz(minind) = false;
                        pz = ~dz;
                        s = zeros(obj.n,1);
                        s(minind) = 1;
                        y = [obj.A', s]\obj.c;
                        s(minind) = y(end);
                        small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),pz,dz);
                    end
                end
            else
                if max(abs(obj.Pc)) <= obj.tol
                    % cost function is constant
                    % checking feasibility of primal
                    Ag = [obj.A, obj.b - obj.A*ones(obj.n,1)];
                    bg = obj.b;
                    cg = [zeros(obj.n,1); 1];
                    fprintf("Dual linear constraints homogeneous, solving auxiliary problem\n")
                    auxLP = LinearProgram(obj.n+1,obj.m,Ag,bg,cg);
                    solLP = solve(auxLP);
                    value = solLP.value;
                    x = solLP.x;
                    delete(auxLP);
                    delete(solLP);
                    if value <= obj.tol
                        % primal feasible
                        s = zeros(obj.n,1);
                        value = obj.c'*x + obj.c0;
                        small_solution = LinearSolution(true,true,x,s,value,false(obj.n,1),false(obj.n,1),true(obj.n,1));
                    else
                        % primal infeasible
                        small_solution = LinearSolution(false);
                    end
                elseif obj.n == 2
                    % approximating cone exact
                    % now this never happens because n <= 3 is catched
                    % above
                    obj.A = obj.A/obj.b;
                    obj.b = 1;
                    fprintf("Dimension of problem decreased to two\n")
                    if max(obj.A) <= 0
                        % primal infeasible
                        small_solution = LinearSolution(false);
                        fprintf("Problem infeasible\n")
                    elseif (obj.A(1) > 0) && (det([obj.A', obj.c]) > 0)
                        x = [1/obj.A(1); 0];
                        value = obj.c(1)/obj.A(1);
                        s = [0; obj.c(2) - obj.A(2)*value];
                        pz = [false; true];
                        dz = [true; false];
                        small_solution = LinearSolution(true,true,x,s,value,false(2,1),pz,dz);
                        fprintf("Problem solved\n")
                    elseif (obj.A(2) > 0) && (det([obj.A', obj.c]) < 0)
                        x = [0; 1/obj.A(2)];
                        value = obj.c(2)/obj.A(2);
                        s = [obj.c(1) - obj.A(1)*value; 0];
                        dz = [false; true];
                        pz = [true; false];
                        small_solution = LinearSolution(true,true,x,s,value,false(2,1),pz,dz);
                        fprintf("Problem solved\n")
                    else
                        % primal unbounded
                        small_solution = LinearSolution(true,false);
                        fprintf("Problem unbounded\n")
                    end
                else
                    % checking for infeasibility or unboundedness from
                    % approximation
                    % now we know that the primal and the dual are
                    % inhomogeneous
                    if (obj.AwIAw + obj.addon_c < 1) || ((obj.AwIAw <= 1) && (obj.wPc <= 0))
                        % dual half-space mu >= 0 does not intersect outer dual cone, primal unbounded
                        small_solution = LinearSolution(true,false);
                        fprintf("Problem unbounded\n")
                    elseif (obj.nAwIAw + obj.addon_b < 1) || ((obj.nAwIAw <= 1) && (obj.Awb <= 0))
                        % primal infeasible
                        small_solution = LinearSolution(false);
                        fprintf("Problem infeasible\n")
                    else
                        solved = false;
                        solution = [];
                    end
                end
            end
            if solved
                solution = restore_big_solution(obj,small_solution);
                delete(small_solution);
            end
        end
        
        function compute_quantities(obj,w)
            % computation of quantities necessary to make the next step
            % the main parameter is the scaling point w
            % this includes the quantities necessary for the check of x_i =
            % 0 or s_i = 0
            % computes auxiliary quantities w2, w2c, Aw2c, Aw,
            % inv_vectors, Pc, wPc, cAWAc, addon_c, Awb, bAWAb, addon_b, AwIAw
            % assumes that m > 0
            obj.w = w;
            obj.iw = 1./w;
            obj.test_x = (1./w < obj.test_threshold);
            obj.test_s = (w < obj.test_threshold);
            obj.w2 = w.^2;
            obj.iw2 = obj.iw.^2;
            obj.w2c = obj.w2.*obj.c;
            obj.Aw2c = obj.A*obj.w2c;
            obj.Aw = obj.A*w;
            obj.sumx = sum(obj.test_x);
            obj.sums = sum(obj.test_s);
            obj.inv_vectors = (obj.A*diag(obj.w2)*obj.A')\[obj.Aw, obj.b, obj.Aw2c, obj.A(:,obj.test_x)*diag(w(obj.test_x)), obj.A(:,obj.test_s)*diag(w(obj.test_s))];
            obj.Pc = obj.c - obj.A'*obj.inv_vectors(:,3);
            obj.wPc = w'*obj.Pc;
            obj.cAWAc = obj.w2c'*obj.Pc;
            obj.addon_c = obj.wPc^2/obj.cAWAc;
            obj.Awb = obj.Aw'*obj.inv_vectors(:,2);
            obj.bAWAb = obj.b'*obj.inv_vectors(:,2);
            obj.addon_b = obj.Awb^2/obj.bAWAb;
            obj.AwIAw = obj.Aw'*obj.inv_vectors(:,1);
            obj.nAwIAw = obj.n - obj.AwIAw;
        end
        
        function update_objective_bounds(obj)
            % strengthens primal and dual objective bounds from inner and
            % outer cone approximations if possible
            % inner cone: a*n/(a+1) = 1
            % upper bound on optimum
            d = 1 - obj.AwIAw;
            if obj.addon_b >= -d
                if (abs(d) < obj.tol) && (obj.Awb > 0)
                    upp_bnd = obj.Awb*obj.wPc*(1/2*(obj.bAWAb/obj.Awb^2 - obj.cAWAc/obj.wPc^2) - d/8*(obj.bAWAb/obj.Awb^2 + obj.cAWAc/obj.wPc^2)^2) + obj.Aw2c'*obj.inv_vectors(:,2) + obj.c0;
                elseif (d > 0) || (obj.Awb > 0)
                    upp_bnd = (sqrt((obj.wPc^2 - d*obj.cAWAc)*(obj.Awb^2 + d*obj.bAWAb)) - obj.Awb*obj.wPc)/d + obj.Aw2c'*obj.inv_vectors(:,2) + obj.c0;
                else
                    upp_bnd = +Inf;
                end
                if ~isinf(upp_bnd) && (obj.upper_bound > upp_bnd)
                    obj.upper_bound = upp_bnd;
                end
            end
            % outer cone: a*n/(a+1) = n-1
            % lower bound on optimum
            d = obj.nAwIAw - 1;
            if obj.addon_c >= d
                if (abs(d) < obj.tol) && (obj.wPc > 0)
                    low_bnd = obj.Awb*obj.wPc*(1/2*(obj.bAWAb/obj.Awb^2 - obj.cAWAc/obj.wPc^2) - d/8*(obj.bAWAb/obj.Awb^2 + obj.cAWAc/obj.wPc^2)^2) + obj.Aw2c'*obj.inv_vectors(:,2) + obj.c0;
                elseif (d < 0) || (obj.wPc > 0)
                    low_bnd = (sqrt((obj.wPc^2 - d*obj.cAWAc)*(obj.Awb^2 + d*obj.bAWAb)) - obj.Awb*obj.wPc)/d + obj.Aw2c'*obj.inv_vectors(:,2) + obj.c0;
                else
                    low_bnd = -Inf;
                end
                if ~isinf(low_bnd) && (obj.lower_bound < low_bnd)
                    obj.lower_bound = low_bnd;
                end
            end
        end
        
        function [x,s] = step(obj)
            if obj.n/2 > obj.AwIAw + obj.addon_c
                % dual totally infeasible for canonical cone
                [x,s] = dual_orth_target_point(obj);
            elseif obj.n/2 < obj.AwIAw - obj.addon_b
                % primal totally infeasible for canonical cone
                [x,s] = primal_orth_target_point(obj);
            else
                diff_crit = obj.n/2-obj.AwIAw;
                if ((diff_crit < 0) && (obj.Awb > 0)) || ((diff_crit > 0) && (obj.wPc > 0))
                    % regular solution exists for approximating cone
                    mu = obj.n*sqrt(max(0,obj.Awb^2 + diff_crit*obj.bAWAb));
                    tau = sqrt(max(0,obj.wPc^2 - diff_crit*obj.cAWAc));
                elseif diff_crit > 0
                    % primal regularly unbounded for approximation
                    mu = -1;
                    tau = 0;
                else
                    % dual regularly unbounded for approximation
                    mu = 0;
                    tau = -1;
                end
                % computing x and s candidates
                AtI = obj.A'*obj.inv_vectors;
                PiTw = obj.w - obj.w2.*AtI(:,1);
                a1wPicDen = obj.wPc/diff_crit;
                ms = obj.Pc + a1wPicDen*AtI(:,1);
                mx = (a1wPicDen*PiTw - obj.w2.*obj.Pc)/obj.n;
                a1wIbDen = obj.Awb/diff_crit;
                ts = -obj.n*(AtI(:,2) + a1wIbDen*AtI(:,1));
                tx = -a1wIbDen*PiTw + obj.w2.*AtI(:,2);
                x = mx*mu + tx*tau;
                s = ms*mu + ts*tau;
            end
        end
        
        function [x,s] = dual_orth_target_point(obj)
            % computes target point on orthogonal complement of dual feasible subspac
            x = obj.w - diag(obj.w2)*obj.A'*obj.inv_vectors(:,1) - obj.wPc*(diag(obj.w2)*obj.Pc)/obj.cAWAc;
            s = 2*(obj.n - obj.AwIAw - obj.addon_c)*obj.iw - obj.n*obj.iw2.*x;
        end
        
        function [x,s] = primal_orth_target_point(obj)
            % computes target point on orthogonal complement of primal feasible subspace
            s = obj.A'*obj.inv_vectors(:,1) - obj.Awb*(obj.A'*obj.inv_vectors(:,2))/obj.bAWAb;
            x = 2*(obj.AwIAw - obj.addon_b)*obj.w - obj.n*obj.w2.*s;
        end
        
        function [lb,S] = x0bounds(obj)
            % computes lower bounds on objective under assumption x_i = 0
            lb = zeros(1,sum(obj.test_x));
            S = zeros(obj.n,sum(obj.test_x));
            ind = find(obj.test_x);
            R = obj.A(:,obj.test_x)*diag(obj.w(obj.test_x));
            XAw = obj.inv_vectors(:,4:end)'*obj.Aw;
            for k = 1:length(ind)
                xi = obj.inv_vectors(:,k+3);
                xr = R(:,k)'*xi;
                Aw2c_ciwi = (xi'*obj.Aw2c - obj.w(ind(k))*obj.c(ind(k)));
                AwIAw1 = obj.AwIAw + (XAw(k)^2 - 2*XAw(k) + xr)/(1-xr);
                Awb1 = obj.Awb + (xi'*obj.b)*(XAw(k) - 1)/(1-xr);
                bAWAb1 = obj.bAWAb + (xi'*obj.b)^2/(1-xr);
                cAWAc1 = obj.cAWAc - Aw2c_ciwi^2/(1-xr);
                wPc1 = obj.wPc - Aw2c_ciwi*(XAw(k) - 1)/(1-xr);
                d = obj.n - 2 - AwIAw1;
                inv11 = obj.inv_vectors(:,1) + (xi'*obj.Aw - 1)/(1-xr)*xi;
                inv21 = obj.inv_vectors(:,2) + (xi'*obj.b)/(1-xr)*xi;
                inv31 = obj.inv_vectors(:,3) + Aw2c_ciwi/(1-xr)*xi;
                if d + Awb1^2/bAWAb1 < 0
                    % dual problem with s_i released is unbounded
                    lb(k) = +Inf;
                    % computing certificate in H_P^{\perp}
                    y = - inv11 + Awb1/bAWAb1*inv21;
                    S(:,k) = -obj.A*y;
                elseif (d < 0) && (Awb1 <= 0)
                    % dual problem with s_i released is unbounded
                    lb(k) = +Inf;
                    % computing certificate with mu = 0, tau = -1
                    y = -(obj.n-1)*( inv21 + Awb1/d*inv11 );
                    S(:,k) = -obj.A'*y;
                elseif wPc1^2/cAWAc1 >= d
                    if (abs(d) < obj.tol) && (wPc1 > 0) && (Awb1 > 0)
                        lb(k) = Awb1*wPc1*(1/2*(bAWAb1/Awb1^2 - cAWAc1/wPc1^2) - d/8*(bAWAb1/Awb1^2 + cAWAc1/wPc1^2)^2) + obj.Aw2c'*obj.inv_vectors(:,2) + (xi'*obj.b)*Aw2c_ciwi/(1-xr) + obj.c0;
                        if lb(k) > obj.upper_bound + obj.tol
                            y = (obj.n-1)*(inv31*sqrt(Awb1^2 + d*bAWAb1) + inv21*sqrt(wPc1^2 - d*cAWAc1) - ((bAWAb1/Awb1 - d/4*bAWAb1^2/Awb1^3)*wPc1 + (cAWAc1/wPc1 + d/4*cAWAc1^2/wPc1^3)*Awb1)/2*inv11);
                            mu = (obj.n-1)*sqrt(Awb1^2 + d*bAWAb1);
                            S(:,k) = -obj.A'*y + mu*obj.c;
                        end
                    elseif (d < 0) || (wPc1 > 0)
                        % regular solution exists
                        lb(k) = (sqrt((wPc1^2 - d*cAWAc1)*(Awb1^2 + d*bAWAb1)) - Awb1*wPc1)/d + obj.Aw2c'*obj.inv_vectors(:,2) + (xi'*obj.b)*Aw2c_ciwi/(1-xr) + obj.c0;
                        if lb(k) > obj.upper_bound + obj.tol
                            mu = (obj.n-1)*sqrt(Awb1^2 + d*bAWAb1);
                            tau = sqrt(wPc1^2 - d*cAWAc1);
                            mu_y = inv31 + (obj.Aw'*obj.inv_vectors(:,3) - obj.c'*obj.w + Aw2c_ciwi*(xi'*obj.Aw - 1)/(1-xr))/d*inv11;
                            tau_y = (obj.n-1)*(inv21 + Awb1/d*inv11);
                            y = mu_y*mu + tau_y*tau;
                            S(:,k) = -obj.A'*y + mu*obj.c;
                        end
                    else
                        % there exists a primal tau = 0, mu = -1 recession
                        % direction
                        lb(k) = -Inf;
                    end
                else
                    % for the reduced problem, H_D^{\perp} intersects the
                    % interior of the primal cone
                    % dual infeasible
                    lb(k) = -Inf;
                end
            end
        end
        
        function [ub,X] = s0bounds(obj)
            % computes upper bounds on objective under assumption s_i = 0
            % only indices i in test_s are considered
            % returns also certificates x in matrix X
            ub = zeros(1,sum(obj.test_s));
            X = zeros(obj.n,sum(obj.test_s));
            ind = find(obj.test_s);
            R = obj.A(:,obj.test_s)*diag(obj.w(obj.test_s));
            XAw = obj.inv_vectors(:,4+obj.sumx:end)'*obj.Aw;
            for k = 1:length(ind)
                xi = obj.inv_vectors(:,k+3+obj.sumx);
                xr = R(:,k)'*xi;
                Aw2c_ciwi = (xi'*obj.Aw2c - obj.w(ind(k))*obj.c(ind(k)));
                AwIAw1 = obj.AwIAw - XAw(k)^2/xr;
                Awb1 = obj.Awb - (xi'*obj.b)*XAw(k)/xr;
                bAWAb1 = obj.bAWAb - (xi'*obj.b)^2/xr;
                cAWAc1 = obj.cAWAc + Aw2c_ciwi^2/xr;
                wPc1 = obj.wPc + Aw2c_ciwi*XAw(k)/xr;
                d = 1 - AwIAw1;
                w_i = zeros(obj.n,1);
                w_i(ind(k)) = obj.w(ind(k));
                P_e_i = obj.w2.*(obj.A'*xi) - w_i;
                PiW1 = obj.w - obj.w2.*(obj.A'*obj.inv_vectors(:,1)) + XAw(k)/xr*P_e_i;
                PiW2c1 = obj.w2.*obj.c - obj.w2.*(obj.A'*obj.inv_vectors(:,3)) + Aw2c_ciwi/xr*P_e_i;
                if d > wPc1^2/cAWAc1
                    % primal problem with x_i released is unbounded
                    ub(k) = -Inf;
                    % computing certificate in H_D^{\perp}
                    X(:,k) = PiW1 - wPc1/cAWAc1*PiW2c1;
                elseif (d > 0) && (wPc1 <= 0)
                    % primal problem with x_i released is unbounded
                    ub(k) = -Inf;
                    % computing certificate with tau = 0, mu = -1
                    X(:,k) = 1/(obj.n-1)*( - wPc1/d*PiW1 + PiW2c1 );
                elseif d >= -Awb1^2/bAWAb1
                    if (abs(d) < obj.tol) && (Awb1 > 0) && (wPc1 > 0)
                        % regular solution exists
                        % compute series with respect to d to avoid
                        % numerical errors caused by division by d
                        ub(k) = Awb1*wPc1*(1/2*(bAWAb1/Awb1^2 - cAWAc1/wPc1^2) - d/8*(bAWAb1/Awb1^2 + cAWAc1/wPc1^2)^2) + obj.Aw2c'*obj.inv_vectors(:,2) - (xi'*obj.b)*Aw2c_ciwi/xr + obj.c0;
                        if ub(k) < obj.lower_bound - obj.tol
                            X(:,k) = - Awb1*PiW2c1 + wPc1*obj.w2.*(obj.A'*obj.inv_vectors(:,2)) - (xi'*obj.b)/xr*wPc1*P_e_i + (1/2*bAWAb1/Awb1 - d/8*bAWAb1^2/Awb1^3)*( wPc1*PiW1 - d*PiW2c1 ) + (- 1/2*cAWAc1/wPc1 - d/8*cAWAc1^2/wPc1^3)*(-Awb1*PiW1 + d*obj.w2.*(obj.A'*obj.inv_vectors(:,2)) - d*(xi'*obj.b)/xr*P_e_i);
                        end
                    elseif (d > 0) || (Awb1 > 0)
                        % regular solution exists
                        ub(k) = (sqrt((wPc1^2 - d*cAWAc1)*(Awb1^2 + d*bAWAb1)) - Awb1*wPc1)/d + obj.Aw2c'*obj.inv_vectors(:,2) - (xi'*obj.b)*Aw2c_ciwi/xr + obj.c0;
                        if ub(k) < obj.lower_bound - obj.tol
                            % computing certificate
                            mu = (obj.n-1)*sqrt(Awb1^2 + d*bAWAb1);
                            tau = sqrt(wPc1^2 - d*cAWAc1);
                            mu_x = 1/(obj.n-1)*( wPc1/d*PiW1 - PiW2c1 );
                            tau_x = -Awb1/d*PiW1 + obj.w2.*(obj.A'*obj.inv_vectors(:,2)) - (xi'*obj.b)/xr*P_e_i;
                            X(:,k) = mu*mu_x + tau*tau_x;
                        end
                    else
                        % there exists a dual mu = 0, tau = -1 recession
                        % direction
                        ub(k) = +Inf;
                    end
                else
                    ub(k) = +Inf;
                end
            end
        end
        
        function solution = solve(obj)
          
            tic,
            countIter = 1;
            x = ones(obj.n,1);
            s = ones(obj.n,1);
            fprintf("NumIter | primal feas gap | dual feas gap | step length | lower bound | upper bound | gap | NumVars | time\n")
            while true
                w = scaling_point(obj,x,s);
                if obj.m > 0
                    compute_quantities(obj,w);
                end
                [solved,solution] = check_special_cases(obj);
                if ~solved
                    update_objective_bounds(obj);
                    if obj.upper_bound == +Inf
                        obj.test_x = false(obj.n,1);
                    end
                    if obj.lower_bound == -Inf
                        obj.test_s = false(obj.n,1);
                    end
                    if any(obj.test_x)
                        % leave those test_x which entail x_i > 0
                        [lb,S] = x0bounds(obj);
                        ind = find(obj.test_x);
                        obj.test_x(ind(lb <= obj.upper_bound + obj.tol)) = false;
                        S(:,lb <= obj.upper_bound + obj.tol) = [];
                    end
                    if any(obj.test_s)
                        % leave those test_s which entail x_i = 0
                        [ub,X] = s0bounds(obj);
                        ind = find(obj.test_s);
                        obj.test_s(ind(ub >= obj.lower_bound - obj.tol)) = false;
                        X(:,ub >= obj.lower_bound - obj.tol) = [];
                    end
                    [xnew,snew] = step(obj);
                    [x,s,step_length] = make_step(obj,x,s,xnew,snew,obj.ga);
                    % adjust x and s by adding conic combination of
                    % certificates
                    if any(obj.test_s)
                        u = -X(obj.test_s,:)\x(obj.test_s);
                        x = x + X*u;
                    end
                    if any(obj.test_x)
                        u = -S(obj.test_x,:)\s(obj.test_x);
                        s = s + S*u;
                    end
                    % eliminating entries with x_i = 0 and s_i = 0
                    while any(obj.test_x)
                        % reduce primal variables which are free
                        ind = find(obj.test_x,1);
                        obj.test_x(ind) = [];
                        obj.test_s(ind) = [];
                        x(ind) = [];
                        s(ind) = [];
                        reduce_system_dual_zero(obj,ind);
                    end
                    while any(obj.test_s)
                        % reduce primal variables which are zero
                        ind = find(obj.test_s,1);
                        obj.test_s(ind) = [];
                        obj.test_x(ind) = [];
                        x(ind) = [];
                        s(ind) = [];
                        reduce_system_primal_zero(obj,ind);
                    end
                    fprintf("%5d | %5e | %5e | %5e | %5e | %5e | %5e | %5d | %5e \n",countIter,obj.AwIAw-obj.addon_b,obj.nAwIAw-obj.addon_c,step_length,obj.lower_bound,obj.upper_bound,obj.upper_bound-obj.lower_bound,obj.n,toc)
                    if obj.n > 0
                        % normalize x and s
                        coef = sqrt(obj.n/(x'*s));
                        x = x*coef;
                        s = s*coef;
                    else
                        small_solution = LinearSolution(true,true,zeros(0,1),zeros(0,1),0,false(1,0),false(1,0),true(1,0));
                        solution = restore_big_solution(obj,small_solution);
                        return;
                    end
                    countIter = countIter + 1;
                else
                    return;
                end
            end
        end
        
        function [x,s,step_length] = make_step(obj,x,s,xnew,snew,ga)
            %update the primal-dual pair (x,s) to (xnew,snew)
            xnew = xnew*obj.n/(xnew'*s);
            snew = snew*obj.n/(x'*snew);
            % finding where (x,s) leave the central path
            % neighbourhood
            coef0 = x.*s - ga*ones(obj.n,1);
            coef1 = (xnew-x).*s + x.*(snew-s);
            coef2 = (xnew-x).*(snew-s) + ga*(1 - xnew'*snew/obj.n);
            if any((abs(coef0) < 10^(-13)) & (coef1 <= 0))
                fprintf("No progress possible\n")
            end
            disc = coef1.^2 - 4*coef0.*coef2;
            disc_pos = (disc > 0);
            coef1 = coef1(disc_pos);
            coef2 = coef2(disc_pos);
            disc = sqrt(disc(disc_pos));
            all_roots = [(-coef1-disc)./(2*coef2); (-coef1+disc)./(2*coef2)];
            all_roots(all_roots < 10^(-14)) = [];
            step_length = min([1; all_roots]);
            x = x + step_length*(xnew - x);
            s = s + step_length*(snew - s);
        end
        
        function solution = restore_big_solution(obj,small_solution)
            % small solution is the solution of the problem on the current

            if ~small_solution.feas
                solution = LinearSolution(false);
            elseif ~small_solution.bounded
                solution = LinearSolution(true,false);
            else
                % feasible and bounded
                x = zeros(obj.N,1);
                s = zeros(obj.N,1);
                x(obj.undetermined) = small_solution.x;
                s(obj.undetermined) = small_solution.s;
                x(obj.dual_zero) = obj.Aoriginal(:,obj.dual_zero)\(obj.boriginal - obj.Aoriginal(:,obj.undetermined)*small_solution.x);
                y = obj.Aoriginal(:,~obj.primal_zero)'\(obj.coriginal(~obj.primal_zero) - s(~obj.primal_zero));
                s(obj.primal_zero) = obj.coriginal(obj.primal_zero) - obj.Aoriginal(:,obj.primal_zero)'*y;
                value = obj.coriginal'*x;
                ind = find(obj.undetermined);
                if ~isempty(small_solution.primal_zero)
                    obj.undetermined(ind(small_solution.primal_zero)) = false;
                    obj.primal_zero(ind(small_solution.primal_zero)) = true;
                end
                if ~isempty(small_solution.dual_zero)
                    obj.undetermined(ind(small_solution.dual_zero)) = false;
                    obj.dual_zero(ind(small_solution.dual_zero)) = true;
                end
                solution = LinearSolution(true,true,x,s,value,obj.undetermined,obj.primal_zero,obj.dual_zero);
            end
        end
        
        function reduce_system_dual_zero(obj,i)
            [~,maxind] = max(abs(obj.A(:,i)));
            av = obj.A(maxind,:);
            bv = obj.b(maxind)/av(i);
            av = av/av(i);
            obj.b = obj.b - obj.A(:,i)*bv;
            obj.A = obj.A - obj.A(:,i)*av;
            obj.A(maxind,:) = [];
            obj.b(maxind) = [];
            c0 = obj.c(i)*bv;
            obj.c = obj.c - obj.c(i)*av';
            obj.A(:,i) = [];
            obj.c(i) = [];
            obj.c0 = obj.c0 + c0;
            ind = find(obj.undetermined);
            obj.undetermined(ind(i)) = false;
            obj.dual_zero(ind(i)) = true;
            obj.n = obj.n - 1;
            obj.m = obj.m - 1;
            obj.ga = (obj.n-2*sqrt(obj.n-1))/(2*(obj.n-1));
        end

        function reduce_system_primal_zero(obj,i)
            % consider the system Ax = b and cost c'*x
            % suppose the variable x_i is zero
            obj.A(:,i) = [];
            obj.c(i) = [];
            obj.n = obj.n - 1;
            obj.ga = (obj.n-2*sqrt(obj.n-1))/(2*(obj.n-1));
            ind = find(obj.undetermined);
            obj.undetermined(ind(i)) = false;
            obj.primal_zero(ind(i)) = true;
        end
        
        function [mx,ms,tx,ts] = mu_tau_coefs(obj,w,a)
            % computes solution of system
            % Ax = t*b
            % s + A'y = m*c
            % w*s = ( (a + 1){\bf 1} - nI )(x/w)
            % solution is of the form
            % x = mx*m + tx*t
            % s = ms*m + ts*t
            w2 = w.^2;
            w2c = w2.*obj.c;
            Aw2c = obj.A*w2c;
            Aw = obj.A*w;
            inv_vectors = (obj.A*diag(w2)*obj.A')\[Aw, obj.b, Aw2c];
            AtI = obj.A'*inv_vectors;
            denom = a*obj.n - (a+1)*(w'*AtI(:,1));
            PiTw = w.*(1 - w.*AtI(:,1));
            Pic = obj.c - AtI(:,3);
            a1Den = (a+1)/denom;
            a1wPicDen = a1Den*(w'*Pic);
            ms = Pic + a1wPicDen*AtI(:,1);
            mx = (a1wPicDen*PiTw - w2c + w2.*AtI(:,3))/obj.n;
            AwI2 = w'*AtI(:,2);
            a1wIbDen = a1Den*AwI2;
            ts = -obj.n*(AtI(:,2) + a1wIbDen*AtI(:,1));
            tx = -a1wIbDen*PiTw + w2.*AtI(:,2);
        end
        
        function [x,s,nesterovTime] = nesterovSolve(obj)
            tic;
            fprintf("Solving linear program\n")
            fprintf("Finding feasible pair\n")
            fprintf("NumIter | Primal infeas | Dual infeas \n")
            x = ones(obj.n,1);
            s = ones(obj.n,1);
            ndP = 1;
            ndD = 1;
            countIter = 0;
            while max(ndP,ndD) > 10^(-12)
                [wp,wd] = scaling_point(obj,x,s);
                [feas_primal,xnew] = min_Q_on_affine_space(obj,wp,obj.A,obj.b,obj.c);
                [feas_dual,snew] = min_Q_on_affine_space(obj,wd,obj.Aperp,obj.cperp,obj.Ainvb);
                dx = xnew - wp;
                ds = snew - wd;
                coefb = dx.*wd+ds.*wp;
                dxds = 2*dx.*ds;
                discr = coefb.^2 - 2*dxds*0.999;
                coefb = coefb(discr >= 0);
                dxds = dxds(discr >= 0);
                discr = discr(discr >= 0);
                if ~isempty(discr)
                    mu = [(-coefb+sqrt(discr))./dxds; (-coefb-sqrt(discr))./dxds];
                    mu = mu((mu > 0) & (mu < 1));
                    if isempty(mu)
                        x = xnew;
                        s = snew;
                    else
                        mu = min(mu);
                        x = wp + mu*dx;
                        s = wd + mu*ds;
                    end
                else
                    x = xnew;
                    s = snew;
                end
                dP = primal_feas_gap(obj,x,1);
                dD = dual_feas_gap(obj,s,1);
                ndP = norm(dP);
                ndD = norm(dD);
                countIter = countIter + 1;
                fprintf("%5d | %5e | %5e \n",countIter,ndP,ndD)
            end
            fprintf("Primal-dual feasible pair achieved\n")
            fprintf("Starting main phase\n")
            fprintf("NumIter | Dual obj | Primal obj | Gap | NumConstr \n")
            t = 1;
            opt = false;
            countIter = 1;
            while ~opt
                [x,s,t] = nesterovStep(obj,x,s,t,countIter);
                countIter = countIter + 1;
            end
            nesterovTime = toc;
            fprintf("Time spent: %5e sec\n",nesterovTime);
        end
        
        function dP = primal_feas_gap(obj,x,t)
            dP = obj.A*x - t*obj.b;
        end
        
        function dD = dual_feas_gap(obj,s,t)
            dD = obj.Aperp*s - t*obj.cperp;
        end
        
        function xmin = estimate_variables(obj,x,C,c)
            % bounds x_i from below on outer approximation
            % {z = x + C*y}, c'*z <= c'*x, z'*(1/x*1/x' -
            % diag(1/x^2))*z >= 0
            % xmin is a vector of lower bounds on the components of x
            % if c is not given, we go without the corresponding constraint
            xi = 1./x;
            Axi = C'*xi;
            if nargin == 4
                Ac = C'*c;
                value = c'*x;
                C0 = C*null(Ac');
                xmin0 = estimate_variables(obj,x,C0);
            end
            H = C'*(diag(xi.^2))*C;
            iH = inv(H);
            iHAxi = iH*Axi;
            critEPH = Axi'*iHAxi;
            xmin = zeros(obj.n,1);
            if critEPH < 1 - 10^(-9)
                % elliptic case
                Hess = iH + iHAxi*iHAxi'/(1-critEPH);
                projektor = -C*Hess*C';
                lasqNum = -(obj.n-1)*(obj.n-critEPH)/(1-critEPH);
                lasqDen = diag(projektor);
                la = sqrt(lasqNum./lasqDen);
                Z = projektor*(diag(la) - (obj.n-1)*xi*ones(1,obj.n)) + x*ones(1,obj.n);
%                 lasq_c_Den = -Ac'*Hess*Ac;
%                 la_c = sqrt(lasqNum/lasq_c_Den);
%                 zcmin = projektor*(la_c*c - (obj.n-1)*xi) + x;
%                 zcmax = projektor*(-la_c*c - (obj.n-1)*xi) + x;
                xmin = diag(Z);
                if nargin == 4
                    min_not_feas = (c'*Z > value);
                    xmin(min_not_feas) = xmin0(min_not_feas);
                end
            elseif critEPH < 1 + 10^(-9)
                % parabolic case
                % hatw = iHAxi
                eChatw = C*iHAxi;
                if nargin == 4
                    cChatw = Ac'*iHAxi;
                end
                for k = 1:obj.n
                    if eChatw > 0
                        w = [1/(obj.n-1)*C'*(xi*xi'-diag(xi.^2))*C; C(k,:)/eChatw(k) + Axi']\[C(k,:)'/eChatw(k)-Axi; -obj.n];
                        z = C*w + x;
                        xmin(k) = z(k);
                        if (nargin == 4) && (c'*z > value)
                            xmin(k) = xmin0(k);
                        end
                    else
                        xmin(k) = -Inf;
                        if (nargin == 4) && (cChatw > 0)
                            xmin(k) = xmin0(k);
                        end                                
                    end
                end
            else
                % hyperbolic case
                Hess = -iH + iHAxi*iHAxi'/(critEPH-1);
                projektor = C*Hess*C';
                % w0 = -(obj.n-1)*Hess*Axi;
                % hatw = iHAxi
                dProj = diag(projektor);
                for k = 1:obj.n
                    if (dProj(k) > 0) && (C(k,:)*iHAxi > 0)
                        lasqNum = (obj.n-1)*(obj.n-critEPH)/(critEPH-1);
                        lasqDen = dProj(k);
                        la = sqrt(lasqNum/lasqDen);
                        ek = zeros(obj.n,1);
                        ek(k) = 1;
                        z = projektor*(la*ek - (obj.n-1)*xi) + x;
                        xmin(k) = z(k);
                        if (nargin == 4) && (c'*z > value)
                            xmin(k) = xmin0(k);
                        end
                    else
                        xmin(k) = -Inf;
                        if nargin == 4
                            G = [c'*projektor*c, projektor(k,:)*c; projektor(k,:)*c, dProj(k)];
                            if (G(1,1) > 0) && (Ac'*iHAxi > 0)
                                xmin(k) = xmin0(k);
                            elseif (G(2,2) < 0) && (G(1,1) < 0) && (G(1,2) >= sqrt(G(1,1)*G(2,2))) && (G(1,2)*C(k,:)*iHAxi >= G(2,2)*Ac'*iHAxi)
                                xmin(k) = xmin0(k);
                            end
                        end
                    end
                end
            end
        end
        
        function [x,s,t] = nesterovStep(obj,x,s,t,countIter)
            dP = obj.A*x - t*obj.b;
            dD = obj.Aperp*s - t*obj.cperp;
            dcost = obj.c'*x + obj.Ainvb'*s - t*obj.bAc - obj.n;
            dstep = obj.M\[-dP; -dD; -dcost]; % correction to maintain feasibility
            x = x + dstep(1:obj.n);
            s = s + dstep(obj.n+1:2*obj.n);
            t = t + dstep(2*obj.n+1);
            xmin = estimate_variables(obj,x/t,obj.Aperp',obj.c);
            smin = estimate_variables(obj,s/t,obj.A',obj.Ainvb);
            H = obj.Mperp'*diag([1./x.^2; 1./s.^2; 0])*obj.Mperp;
            g = -[1./x; 1./s; 0]'*obj.Mperp;
            dir_forward = -obj.Mperp*(H\g');
            d = [1./x; 1./s; 0]'*dir_forward/(2*obj.n);
            ind_neg = (dir_forward < 0);
            v = [x; s; t];
            a = min(-v(ind_neg)./dir_forward(ind_neg));
            alpha_coef = a*(d*a + 1)/(a + 1);
            v = v + alpha_coef*dir_forward;
            x = v(1:obj.n);
            s = v(obj.n+1:2*obj.n);
            t = v(2*obj.n+1);
            pobj = obj.c'*x/t;
            dobj = -obj.Ainvb'*s/t+obj.bAc;
            fprintf("%d | %5e | %5e | %5e | \n",countIter,dobj,pobj,pobj-dobj)
        end
        
        function [wp,wd] = scaling_point(obj,x,s)
            wp = sqrt(x./s);
            wd = 1./wp;
        end
        
        function [x,s,y,t,h,hmax0] = homogeneousStep(obj,x,s,y,t)
            alpha = obj.n-sqrt(obj.n); % penalty weight at proximity function
            is = 1./s;
            Ax = obj.A*x;
            As = obj.A*is;
            x2 = x.^2;
            s2 = is.^2;
            cx2 = obj.c.*x2;
            cs2 = obj.c.*s2;
            Acx2 = obj.A*cx2;
            Acs2 = obj.A*cs2;
            AXA = obj.A*diag(x2)*obj.A';
            ASA = obj.A*diag(s2)*obj.A';
            Dx = Ax - t*obj.b;
            Ds = s + obj.A'*y - t*obj.c;
            D = obj.c'*x - obj.b'*y - obj.n;
            ADs = obj.A*(s2.*Ds);
            AXAinv = AXA\[Ax, Acx2, obj.b, Dx];
            ASAinv = ASA\[As, Acs2, obj.b, ADs];
            ld = [obj.c'*cx2 - Acx2'*AXAinv(:,2) + obj.b'*ASAinv(:,3), Acx2'*AXAinv(:,3) - obj.b'*ASAinv(:,2); Acs2'*ASAinv(:,3) - obj.b'*AXAinv(:,2), obj.b'*AXAinv(:,3) + obj.c'*cs2 - Acs2'*ASAinv(:,2) - alpha/t^2]\[Acx2'*AXAinv(:,1) - obj.c'*x - obj.b'*ASAinv(:,1); -Acs2'*ASAinv(:,1) + obj.b'*AXAinv(:,1) + obj.c'*is - alpha/t];
            la = ld(1);
            dt = ld(2);
            w = -la*AXAinv(:,2) - AXAinv(:,1) + dt*AXAinv(:,3);
            dy = -ASAinv(:,1) + dt*ASAinv(:,2) - la*ASAinv(:,3);
            dx = x2.*(obj.A'*w) + la*cx2 + x;
            ds = -obj.A'*dy + dt*obj.c;
            dxs = [dx; ds]./[x; s];
            hmax = -1/min(dxs);
            hmax0 = hmax;
            hmin = 0;
            while hmax-hmin > 10^(-10)
                h = (hmax+hmin)/2;
                df = -sum(dxs./(1+h*dxs))+alpha*dt/(t+h*dt); % second term is a proximity penalty
                if df < 0
                    hmin = h;
                else
                    hmax = h;
                end
            end
        end
    end  

    
end