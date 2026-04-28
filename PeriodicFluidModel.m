classdef PeriodicFluidModel < handle
% PERIODICFLUIDMODEL  Cyclic stochastic fluid model with time-varying rates.
%
% Models a single-buffer fluid system {(X(t), J(t)) : t >= 0} with a
% periodic background CTMC having generator T(t) = T(t+1).  Computes the
% asymptotic periodic distribution and related quantities following the
% LST approach of Margolius & O'Reilly (2016) and the companion paper.
%
% USAGE
%   model = PeriodicFluidModel(T_fun, C_diag, S_plus, S_minus)
%   model = PeriodicFluidModel(T_fun, C_diag, S_plus, S_minus, Name, Value)
%
% REQUIRED INPUTS
%   T_fun   - Function handle @(t) returning an m x m generator matrix.
%             For the constant-rate case, pass @(t) T_const.
%   C_diag  - 1 x m vector of fluid flow rates c_i.
%   S_plus  - Index vector of positive-rate phases  (subset of 1:m)
%   S_minus - Index vector of negative-rate phases  (subset of 1:m)
%             S0 is inferred as the complement of S_plus and S_minus.
%
% OPTIONAL NAME-VALUE PAIRS
%   'n'          - Time discretisation points in [0,1].  Default: 100.
%   'WfMethod'   - 'expm' (default) | 'ode'
%   'ResMethod'  - 'adjugate' (default) | 'eigvec'
%   'imag_max'   - Upper imaginary cutoff for zero search.  Default: 10*pi.
%   'K_rect'     - Full Im-spacings covered by rectangle.  Default: 3.
%
% LAYOUT CONVENTION
%   Row vectors of length m*n use PHASE-MAJOR layout:
%     columns (q-1)*n+1 : q*n  hold values for phase q at t_1,...,t_n.
%   Applies to P0 (1 x m*n) and each row of E (N x m*n) where N=numel(eigenvalues).
%
% CONJUGATE PAIRS
%   eigenvalues stores Im >= 0 representatives only.  Contributions from
%   complex eigenvalues use  2*real(...)  throughout; conjugates are never
%   stored or iterated over.
%
% WORKFLOW  (Algorithm, Appendix C of the companion paper)
%   Step 1  model.computeGamma()          % gamma(t)
%   Step 2  model.findZeros()             % zeros of det(I - Wf(0,1,lambda))
%   Step 3  model.computeP0()             % P(u,0)
%   Step 4  model.computeEigenfunctions() % phi(t, lambda_k)
%   Step 5  model.computeDistribution()   % CDF, density, mean, variance
%   -or-    model.run()                   % all five in sequence
%
% REFERENCES
%   Margolius & O'Reilly (2016), Queueing Systems 82(1-2):43-73.
%   [companion paper, 2026]

    %% ----------------------------------------------------------------
    %  PROPERTIES - Model specification
    %% ----------------------------------------------------------------
    properties (SetAccess = private)
        T_fun       % @(t) -> m x m generator
        C_diag      % 1 x m fluid rate vector
        C           % m x m diagonal rate matrix
        m           % number of phases
        S_plus      % positive-rate phase indices
        S_minus     % negative-rate phase indices
        S_zero      % zero-rate phase indices (inferred)

        n           % discretisation points in [0,1]
        dt          % 1/n
        t_grid      % 1 x n:  t_j = (j-1)/n

        WfMethod    % 'expm' | 'ode'
        ResMethod   % 'adjugate' | 'eigvec'
        imag_max    % upper Im cutoff for zero search (default 10*pi)
        K_rect      % number of full Im-spacings covered by rectangle (default 3)
    end

    %% ----------------------------------------------------------------
    %  PROPERTIES - Computed results
    %% ----------------------------------------------------------------
    properties (SetAccess = private)
        % Step 1
        gamma       % n x m  (plain row-per-time layout, NOT phase-major)
        mu          % scalar mean drift
        isStable    % logical

        % Step 2
        eigenvalues % N x 1 complex (Im >= 0 only); N determined by imag_max

        % Step 3
        P0          % 1 x (m*n), phase-major

        % Step 4
        E           % N x (m*n), phase-major; complex.
                    % stores phi for real lk, 2*phi for complex lk (Im>0).
                    % real(exp(lk*x)*E(k,:)) gives correct contribution.

    end

    %% ----------------------------------------------------------------
    %  PROPERTIES - Internal cache
    %% ----------------------------------------------------------------
    % (No cached internal state: H_matrix is not stored between steps)

    %% ----------------------------------------------------------------
    %  CONSTRUCTOR
    %% ----------------------------------------------------------------
    methods (Access = public)

        function obj = PeriodicFluidModel(T_fun, C_diag, varargin)
        % PERIODICFLUIDMODEL  Construct model, validate inputs, set defaults.
        %
        % S_plus, S_minus, S_zero are inferred from signs of C_diag.
        % C_diag must be organised as [positive rates, negative rates, zero rates].
            validateattributes(T_fun,  {'function_handle'}, {}, 'PeriodicFluidModel','T_fun');
            validateattributes(C_diag, {'numeric'},{'vector','real'}, 'PeriodicFluidModel','C_diag');

            obj.T_fun  = T_fun;
            obj.C_diag = C_diag(:).';
            obj.m      = numel(C_diag);
            obj.C      = diag(obj.C_diag);

            % Infer phase partition from signs of C_diag
            obj.S_plus  = find(obj.C_diag > 0);
            obj.S_minus = find(obj.C_diag < 0);
            obj.S_zero  = find(obj.C_diag == 0);

            % Enforce canonical ordering: [+...+ -...- 0...0]
            % This is a hard requirement -- if C is misordered we cannot
            % know whether T is consistently ordered, and the results would
            % be silently wrong.
            bad_order = false;
            if ~isempty(obj.S_plus) && ~isempty(obj.S_minus)
                bad_order = bad_order || (max(obj.S_plus) > min(obj.S_minus));
            end
            if ~isempty(obj.S_minus) && ~isempty(obj.S_zero)
                bad_order = bad_order || (max(obj.S_minus) > min(obj.S_zero));
            end
            if ~isempty(obj.S_plus) && ~isempty(obj.S_zero)
                bad_order = bad_order || (max(obj.S_plus) > min(obj.S_zero));
            end
            if bad_order
                error('PeriodicFluidModel:ordering', ...
                    ['C_diag must be in canonical order [+...+ -...- 0...0].\n' ...
                     'Got: [%s]\n' ...
                     'S_plus=%s, S_minus=%s, S_zero=%s.\n' ...
                     'Reorder both C_diag and the rows/columns of T(t) consistently.'], ...
                    num2str(obj.C_diag), num2str(obj.S_plus), ...
                    num2str(obj.S_minus), num2str(obj.S_zero));
            end

            p = inputParser();
            addParameter(p,'n',         100,        @(x) isscalar(x) && x>0);
            addParameter(p,'WfMethod',  'expm',     @(x) ismember(x,{'expm','ode'}));
            addParameter(p,'ResMethod', 'adjugate', @(x) ismember(x,{'adjugate','eigvec'}));
            addParameter(p,'imag_max',  10*pi,      @(x) isscalar(x) && x>0);
            addParameter(p,'K_rect',    3,          @(x) isscalar(x) && x>0);
            parse(p, varargin{:});

            obj.n         = p.Results.n;
            obj.dt        = 1 / obj.n;
            obj.t_grid    = (0:obj.n-1) / obj.n;
            obj.WfMethod  = p.Results.WfMethod;
            obj.ResMethod = p.Results.ResMethod;
            obj.imag_max  = p.Results.imag_max;
            obj.K_rect    = p.Results.K_rect;
        end

    end

    %% ----------------------------------------------------------------
    %  PUBLIC - Five algorithm steps
    %% ----------------------------------------------------------------
    methods (Access = public)

        function obj = run(obj)
        % RUN  Execute all five steps in sequence.
            obj = obj.computeGamma();
            obj = obj.findZeros();
            obj = obj.computeP0();
            obj = obj.computeEigenfunctions();
            obj = obj.computeDistribution();
        end

        % ----- Step 1 -----
        function restoreFromStruct(obj, S)
        % RESTOREFROMSTRUCT  Restore computed quantities from a loaded struct.
        % Called by loadModel to populate private properties from a .mat file.
            if isfield(S, 'gamma_data') && ~isempty(S.gamma_data)
                obj.gamma    = S.gamma_data;
                obj.mu       = obj.dt * sum(S.gamma_data * obj.C_diag.');
                obj.isStable = (obj.mu < 0);
            end
            if isfield(S, 'eigenvalues') && ~isempty(S.eigenvalues)
                obj.eigenvalues = S.eigenvalues;
            end
            if isfield(S, 'P0') && ~isempty(S.P0)
                obj.P0 = S.P0;
            end
            if isfield(S, 'E') && ~isempty(S.E)
                obj.E = S.E;
            end
        end

        function obj = computeGamma(obj)
        % COMPUTEGAMMA  Solve for gamma(t), asymptotic periodic distribution
        % of the background process J(t).
        %
        % Solves:  d/dt gamma(t) = gamma(t) T(t),  gamma(t)*1 = 1,
        %          with periodic boundary condition gamma(0) = gamma(1).
        %
        % Method (ported from CyclicFluidQueuev17.get.ODEsol):
        %   1. Integrate dp/dt = T(t)'*p from t=0 to t=1000 to wash out
        %      transients; use the endpoint as initial condition.
        %   2. Integrate one period [0,1] with that IC; interpolate onto
        %      t_grid and row-normalise to enforce sum=1 at each time point.
        %   3. Compute mu = (1/n) * sum_j gamma(t_j)*C*1  (midpoint rule).
        %
        % Populates: obj.gamma (n x m), obj.mu, obj.isStable.

            T     = obj.T_fun;
            m_loc = obj.m;

            % ---- Step 1: wash out transients --------------------------------
            % Integrate the COLUMN ODE  dp/dt = T(t)' * p
            % (gamma row-vector ODE transposed to column form).
            IC = ones(m_loc, 1) / m_loc;   % uniform start
            rhs = @(t, p) T(t).' * p;

            opts_coarse = odeset('RelTol',1e-8,'AbsTol',1e-10);
            [~, p_long] = ode45(rhs, [0, 1000], IC, opts_coarse);
            IC_refined  = p_long(end, :).';
            IC_refined  = abs(IC_refined) / sum(abs(IC_refined));  % ensure normalised

            % ---- Step 2: integrate one period with refined IC ---------------
            opts_fine = odeset('RelTol',1e-10,'AbsTol',1e-12);
            [t_sol, p_sol] = ode45(rhs, [0, 1], IC_refined, opts_fine);

            % Interpolate onto t_grid (row = time, col = phase)
            pq = interp1(t_sol, p_sol, obj.t_grid.');   % n x m

            % Row-normalise: enforce gamma(t)*1 = 1 at each t_j
            row_sums = sum(pq, 2);                       % n x 1
            pq       = pq ./ row_sums;

            obj.gamma = pq;

            % ---- Step 3: stability ------------------------------------------
            % mu = integral_0^1 gamma(t)*C*1 dt  (midpoint-rule, dt=1/n)
            % gamma rows are at t_j = (j-1)/n; midpoint rule gives equal weights.
            obj.mu       = obj.dt * sum(pq * obj.C_diag.');   % scalar
            obj.isStable = (obj.mu < 0);

            if obj.isStable
                fprintf('computeGamma complete. mu = %.6f  (stable).\n', obj.mu);
            else
                warning('PeriodicFluidModel:unstable', ...
                    'Model is NOT stable: mu = %.6f >= 0.', obj.mu);
            end
        end

        % ----- Step 2 -----
        function obj = findZeros(obj)
        % FINDZEROS  Find ALL zeros of det(I - Wf(0,1,lambda)) with
        %            Re(lambda) <= 0 and Im(lambda) <= imag_max.
        %
        % Architecture:
        %   1. Rectangle R = [sigma, 0] x [eps_im, B] where
        %        B = (K_rect + 0.5) * max_j(2*pi/c_j).
        %      Argument principle gives exact count N_rect in R.
        %      ALL rectangle zeros are kept (argument-principle verified).
        %   2. Real zeros found by sign-change scan; count checked against
        %      |S_+| + 1 (expected: |S_+| negative reals + lambda=0).
        %   3. Above B up to imag_max: asymptotic seeds polished by fsolve.
        %      Each polished zero checked for proximity to its seed
        %      (guard against starting asymptotic regime too early).
        %      Per-family spacing checked for gaps after collection.
        %   4. No truncation: all verified zeros are stored.
        %
        % Populates: obj.eigenvalues (N x 1 complex, Im >= 0).

            % ---- Stability check -------------------------------------------
            if ~isempty(obj.gamma)
                if ~obj.isStable
                    error('PeriodicFluidModel:unstable', ...
                        'Model is not stable (mu = %.4f >= 0).', obj.mu);
                end
            else
                warning('PeriodicFluidModel:stabilityNotChecked', ...
                    ['gamma not yet computed; stability assumed. ' ...
                     'Call computeGamma first to verify.']);
            end

            % ---- Fixed parameters ------------------------------------------
            tol_zero  = 1e-6;   % acceptance tolerance on |detFunc|
            tol_dedup = 1e-3;   % uniquetol radius for deduplication
            n_grid    = 40;     % coarse minSVD grid points per axis
            eps_im    = 1e-6;   % shift bottom edge off real axis

            % ---- Asymptotic family data ------------------------------------
            [re_asym, im_step] = obj.asymptoticFamilies();
            nPlus  = numel(re_asym);
            sigma  = min(re_asym) - 2;   % initial left edge; expands if needed
            B      = (obj.K_rect + 0.5) * max(im_step);

            if obj.imag_max <= B
                warning('PeriodicFluidModel:imagMaxTooSmall', ...
                    'imag_max (%.4g) <= B (%.4g); raising imag_max to 2*B.', ...
                    obj.imag_max, B);
                obj.imag_max = 2 * B;
            end

            fprintf('findZeros: sigma=%.3f  B=%.3f  imag_max=%.3f\n', ...
                sigma, B, obj.imag_max);

            % ---- Real zeros ------------------------------------------------
            % findRealZeros expands leftward until |S_+|+1 zeros are found.
            real_zeros = obj.findRealZeros(sigma, tol_zero, tol_dedup);

            % ---- Complex zeros in rectangle --------------------------------
            rect   = [sigma, 0, eps_im, B];
            N_rect = obj.countZerosInRectangle(rect);
            fprintf('  Argument principle: %d complex zeros in rectangle.\n', N_rect);

            complex_zeros = obj.findZerosInRectangle(rect, N_rect, ...
                                n_grid, tol_zero, tol_dedup);

            % ---- Asymptotic zeros above B ----------------------------------
            % Pass verified rectangle zeros so asymptotic search can
            % compute an accurate shortfall count.
            rect_zeros = [real_zeros(:); complex_zeros(:)];
            asym_zeros = obj.findZerosAsymptotic(B, obj.imag_max, ...
                             im_step, re_asym, tol_zero, tol_dedup, rect_zeros);

            % ---- Detect coincident asymptotic families ----------------------
            % Pairs of families with same asymptotic real part will produce
            % genuinely distinct but close zeros at coincident heights.
            % We record these so spacing/dedup checks handle them correctly.
            % tol_re: threshold for 'same asymptotic real part'.
            % Use half the minimum gap between distinct re_asym values,
            % or 0.5 if all families share the same column.
            re_gaps = abs(diff(sort(unique(re_asym))));
            if isempty(re_gaps)
                tol_re = 0.5;
            else
                tol_re = min(re_gaps) / 2;
            end
            [coin_pairs, coin_heights] = obj.findCoincidentFamilies( ...
                re_asym, im_step, obj.imag_max, tol_re);

            % ---- Assemble all zeros ----------------------------------------
            % Rectangle zeros are kept unconditionally (arg-principle verified).
            % Asymptotic zeros added only if not already in rectangle set.
            new_asym  = obj.removeIfNearExisting(asym_zeros, rect_zeros, tol_dedup);
            all_zeros = [rect_zeros(:); new_asym(:)];
            fprintf('  rect_zeros: %d  (real: %d + complex: %d)\n', ...
                numel(rect_zeros), numel(real_zeros), numel(complex_zeros));
            fprintf('  asym_zeros: %d  new after removeIfNear: %d\n', ...
                numel(asym_zeros), numel(new_asym));
            all_zeros = obj.sortEigenvalues(all_zeros);

            % ---- Total count verification -----------------------------------
            % Expected total = real zeros + all complex zeros up to imag_max.
            % The rectangle/asymptotic split is an internal search detail;
            % the total expected count is independent of B.
            %   Real:    |S_+| + 1  (exact)
            %   Complex: sum_j floor(imag_max / im_step(j))
            n_real_exp  = numel(obj.S_plus) + 1;
            n_cx_exp    = sum(floor(obj.imag_max ./ im_step));
            n_expected  = n_real_exp + n_cx_exp;
            n_found     = numel(all_zeros);

            fprintf('  Total zeros: found %d, expected %d\n', ...
                n_found, n_expected);
            fprintf('    (real: %d,  complex up to imag_max: %d)\n', ...
                n_real_exp, n_cx_exp);

            if n_found > n_expected
                % More than expected: deduplicate with tight tolerance
                tol_tight = tol_zero * 10;
                all_zeros = obj.deduplicateZeros(all_zeros, tol_tight);
                fprintf('  After dedup: %d zeros\n', numel(all_zeros));
            elseif n_found < n_expected
                warning('PeriodicFluidModel:totalZeroCount', ...
                    ['Total zero count short: found %d, expected ~%d. ' ...
                     'Some zeros may be missing.'], n_found, n_expected);
            end

            % ---- Per-family spacing check -----------------------------------
            obj.checkFamilySpacing(all_zeros, re_asym, im_step, ...
                coin_pairs, coin_heights, B, tol_dedup);

            all_zeros       = obj.sortEigenvalues(all_zeros);
            obj.eigenvalues = all_zeros;
            fprintf('findZeros complete: %d eigenvalues stored.\n', ...
                numel(obj.eigenvalues));
        end

        % ----- Step 3 -----
        function obj = computeP0(obj)
        % COMPUTEP0  Solve for boundary mass P(u,0) (Algorithm steps 3a-3c).
        %
        % 3a. Build H = (dt) * sum_k contribution_k  where:
        %       contribution_k = buildBlockMatrix(Winc_k) * Rarray_k  (real lk)
        %                      = 2*real(...)                          (complex lk)
        %     H is (m*n) x (m*n), phase-major.
        %     Each column j of the block matrix is pre-multiplied by C.
        %
        % 3b. Extract H_{--}.  Solve augmented system:
        %       p_- * [H_{-+}  (I-H_{--})  1_{n*|S-|}] = [0...0  1]
        %     Recover p_0 = p_- * H_{-0}.
        %
        % 3c. Normalise: solve gamma(t) = P(t,0) * F for scalar alpha,
        %     where F is the lambda=0 contribution matrix.
        %
        % Populates: obj.P0 (1 x m*n, phase-major).
            % ---- Auto-compute prerequisites --------------------------------
            if isempty(obj.gamma)
                fprintf('computeP0: gamma not computed, running computeGamma.\n');
                obj.computeGamma();
            end
            if isempty(obj.eigenvalues)
                fprintf('computeP0: eigenvalues not computed, running findZeros.\n');
                obj.findZeros();
            end

            m  = obj.m;
            n  = obj.n;
            N  = m * n;    % total size

            % ----------------------------------------------------------------
            % LAYOUT NOTE: all internal computation uses TIME-MAJOR layout:
            %   block (i,j) occupies rows (i-1)*m+1:i*m, cols (j-1)*m+1:j*m
            %   i = initial time index, j = final time index
            % Only at the end do we convert to phase-major for storage.
            % ----------------------------------------------------------------

            % Helper: time-major index sets for each phase group
            % In time-major: for phase q at time i, row = (i-1)*m + q
            % So phase group S occupies rows [S, m+S, 2m+S, ...] interleaved.
            % We build the permutation from time-major to phase-major once.
            perm_t2p = reshape(reshape(1:N, m, n).', [], 1);  % time->phase major
            perm_p2t = zeros(N,1);
            perm_p2t(perm_t2p) = (1:N).';                     % phase->time major

            % ---- 3a: Build H in time-major layout ---------------------------
            fprintf('computeP0: building H matrix (%dx%d)...\n', N, N);
            H_tm = zeros(N, N);   % time-major

            for ki = 1:numel(obj.eigenvalues)
                lk   = obj.eigenvalues(ki);
                Winc = obj.computeWf_inc(lk);

                % create_block_matrix returns time-major layout directly
                BigW_tm = obj.create_block_matrix(Winc);   % N x N, time-major

                % Pre-multiply block-rows by C: block-row i (rows (i-1)*m+1:i*m)
                % gets left-multiplied by C (m x m).
                % Equivalent to kron(eye(n), C) * BigW_tm in time-major.
                BigCW_tm = kron(eye(n), obj.C) * BigW_tm;

                % Residue array at each t_j
                Rarray = obj.computeResidueArray(lk, Winc);

                % Post-multiply block-column j by Rarray(:,:,j):
                % block-col j occupies cols (j-1)*m+1:j*m.
                A_k_tm = zeros(N, N);
                for j = 1:n
                    col_idx = (j-1)*m + (1:m);
                    A_k_tm(:, col_idx) = BigCW_tm(:, col_idx) * Rarray(:,:,j);
                end

                % Accumulate with dt weight and conjugate handling
                if imag(lk) == 0
                    H_tm = H_tm + obj.dt * real(A_k_tm);
                else
                    H_tm = H_tm + obj.dt * 2 * real(A_k_tm);
                end

                if mod(ki, 10) == 0
                    fprintf('  Eigenvalue %d / %d done.\n', ki, numel(obj.eigenvalues));
                end
            end

            % Convert H to phase-major for block extraction
            H = H_tm(perm_t2p, perm_t2p);

            % ---- 3b: Solve for alpha * p_-(t) --------------------------------
            % In phase-major: phase q occupies rows/cols (q-1)*n+1:q*n.
            nm = numel(obj.S_minus) * n;

            idx_plus  = cell2mat(arrayfun(@(q) (q-1)*n+(1:n), ...
                obj.S_plus,  'UniformOutput', false));
            idx_minus = cell2mat(arrayfun(@(q) (q-1)*n+(1:n), ...
                obj.S_minus, 'UniformOutput', false));
            idx_zero  = cell2mat(arrayfun(@(q) (q-1)*n+(1:n), ...
                obj.S_zero,  'UniformOutput', false));

            H_mm = H(idx_minus, idx_minus);
            H_mp = H(idx_minus, idx_plus);
            H_m0 = H(idx_minus, idx_zero);
            np   = numel(idx_plus);

            % Solve p_- * [H_{-+}  (I-H_{--})  1] = [0  0  1]
            % Transposed: [H_{-+}'; (I-H_{--})'; 1] * p_-' = [0; 0; 1]
            lhs_3b = [H_mp.' ; (eye(nm) - H_mm).' ; ones(1, nm)];
            rhs_3b = [zeros(np, 1); zeros(nm, 1); 1];
            p_minus_alpha = (lhs_3b \ rhs_3b).';   % 1 x nm

            p_zero_alpha = p_minus_alpha * H_m0;    % 1 x nz

            % Assemble alpha*P0 in phase-major
            alpha_P0 = zeros(1, N);
            alpha_P0(idx_minus) = p_minus_alpha;
            alpha_P0(idx_zero)  = p_zero_alpha;

            % ---- 3c: Normalise using gamma(t) = P(t,0) * F ------------------
            fprintf('computeP0: computing normalisation F...\n');
            Winc0   = obj.computeWf_inc(0);
            BigW0_tm = obj.create_block_matrix(Winc0);
            BigCW0_tm = kron(eye(n), obj.C) * BigW0_tm;
            Rarray0 = obj.computeResidueArray(0, Winc0);

            F_tm = zeros(N, N);
            for j = 1:n
                col_idx = (j-1)*m + (1:m);
                F_tm(:, col_idx) = BigCW0_tm(:, col_idx) * Rarray0(:,:,j);
            end
            F_tm = obj.dt * real(F_tm);

            % Convert F to phase-major
            F = F_tm(perm_t2p, perm_t2p);

            % gamma in phase-major layout: 1 x N
            gamma_pm = zeros(1, N);
            for q = 1:m
                gamma_pm((q-1)*n + (1:n)) = obj.gamma(:,q).';
            end

            % Solve alpha * (alpha_P0 * F) = gamma_pm for scalar alpha
            lhs_alpha = alpha_P0 * F;
            alpha = real((lhs_alpha * gamma_pm.') / (lhs_alpha * lhs_alpha.'));

            obj.P0 = alpha * alpha_P0;

            fprintf('computeP0 complete. ||P0||=%.6f\n', norm(real(obj.P0)));
        end

        % ----- Step 4 -----
        function obj = computeEigenfunctions(obj)
        % COMPUTEEIGENFUNCTIONS  Compute phi(t, lambda_k) for each eigenvalue.
        %
        % phi(t, lambda_k) = integral_{t-1}^{t} P(u,0) C Wf(u,t,lk) du
        %                    * Res_{lk}[(I - Wf(t-1,t,lambda))^{-1}]
        %
        % Discretised: E(k,:) = P0_tm * A_k  (converted to phase-major)
        % where A_k is the same (m*n)x(m*n) time-major block matrix used
        % in computeP0, with (i,j) block = C * Wf(t_i,t_j,lk) * Res_{lk,j}.
        %
        % Real lk:    E(k,:) = real(P0_tm * A_k)  -> stored directly
        % Complex lk: E(k,:) = 2*real(P0_tm * A_k) -> captures conjugate pair
        %
        % A_k is recomputed fresh for each eigenvalue (storing all would
        % require numel(eigenvalues) copies of an (m*n)^2 matrix).
        %
        % Populates: obj.E (N x m*n, phase-major) where N=numel(eigenvalues).
            if isempty(obj.P0)
                fprintf('computeEigenfunctions: P0 not computed, running computeP0.\n');
                obj.computeP0();
            end

            m  = obj.m;
            n  = obj.n;
            N  = m * n;

            % Permutation: time-major -> phase-major (same as computeP0)
            perm_t2p = reshape(reshape(1:N, m, n).', [], 1);

            % Skip lambda=0: that eigenfunction is gamma(t), already stored.
            eig_idx = find(real(obj.eigenvalues) < 0);
            nK_eff  = numel(eig_idx);
            obj.E   = zeros(nK_eff, N);
            fprintf('computeEigenfunctions: computing %d eigenfunctions...\n', nK_eff);

            for ki_eff = 1:nK_eff
                ki   = eig_idx(ki_eff);
                lk   = obj.eigenvalues(ki);
                Winc = obj.computeWf_inc(lk);

                % Build A_k in time-major layout (identical to computeP0):
                %   block(i,j) = C * Wf(t_i,t_j,lk) * Res_{lk}(t_j)
                BigW_tm  = obj.create_block_matrix(Winc);
                BigCW_tm = kron(eye(n), obj.C) * BigW_tm;
                Rarray   = obj.computeResidueArray(lk, Winc);

                A_k_tm = zeros(N, N);
                for j = 1:n
                    col_idx = (j-1)*m + (1:m);
                    A_k_tm(:, col_idx) = BigCW_tm(:, col_idx) * Rarray(:,:,j);
                end

                % Convert A_k to phase-major
                A_k_pm = A_k_tm(perm_t2p, perm_t2p);

                % phi(t,lk) = P0 * A_k  (both phase-major)
                % P0 is 1x(m*n) phase-major; result is 1x(m*n) phase-major
                phi_pm = obj.P0 * A_k_pm;

                % Store: phi for real lk, 2*phi for complex lk
                if imag(lk) == 0
                    obj.E(ki_eff,:) = obj.dt * phi_pm;
                else
                    obj.E(ki_eff,:) = obj.dt * 2 * phi_pm;
                end

                if mod(ki_eff, 10) == 0
                    fprintf('  Eigenfunction %d / %d done.\n', ki_eff, nK_eff);
                end
            end

            fprintf('computeEigenfunctions complete.\n');
        end

        % ----- Step 5: query methods (compute on demand, no storage) -----

        function EX = meanFluid(obj)
        % MEANFLUID  E[X(t)], 1 x n vector over the period.
        %
        % E[X(t)] = sum_k (1/lambda_k) * phi(t,lambda_k) * 1
        %
        % For real lk:    contribution = (1/lk) * E(k,:) * 1
        % For complex lk: E(k,:) = 2*real(phi), so we need
        %                 2*real((1/lk) * phi) = 2*real((1/lk) * E(k,:)/2)
        %                 = real((1/lk) * E(k,:))
        % In both cases the contribution is real((1/lk) * E(k,:)) * 1,
        % since for real lk, real((1/lk)*E) = (1/lk)*E.
            if isempty(obj.E), obj.computeEigenfunctions(); end
            lam  = obj.eigenvalues(real(obj.eigenvalues) < 0);
            % Phase-sum: reshape E to (nK x n x m), sum over phase dim
            E_3d = reshape(obj.E, size(obj.E,1), obj.n, obj.m);
            Esum = reshape(sum(E_3d, 3), size(obj.E,1), obj.n);  % nK x n
            E_3d = reshape(obj.E, size(obj.E,1), obj.n, obj.m);
            Esum = reshape(sum(E_3d, 3), size(obj.E,1), obj.n);  % nK x n
            EX   = sum(real((1./lam(:)) .* Esum), 1);             % 1 x n
        end

        function VX = varFluid(obj)
        % VARFLUID  Var[X(t)], 1 x n vector over the period.
        %
        % Var[X(t)] = -sum_k (2/lambda_k^2) * phi(t,lambda_k)*1 - E[X]^2
        % Same real/complex handling as meanFluid: use real((-2/lk^2)*E)*1.
            if isempty(obj.E), obj.computeEigenfunctions(); end
            lam  = obj.eigenvalues(real(obj.eigenvalues) < 0);
            E_3d = reshape(obj.E, size(obj.E,1), obj.n, obj.m);
            Esum = reshape(sum(E_3d, 3), size(obj.E,1), obj.n);
            E_3d = reshape(obj.E, size(obj.E,1), obj.n, obj.m);
            Esum = reshape(sum(E_3d, 3), size(obj.E,1), obj.n);
            EX2  = sum(real((-2./lam(:).^2) .* Esum), 1);
            VX   = EX2 - obj.meanFluid().^2;
        end

        function sdX = stdFluid(obj)
        % STDFLUID  std(X(t)), 1 x n vector.
            sdX = sqrt(max(0, obj.varFluid()));
        end

        function Ptx = cdf(obj, x_vals)
        % CDF  Joint CDF P(t, X<=x, J=j) for all phases and times.
        %
        % INPUTS
        %   x_vals - Nx x 1 column vector of x values.
        %            Default: (0:0.02:10)'.
        % RETURNS
        %   Ptx - Nx x (m*n), phase-major layout.
        %
        % P(t,x) = sum_k exp(lk*x) * phi(t,lk) + gamma(t)
            if isempty(obj.E), obj.computeEigenfunctions(); end
            if nargin < 2, x_vals = (0:0.02:10).'; end
            x_vals   = x_vals(:);
            lam      = obj.eigenvalues(real(obj.eigenvalues) < 0);
            B        = exp(x_vals * lam(:).');          % Nx x nK
            gamma_pm = obj.gammaPhMajor();               % 1 x (m*n)
            % E stores phi for real lk, 2*phi for complex lk.
            % real(exp(lk*x) * E(k,:)) gives correct contribution in both cases.
            B   = exp(x_vals * lam(:).');        % Nx x nK
            Ptx = real(B * obj.E) + repmat(gamma_pm, numel(x_vals), 1);
        end

        function ptx = density(obj, x_vals)
        % DENSITY  Density d/dx P(t,x) for all phases and times.
        %
        % RETURNS ptx - Nx x (m*n), same layout as cdf.
            if isempty(obj.E), obj.computeEigenfunctions(); end
            if nargin < 2, x_vals = (0:0.02:10).'; end
            x_vals = x_vals(:);
            lam    = obj.eigenvalues(real(obj.eigenvalues) < 0);
            B      = exp(x_vals * lam(:).');             % Nx x nK
            % d/dx P = sum_k lk * exp(lk*x) * phi(t,lk)
            B   = exp(x_vals * lam(:).');
            ptx = real(B * (lam(:) .* obj.E));
        end

    end

    %% ----------------------------------------------------------------
    %  PUBLIC - Plots and display
    %% ----------------------------------------------------------------
    methods (Access = public)

        function obj = plotGamma(obj)
        % PLOTGAMMA  gamma_j(t): asymptotic periodic phase distribution over one period.
            if isempty(obj.gamma), obj = obj.computeGamma(); end
            figure;
            plot(obj.t_grid, obj.gamma, 'LineWidth', 2);
            xlabel('Time within period');  ylabel('Probability');
            title('Phases of background process');
            legend(arrayfun(@(j) sprintf('\\gamma_%d(t)',j), 1:obj.m, 'UniformOutput',false));
            grid on; grid minor;
        end

        function obj = plotDrift(obj)
        % PLOTDRIFT  Instantaneous drift gamma(t)*C*1 over one period.
            if isempty(obj.gamma), obj = obj.computeGamma(); end
            figure;
            plot(obj.t_grid, obj.gamma * obj.C_diag.', 'LineWidth', 2);
            xlabel('Time within period');  ylabel('Fluid level');
            title('Drift');  grid on; grid minor;
        end

        function obj = plotP0(obj)
        % PLOTP0  P_j(t,0) for all phases.
            if isempty(obj.P0), obj = obj.computeP0(); end
            figure; hold on;
            for q = 1:obj.m
                idx = (q-1)*obj.n + (1:obj.n);
                plot(obj.t_grid, real(obj.P0(idx)), 'LineWidth', 2);
            end
            xlabel('Time within period');  ylabel('Probability');
            title('Probability mass at zero');
            legend(arrayfun(@(j) sprintf('P_%d(t,0)',j), 1:obj.m, 'UniformOutput',false));
            grid on; grid minor; hold off;
        end

        function T = Toft(obj, t)
        % TOFT  Generator matrix at time t. Compatibility method for mysimc.
            T = obj.T_fun(t);
        end

        function c = flowRates(obj)
        % FLOWRATES  Fluid rate vector. Compatibility method for mysimc.
            c = obj.C_diag;
        end

        function obj = plotDensity(obj, x_vals, phases)
        % PLOTDENSITY  Contourf plots of d/dx P_j(t,x) for each phase.
        %
        % Reproduces the density figures in the companion paper (Tables 4,7).
        % One subplot per phase, x-axis = time within period, y-axis = fluid
        % level, color = density value.
        %
        % INPUTS
        %   x_vals - column vector of fluid levels.  Default: linspace(0,2,100)'.
        %   phases - which phases to plot.  Default: all phases.
            if isempty(obj.E), obj.computeEigenfunctions(); end
            if nargin < 2 || isempty(x_vals), x_vals = linspace(0, 2, 100).'; end
            if nargin < 3 || isempty(phases),  phases  = 1:obj.m; end

            x_vals = x_vals(:);
            ptx    = obj.density(x_vals);   % Nx x (m*n)

            np = numel(phases);
            % If single phase, use full figure; otherwise subplots
            if np == 1
                q   = phases(1);
                idx = (q-1)*obj.n + (1:obj.n);
                D   = real(ptx(:, idx));
                contourf(obj.t_grid, x_vals, D, 20, 'LineColor', 'none');
                colorbar;
                xlabel('Time within the period');
                ylabel('Fluid level');
                title(sprintf('Fluid level density phase %d', q));
            else
                nc = min(3, np);
                nr = ceil(np / nc);
                figure('Position', [100 100 300*nc 250*nr]);
                for pi_idx = 1:np
                    q   = phases(pi_idx);
                    idx = (q-1)*obj.n + (1:obj.n);
                    D   = real(ptx(:, idx));
                    subplot(nr, nc, pi_idx);
                    contourf(obj.t_grid, x_vals, D, 20, 'LineColor', 'none');
                    colorbar;
                    xlabel('Time within the period');
                    ylabel('Fluid level');
                    title(sprintf('Fluid level density phase %d', q));
                end
            end
        end

        function obj = plotMeanFluid(obj)
        % PLOTMEANFLUID  E[X(t)] +/- sigma(t).
            if isempty(obj.E), obj.computeEigenfunctions(); end
            figure;
            EX = obj.meanFluid();
            sdX = obj.stdFluid();
            plot(obj.t_grid, EX,       'b-',  'LineWidth',2); hold on;
            plot(obj.t_grid, EX + sdX, 'b--', 'LineWidth',1);
            plot(obj.t_grid, EX - sdX, 'b:',  'LineWidth',1);
            xlabel('Time within period');  ylabel('Fluid level');
            title('E[X(t)] as a function of time within the cycle');
            legend('E[X(t)]','E[X(t)]+\sigma(t)','E[X(t)]-\sigma(t)');
            grid on; grid minor; hold off;
        end

        function pdistcum = simulate(obj, n_events, n_runs, max_fluid, sim_mesh)
        % SIMULATE  Monte Carlo simulation of the cyclic fluid model.
        %
        % Adapted from doSim.m.  Uses the rejection method (Ross [27])
        % to simulate the time-varying CTMC and fluid level.
        %
        % INPUTS
        %   n_events  - events per run.         Default: 50000
        %   n_runs    - number of independent runs. Default: 20
        %   max_fluid - upper fluid level for CDF bins. Default: 5
        %   sim_mesh  - time discretisation for output. Default: 100
        %
        % RETURNS
        %   pdistcum - (max_fluid*stepx+1) x sim_mesh x m array of
        %              empirical CDFs P(X(t)<=x, J=j) for each phase j.
        %
        % NOTE: requires mysimc.m on the MATLAB path.
            if nargin < 2 || isempty(n_events),  n_events  = 50000; end
            if nargin < 3 || isempty(n_runs),    n_runs    = 20;    end
            if nargin < 4 || isempty(max_fluid), max_fluid = 5;     end
            if nargin < 5 || isempty(sim_mesh),  sim_mesh  = 100;   end

            phases  = obj.m;
            stepx   = 50;          % fluid bins per unit
            deltax  = 1/stepx;
            deltat  = 1/sim_mesh;
            maxfluidbin = max_fluid * stepx + 1;

            % Majorising rate for rejection method
            tempArray = zeros(phases, phases, sim_mesh);
            for j = 1:sim_mesh
                tempArray(:,:,j) = obj.T_fun((j-1)/sim_mesh);
            end
            lam = 1.1 * max(abs(tempArray(:)));
            clear tempArray

            % Initial state: phase 1, fluid 0, time 0
            mycont = [1 0 0];

            daytimeArray    = zeros(sim_mesh, 1);
            accumphaseStack = zeros(maxfluidbin, sim_mesh, phases);

            for ll = 1:n_runs
                fprintf('  Run %d / %d...\n', ll, n_runs);
                [cumRV, fluidLevel] = mysimc(n_events, mycont, obj, lam);

                lasttime = cumRV(end,1);
                xx  = 0:deltat:lasttime;
                eps_sm = 1e-11;
                smidge = cumsum(ones(size(cumRV(:,1)))) * eps_sm;
                x_t    = smidge + cumRV(:,1);

                % Interpolate fluid level and phase onto uniform time grid
                fl = interp1(x_t, fluidLevel, xx).';
                fl = floor(stepx * fl);
                if maxfluidbin < max(fl)
                    fl(fl > maxfluidbin-1) = maxfluidbin-1;
                end
                ph = interp1(x_t, cumRV(:,2), xx, 'previous').';

                % Fold time back onto [0,1) and bin
                roundtime = (xx - floor(xx)).';
                roundtime = floor(sim_mesh * roundtime) / sim_mesh;

                A = [ph roundtime fl];
                A(1,:) = [1 0 0];

                TT = array2table(A);
                G  = groupsummary(TT, ["A1","A2","A3"], ...
                    "IncludeEmptyGroups", true);
                GG = groupsummary(TT, "A2", ...
                    "IncludeEmptyGroups", true);

                theseCounts  = table2array(GG);
                daytimeArray = daytimeArray + theseCounts(:,2);

                AA = table2array(G);
                [myrows,~] = size(AA);
                fluidBuckets = myrows / (phases * sim_mesh);

                phaseStack = zeros(fluidBuckets, sim_mesh, phases);
                for j = 1:phases
                    indp   = AA(:,1) == j;
                    mymat  = AA(indp,:);
                    coltoreshape = mymat(:,4);
                    phaseStack(:,:,j) = reshape(coltoreshape, fluidBuckets, sim_mesh);
                end

                augment   = maxfluidbin - fluidBuckets;
                zeroblock = zeros(augment, sim_mesh, phases);
                phaseStack      = cat(1, phaseStack, zeroblock);
                accumphaseStack = accumphaseStack + phaseStack;
            end

            % Normalise by total time at each time bin
            daytimeArray  = sim_mesh * daytimeArray;
            daytimeRecip  = 1 ./ daytimeArray;
            for j = 1:phases
                accumphaseStack(:,:,j) = accumphaseStack(:,:,j) * diag(daytimeRecip);
            end

            % Cumulative sum over fluid levels -> empirical CDF
            pdistcum = zeros(maxfluidbin, sim_mesh, phases);
            for j = 1:phases
                tempdist = squeeze(accumphaseStack(:,:,j));
                pdistcum(:,:,j) = sim_mesh * cumsum(tempdist, 1);
            end
        end

        function plotCDFwithSim(obj, pdistcum, phase, x_levels, sim_mesh)
        % PLOTCDFWITHSIM  Overlay analytic CDF with simulation results.
        %
        % INPUTS
        %   pdistcum  - output of simulate()
        %   phase     - which phase to plot
        %   x_levels  - x thresholds to show.  Default: [0.2 0.4 0.6 0.8 1 5]
        %   sim_mesh  - time mesh used in simulate().  Default: 100
            if isempty(obj.E), obj.computeEigenfunctions(); end
            if nargin < 3 || isempty(phase),    phase    = 1; end
            if nargin < 4 || isempty(x_levels), x_levels = [0.2 0.4 0.6 0.8 1.0 5.0]; end
            if nargin < 5 || isempty(sim_mesh), sim_mesh = 100; end

            stepx    = 50;
            t_sim    = (0:sim_mesh-1) / sim_mesh;
            idx      = (phase-1)*obj.n + (1:obj.n);
            Ptx      = obj.cdf(x_levels(:));

            figure; hold on;
            plot(obj.t_grid, obj.gamma(:,phase), 'k', 'LineWidth', 2);
            lgd = {['$\gamma_{' num2str(phase) '}(t)$']};
            colors = lines(numel(x_levels));

            for k = 1:numel(x_levels)
                % Analytic (in legend)
                h_an = plot(obj.t_grid, Ptx(k, idx), '-', ...
                    'Color', colors(k,:), 'LineWidth', 2);
                % Simulation dots same color, excluded from legend
                % Empirically: round(x*stepx) gives correct bin alignment
                x_bin = round(x_levels(k) * stepx);
                x_bin = max(1, min(x_bin, size(pdistcum,1)));
                sim_cdf = pdistcum(x_bin, :, phase);
                h_sim = plot(t_sim, sim_cdf, '.', 'Color', colors(k,:));
                set(get(get(h_sim,'Annotation'),'LegendInformation'), ...
                    'IconDisplayStyle','off');
                lgd{end+1} = ['$P(x\leq ' sprintf('%.1f', x_levels(k)) ')$']; %#ok<AGROW>
            end

            xlabel('Time within period');
            ylabel('Probability');
            title(sprintf('Fluid level distribution phase %d', phase));
            leg = legend(lgd);
            set(leg, 'Interpreter', 'latex', 'Location', 'best');
            grid on; grid minor; hold off;
        end

        function obj = plotCDF(obj, phase, x_levels)
        % PLOTCDF  P_j(t, x<=xk) for a given phase and several x thresholds.
            if isempty(obj.E), obj.computeEigenfunctions(); end
            if nargin < 2, phase    = 1; end
            if nargin < 3, x_levels = [0.2 0.4 0.6 0.8 1.0 5.0]; end
            idx = (phase-1)*obj.n + (1:obj.n);
            Ptx = obj.cdf(x_levels(:));
            figure; hold on;
            plot(obj.t_grid, obj.gamma(:,phase), 'k', 'LineWidth',2);
            lgd = {sprintf('\\gamma_%d(t)', phase)};
            for k = 1:numel(x_levels)
                plot(obj.t_grid, Ptx(k, idx), 'LineWidth', 2);
                lgd{end+1} = sprintf('$P(x\\leq %.1f)$', x_levels(k)); %#ok<AGROW>
            end
            xlabel('Time within period');  ylabel('Probability');
            title(sprintf('Fluid level distribution phase %d', phase));
            leg = legend(lgd);  set(leg,'Interpreter','latex');
            grid on; grid minor; hold off;
        end

        function obj = plotPhasePlot(obj, realRange, imagRange, gridSize)
        % PLOTPHASEPLOT  Phase portrait of det(I - Wf(0,1,lambda)).
        % Requires Wegert''s PhasePlot function on the MATLAB path.
        %
        % NOTE: In Wegert''s pltphase.m, the line "axis off" must be
        % commented out for axis labels and tick marks to display correctly.
        % The distributed version of pltphase.m included with this code
        % has this modification applied.
            % eigenvalues not strictly required for portrait; computed if available
            if nargin < 4, gridSize = [300, 600]; end

            realVec = linspace(realRange(1), realRange(2), gridSize(1));
            imagVec = linspace(imagRange(1), imagRange(2), gridSize(2));
            [Re, Im] = meshgrid(realVec, imagVec);
            z    = Re.' + 1i*Im.';
            Ts   = obj.buildTarray();
            myC  = obj.dt * obj.C;
            Ieye = eye(obj.m);
            w    = zeros(gridSize(1), gridSize(2));

            parfor j = 1:gridSize(1)
                w_row = zeros(1, gridSize(2));   % local; avoids parfor slicing error
                for k = 1:gridSize(2)
                    tmp = Ieye;
                    zmC = -z(j,k) * myC;
                    for r = 1:obj.n
                        tmp = tmp * expm(Ts(:,:,r) + zmC);
                    end
                    w_row(k) = det(Ieye - tmp);
                end
                w(j,:) = w_row;
            end

            % Capture current figure before PhasePlot runs, in case
            % a different version creates a new figure internally.
            hfig = gcf;
            PhasePlot(z, w, 'c');
            % Ensure we are working on the correct figure/axes after PhasePlot
            figure(hfig);
            hax = gca;
            obj.formatPiAxis(imagRange, hax);
            xlim(hax, realRange);
            xlabel(hax, '$\Re\{z\}$', 'Interpreter', 'latex');
            ylabel(hax, '$\Im\{z\}$', 'Interpreter', 'latex');
        end

        function disp(obj)
        % DISP  Summary of model state.
            fprintf('PeriodicFluidModel  (m=%d phases, n=%d time points)\n', obj.m, obj.n);
            fprintf('  S+: [%s]   S-: [%s]   S0: [%s]\n', ...
                num2str(obj.S_plus), num2str(obj.S_minus), num2str(obj.S_zero));
            fprintf('  C:  [%s]\n', num2str(obj.C_diag));
            fprintf('  S+: [%s]  S-: [%s]  S0: [%s]\n', ...
                num2str(obj.S_plus), num2str(obj.S_minus), num2str(obj.S_zero));
            fprintf('  WfMethod: %-6s  ResMethod: %-9s  imag_max: %.4g  K_rect: %d\n', ...
                obj.WfMethod, obj.ResMethod, obj.imag_max, obj.K_rect);
            if ~isempty(obj.isStable)
                if obj.isStable
                    fprintf('  Stable: YES   mu = %.6f\n', obj.mu);
                else
                    fprintf('  Stable: NO (mu = %.6f, check inputs)\n', obj.mu);
                end
            end
            fields = {'gamma','eigenvalues','P0','E'};
            labels = {'gamma','zeros','P0','eigenfuncs'};
            done   = labels(~cellfun(@(f) isempty(obj.(f)), fields));
            if isempty(done)
                fprintf('  Nothing computed yet. Call run() or individual step methods.\n');
            else
                fprintf('  Computed: %s\n', strjoin(done,', '));
            end
        end

    end

    %% ----------------------------------------------------------------
    %  PRIVATE - Wf computation
    %% ----------------------------------------------------------------
    methods (Access = private)

        function W = computeWf_inc_expm(obj, lambda)
        % W(:,:,j) = expm( dt*(T(t_mid_j) - lambda*C) ),  t_mid_j=(j-0.5)/n.
            W  = zeros(obj.m, obj.m, obj.n);
            lC = lambda * obj.C;
            for j = 1:obj.n
                W(:,:,j) = expm(obj.dt * (obj.T_fun((j-0.5)/obj.n) - lC));
            end
        end

        function W = computeWf_inc_ode(obj, lambda)
        % Integrate d/dt Wf = Wf*(T(t)-lambda*C) on each sub-interval.
            W    = zeros(obj.m, obj.m, obj.n);
            lC   = lambda * obj.C;
            m2   = obj.m^2;
            rhs  = @(t,w) reshape(reshape(w,obj.m,obj.m)*(obj.T_fun(t)-lC), m2, 1);
            opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
            for j = 1:obj.n
                t0 = (j-1)*obj.dt;  t1 = j*obj.dt;
                [~,ws] = ode45(rhs, [t0,t1], eye(obj.m,obj.m), opts);
                W(:,:,j) = reshape(ws(end,:), obj.m, obj.m);
            end
        end

        function W = computeWf_inc(obj, lambda)
            switch obj.WfMethod
                case 'expm', W = obj.computeWf_inc_expm(lambda);
                case 'ode',  W = obj.computeWf_inc_ode(lambda);
            end
        end

        function M = buildMonodromy(obj, Winc)
        % Wf(0,1,lambda) = W_1 * W_2 * ... * W_n.
            M = eye(obj.m);
            for j = 1:obj.n
                M = M * Winc(:,:,j);
            end
        end

    end

    %% ----------------------------------------------------------------
    %  PRIVATE - Block matrix assembly (ported from CyclicFluidQueuev17)
    %% ----------------------------------------------------------------
    methods (Access = private)

        function BigMat = buildBlockMatrix(obj, Winc)
        % BUILDBLOCKMATRIX  (n*m)x(n*m) cyclic product matrix, phase-major.
        %
        % Calls create_block_matrix (time-major) then permuteBB (-> phase-major).
            BigMat = obj.permuteBB( obj.create_block_matrix(Winc) );
        end

        function M = create_block_matrix(obj, Warray)
        % CREATE_BLOCK_MATRIX  Time-major cyclic block matrix.
        %
        % Block (i,j) = Wf(t_i, t_j):
        %   i<=j : W_i * ... * W_{j-1}     (forward product)
        %   i >j : W_i * ... * W_n * W_1 * ... * W_{j-1}  (wrap-around)
        %   diagonal (i=j): full one-period product (from index i)
        %
        % Port of create_block_matrix from CyclicFluidQueuev17.
            [phases, ~, n] = size(Warray);
            N = n * phases;
            M = zeros(N);

            % First column: blocks computed bottom-to-top.
            block_col        = zeros(phases, phases, n);
            block_col(:,:,n) = Warray(:,:,n);
            for j = n-1:-1:1
                block_col(:,:,j) = Warray(:,:,j) * block_col(:,:,j+1);
            end

            for i = 1:n
                col_idx = (i-1)*phases + (1:phases);
                supdcol = mod(i-2, n) + 1;
                block_col(:,:,supdcol) = Warray(:,:,supdcol);
                M(:,col_idx) = obj.towerToCol(block_col);
                block_col    = pagemtimes(block_col, Warray(:,:,i));
            end
        end

        function BBfinal = permuteBB(obj, BB)
        % PERMUTEBB  Time-major -> phase-major reordering.
        %
        % Port of permuteBB from CyclicFluidQueuev17.
            permOrder = reshape(reshape(1:(obj.m*obj.n), obj.m, obj.n).', [], 1);
            BBfinal   = BB(permOrder, permOrder);
        end

        function BlockTower = colToTower(~, col)
        % COLTOTOWER  (n*phases) x phases  ->  phases x phases x n.
        %
        % Port of ColToTower from CyclicFluidQueuev17.
            phases     = size(col, 2);
            BlockTower = permute(reshape(col, phases, [], phases), [1, 3, 2]);
        end

        function col = towerToCol(~, BlockTower)
        % TOWERTOCOL  phases x phases x n  ->  (n*phases) x phases.
        %
        % Port of towerToCol from CyclicFluidQueuev17.
            col = reshape(permute(BlockTower, [1, 3, 2]), [], size(BlockTower,1));
        end

    end

    %% ----------------------------------------------------------------
    %  PRIVATE - Residue computation (stubs)
    %% ----------------------------------------------------------------
    methods (Access = private)

        function R = computeResidue_adjugate(obj, lambda, Winc)
        % COMPUTERESIDUE_ADJUGATE
        %
        % Res_{lambda_k}[(I-Wf(t-1,t,lambda))^{-1}] = a(lk)*adj(I-Wf(0,1,lk))
        %
        % where a(lk) = 1 / (d/dlambda det(I-Wf(0,1,lambda))|_{lk}).
        %
        % a(lk) is computed once via complex-step finite difference on detFunc.
        % adj(I-Wf(0,1,lk)) is computed from the monodromy matrix.
        % By Lemma 2, det does not depend on t so a(lk) is t-independent.
        % The adjugate DOES depend on t; see computeResidueArray for per-t
        % computation.
            M     = obj.buildMonodromy(Winc);
            A     = eye(obj.m) - M;
            % Complex-step derivative for high accuracy
            h     = 1e-8;
            d_det = (obj.detFunc(lambda + 1i*h) - obj.detFunc(lambda - 1i*h)) / (2i*h);
            if abs(d_det) < eps
                warning('PeriodicFluidModel:zeroDeriv', ...
                    'Near-zero derivative at lambda=%.4f%+.4fi; residue may be inaccurate.', ...
                    real(lambda), imag(lambda));
                R = zeros(obj.m);
                return
            end
            a_lk = 1 / d_det;
            R    = a_lk * adjoint(A);
        end

        function R = computeResidue_eigvec(obj, lambda, Winc)
        % R = -u*v' / (v' * dWf/dlambda * u)  (Corollary 3).
        % Port of compute_residue_projection from CyclicFluidQueuev17.
            R = zeros(obj.m);   % STUB
        end

        function R = computeResidue(obj, lambda, Winc)
            switch obj.ResMethod
                case 'adjugate', R = obj.computeResidue_adjugate(lambda, Winc);
                case 'eigvec',   R = obj.computeResidue_eigvec(lambda, Winc);
            end
        end

        function Rarray = computeResidueArray(obj, lambda, Winc)
        % COMPUTERESIDUEARRAY
        %
        % Rarray(:,:,j) = a(lk) * adj(I - Wf(t_j-1, t_j, lambda_k))
        % for j = 1,...,n.
        %
        % Wf(t_j-1, t_j, lk) is the one-period product starting at t_j:
        %   Wf(t_j-1, t_j) = W_j * W_{j+1} * ... * W_n * W_1 * ... * W_{j-1}
        % i.e. the cyclic product starting from index j.
        %
        % a(lk) is t-independent (Lemma 2) and computed once from Wf(0,1,lk).
            % Compute a(lk) once from monodromy
            h     = 1e-8;
            d_det = (obj.detFunc(lambda + 1i*h) - obj.detFunc(lambda - 1i*h)) / (2i*h);
            if abs(d_det) < eps
                warning('PeriodicFluidModel:zeroDeriv', ...
                    'Near-zero derivative at lambda=%.4f%+.4fi.', ...
                    real(lambda), imag(lambda));
                Rarray = zeros(obj.m, obj.m, obj.n);
                return
            end
            a_lk = 1 / d_det;

            Rarray = zeros(obj.m, obj.m, obj.n);
            for j = 1:obj.n
                % Build Wf(t_j-1, t_j, lk) = cyclic product starting at j
                % = W_j * W_{j+1} * ... * W_n * W_1 * ... * W_{j-1}
                M = eye(obj.m);
                for r = 0:obj.n-1
                    idx = mod(j-1+r, obj.n) + 1;
                    M   = M * Winc(:,:,idx);
                end
                Rarray(:,:,j) = a_lk * adjoint(eye(obj.m) - M);
            end
        end

    end

    %% ----------------------------------------------------------------
    %  PRIVATE - Zero-finding helpers
    %% ----------------------------------------------------------------
    methods (Access = private)

        function [re_asym, im_step] = asymptoticFamilies(obj)
        % ASYMPTOTICFAMILIES  Real parts and Im spacings of asymptotic zero
        % columns (Remarks 3-5), one entry per S_plus phase.
        %
        % re_asym(j) = mean(T_jj(t)) / c_j   (real part of asymptotic column)
        % im_step(j) = 2*pi / c_j             (imaginary spacing within column)
        %
        % Mean diagonal is approximated by the midpoint-rule average over
        % the discretisation used for Tarray.
            nPlus   = numel(obj.S_plus);
            re_asym = zeros(1, nPlus);
            im_step = zeros(1, nPlus);
            % Midpoint-rule average of T_{jj}(t) over one period
            Tbar = zeros(obj.m, 1);
            for r = 1:obj.n
                Td   = diag(obj.T_fun((r - 0.5)/obj.n));
                Tbar = Tbar + Td;
            end
            Tbar = Tbar / obj.n;
            for idx = 1:nPlus
                j          = obj.S_plus(idx);
                cj         = obj.C_diag(j);
                re_asym(idx) = Tbar(j) / cj;
                im_step(idx) = 2*pi    / cj;
            end
        end

        function lam_asym = asymptoticZeros(obj, K)
        % ASYMPTOTICZEROES  Grid of asymptotic estimates, Im >= 0.
        %
        % Returns candidates re_asym(j) + k * im_step(j) * i
        % for each S_plus family j and k = 0, 1, ..., K.
            [re_asym, im_step] = obj.asymptoticFamilies();
            nPlus    = numel(re_asym);
            lam_asym = zeros(nPlus * (K+1), 1);
            idx = 0;
            for j = 1:nPlus
                for k = 0:K
                    idx = idx + 1;
                    lam_asym(idx) = re_asym(j) + k * im_step(j) * 1i;
                end
            end
        end

        % -----------------------------------------------------------------
        function f = detFunc(obj, lambda)
        % det(I - Wf(0,1,lambda)).
            Winc = obj.computeWf_inc(lambda);
            M    = obj.buildMonodromy(Winc);
            f    = det(eye(obj.m) - M);
        end

        function f = minSVD(obj, lambda)
        % Smallest singular value of (I - Wf(0,1,lambda)).
            Winc = obj.computeWf_inc(lambda);
            M    = obj.buildMonodromy(Winc);
            f    = min(svd(eye(obj.m) - M));
        end

        function deriv = estimateDerivative(obj, lambda)
        % |d/dlambda det(I-Wf(0,1,lambda))| by centred finite difference.
            h  = 1e-8;
            dr = (obj.detFunc(lambda+h)    - obj.detFunc(lambda-h))    / (2*h);
            di = (obj.detFunc(lambda+1i*h) - obj.detFunc(lambda-1i*h)) / (2i*h);
            deriv = abs(dr) + abs(di);
        end

        % -----------------------------------------------------------------
        function N = countZerosInRectangle(obj, rect)
        % COUNTZEROSINRECTANGLE  Argument principle on a rectangular contour.
        %
        % rect = [re_min, re_max, im_min, im_max]
        %
        % Traverses the boundary anticlockwise, accumulates the change in
        % arg(detFunc), divides by 2*pi to get the winding number = zero count.
        % Uses 500 points per edge for reliable argument tracking.
            re_min = rect(1);  re_max = rect(2);
            im_min = rect(3);  im_max = rect(4);
            np     = 500;

            % Four edges: bottom, right, top (reversed), left (reversed)
            bottom = linspace(re_min, re_max, np) + 1i*im_min;
            right  = re_max + 1i*linspace(im_min, im_max, np);
            top    = linspace(re_max, re_min, np) + 1i*im_max;
            left   = re_min + 1i*linspace(im_max, im_min, np);

            contour_pts = [bottom, right(2:end), top(2:end), left(2:end)];

            % Evaluate detFunc along contour
            f_vals = arrayfun(@(z) obj.detFunc(z), contour_pts);

            % Accumulate argument change carefully to avoid 2*pi jumps
            darg   = diff(unwrap(angle(f_vals)));
            total  = sum(darg);
            N      = round(total / (2*pi));
        end

        % -----------------------------------------------------------------
        function zeros_found = findZerosInRectangle(obj, rect, N_expected, ...
                                    n_grid, tol_zero, tol_dedup)
        % FINDZEROSINRECTANGLE  Find all complex zeros inside rect.
        %
        % Zeros are accepted only if strictly inside the rectangle bounds
        % (im_min <= Im <= im_max) so that real zeros or zeros from outside
        % the contour cannot consume slots counted by the argument principle.
        %
        % Strategy:
        %   1. Coarse minSVD grid over rectangle interior.
        %   2. Polish each local minimum with fminsearch on minSVD.
        %   3. Accept only if inside rect bounds and quality check passes.
        %   4. Deduplicate.
        %   5. If count < N_expected, subdivide into quadrants and recurse.

            re_min = rect(1);  re_max = rect(2);
            im_min = rect(3);  im_max = rect(4);
            N_complex = N_expected;

            re_vec = linspace(re_min, re_max, n_grid);
            im_vec = linspace(im_min, im_max, n_grid);
            [RE, IM] = meshgrid(re_vec, im_vec);

            % Evaluate minSVD on coarse grid
            svd_grid = zeros(n_grid, n_grid);
            for ii = 1:n_grid
                for jj = 1:n_grid
                    svd_grid(ii,jj) = obj.minSVD(RE(ii,jj) + 1i*IM(ii,jj));
                end
            end

            % Find local minima on the grid
            candidates = obj.extractLocalMinima(RE, IM, svd_grid);

            % Polish each candidate
            zeros_found = [];
            opts = optimset('Display','off','TolX',1e-10,'TolFun',1e-12,'MaxIter',2000);
            for k = 1:size(candidates,1)
                z0  = candidates(k);
                obj_fun = @(v) obj.minSVD(v(1) + 1i*v(2));
                v_opt = fminsearch(obj_fun, [real(z0), imag(z0)], opts);
                z_opt = v_opt(1) + 1i*v_opt(2);

                if real(z_opt) >= re_min && real(z_opt) <= re_max && ...
                   imag(z_opt) >= im_min && imag(z_opt) <= im_max
                    fval  = abs(obj.detFunc(z_opt));
                    deriv = obj.estimateDerivative(z_opt);
                    if fval < tol_zero || (deriv > 0 && fval/deriv < 0.1)
                        zeros_found = [zeros_found; z_opt]; %#ok<AGROW>
                    end
                end
            end

            % Deduplicate
            zeros_found = obj.deduplicateZeros(zeros_found, tol_dedup);

            % Verify count; subdivide if short
            n_found = numel(zeros_found);
            fprintf('  Rectangle [%.2f,%.2f]x[%.2fi,%.2fi]: found %d / %d zeros.\n', ...
                re_min, re_max, im_min, im_max, n_found, N_complex);

            if n_found < N_complex
                fprintf('  Count mismatch -- subdividing rectangle.\n');
                zeros_found = obj.subdivideRectangle(rect, N_complex, ...
                                    zeros_found, n_grid, tol_zero, tol_dedup);
            end
        end

        % -----------------------------------------------------------------
        function zeros_found = subdivideRectangle(obj, rect, N_expected, ...
                                    zeros_so_far, n_grid, tol_zero, tol_dedup)
        % SUBDIVIDIRECTANGLE  Bisect rect into quadrants, recount, recurse.
        %
        % Only recurses into a quadrant if its argument-principle count
        % exceeds the number of already-found zeros in that quadrant.
            re_min = rect(1);  re_max = rect(2);
            im_min = rect(3);  im_max = rect(4);
            re_mid = (re_min + re_max) / 2;
            im_mid = (im_min + im_max) / 2;

            quads = {[re_min, re_mid, im_min, im_mid], ...
                     [re_mid, re_max, im_min, im_mid], ...
                     [re_min, re_mid, im_mid, im_max], ...
                     [re_mid, re_max, im_mid, im_max]};

            zeros_found = zeros_so_far;
            for q = 1:4
                qrect  = quads{q};
                N_q    = obj.countZerosInRectangle(qrect);
                in_q   = zeros_found(real(zeros_found) >= qrect(1) & ...
                                     real(zeros_found) <= qrect(2) & ...
                                     imag(zeros_found) >= qrect(3) & ...
                                     imag(zeros_found) <= qrect(4));
                if N_q > numel(in_q)
                    new_z = obj.findZerosInRectangle(qrect, N_q, ...
                                n_grid, tol_zero, tol_dedup);
                    zeros_found = obj.deduplicateZeros( ...
                                    [zeros_found; new_z], tol_dedup);
                end
            end
        end

        % -----------------------------------------------------------------
        function real_zeros = findRealZeros(obj, sigma_init, tol_zero, tol_dedup)
        % FINDREALZEROS  Find |S_+|+1 real zeros by sign-change scan with
        % optional leftward expansion.
        %
        % Logic:
        %   1. Scan [sigma, -origin_guard] for sign changes in real(detFunc).
        %   2. For each bracket, run fzero then accept if:
        %        (a) |detFunc| < tol_zero  (strict), OR
        %        (b) fval/deriv < 0.1  (fallback for poor expm conditioning)
        %      The fallback is safe here because a sign-change bracket
        %      guarantees a true zero exists -- the ratio test just handles
        %      the case where large |lambda| makes detFunc inaccurate.
        %   3. Expansion is driven by SIGN-CHANGE count, not zero count:
        %      if we have enough brackets we already have all the zeros, so
        %      we never wander leftward unnecessarily.
        %   4. lambda=0 is checked explicitly and always added.
            n_neg_expected = numel(obj.S_plus);   % negative real zeros expected
            n_expected     = n_neg_expected + 1;  % +1 for lambda=0
            origin_guard   = 10 * tol_dedup;
            sigma_step     = 2;
            max_expand     = 20;
            opts           = optimset('Display','off','TolX',1e-12);

            real_zeros = [];
            sigma      = sigma_init;

            % Always add lambda=0 explicitly
            if abs(obj.detFunc(0)) < tol_zero
                real_zeros = [real_zeros; 0];
                fprintf('  Zero at origin confirmed.\n');
            end

            for expand = 0 : floor(max_expand / sigma_step)
                x_grid  = linspace(sigma, -origin_guard, 500);
                f_vals  = arrayfun(@(x) real(obj.detFunc(x)), x_grid);
                changes = find(diff(sign(f_vals)) ~= 0);

                n_brackets = numel(changes);
                fprintf('  Scan [%.2f, %.2f]: %d sign changes found (need %d).\n', ...
                    sigma, -origin_guard, n_brackets, n_neg_expected);

                % Process every bracket regardless -- derivative fallback
                % handles poor conditioning without discarding real zeros.
                new_zeros = [];
                for k = 1:n_brackets
                    try
                        z     = fzero(@(x) real(obj.detFunc(x)), ...
                                      [x_grid(changes(k)), x_grid(changes(k)+1)], opts);
                        fval  = abs(obj.detFunc(z));
                        deriv = obj.estimateDerivative(z);
                        accepted = fval < tol_zero || ...
                                   (deriv > 0 && fval/deriv < 0.1);
                        if accepted
                            if isempty(real_zeros) && isempty(new_zeros) || ...
                                    min(abs(z - [real_zeros(:); new_zeros(:)])) > tol_dedup
                                new_zeros = [new_zeros; z]; %#ok<AGROW>
                                fprintf('  Real zero: %.6f  |det|=%.2e  deriv=%.2e\n', ...
                                    z, fval, deriv);
                            end
                        else
                            fprintf('  Rejected: %.6f  |det|=%.2e  deriv=%.2e  ratio=%.2e\n', ...
                                z, fval, deriv, fval/max(deriv,eps));
                        end
                    catch ME
                        fprintf('  fzero failed in bracket %d: %s\n', k, ME.message);
                    end
                end
                real_zeros = [real_zeros; new_zeros]; %#ok<AGROW>

                % Expansion decision: based on bracket count, not zero count.
                % If we have enough sign changes, all zeros are bracketed --
                % no need to look further left even if a quality check failed.
                if n_brackets >= n_neg_expected
                    fprintf('  Enough brackets found (%d >= %d). Stopping scan.\n', ...
                        n_brackets, n_neg_expected);
                    break;
                end

                % Genuinely short on sign changes -- expand leftward.
                sigma = sigma - sigma_step;
                fprintf('  Only %d/%d brackets. Expanding left to sigma=%.2f\n', ...
                    n_brackets, n_neg_expected, sigma);
            end

            if numel(real_zeros) < n_expected
                warning('PeriodicFluidModel:realZeroCount', ...
                    ['Found %d / %d real zeros. ' ...
                     'Brackets found: check rejected candidates above.'], ...
                    numel(real_zeros), n_expected);
            else
                fprintf('  Real zeros: found %d / %d. OK.\n', ...
                    numel(real_zeros), n_expected);
            end

            % Sort descending (least negative first)
            real_zeros = sort(real(real_zeros(:)), 'descend');
        end

        % -----------------------------------------------------------------
        function asym_zeros = findZerosAsymptotic(obj, B, imag_max, ...
                                    im_step, re_asym, tol_zero, tol_dedup, rect_zeros)
        % FINDZEROSASYMPTOTIC  Zeros above B up to imag_max using seeds + fsolve.
        %
        % For coincident families (same re_asym): seeds perturbed by +/-delta
        % so fsolve can find two distinct zeros.  When both converge to the same
        % point, the found zero is duplicated.  No deduplication is performed
        % here -- deduplication only makes sense when found > expected.
            nPlus    = numel(re_asym);
            opts     = optimset('Display','off','TolX',1e-10,'TolFun',1e-12);
            min_step = min(im_step);

            % ---- Seed perturbations for coincident families ----------------
            re_perturb     = zeros(1, nPlus);
            fallback_delta = 0.1;
            for j1 = 1:nPlus
                for j2 = j1+1:nPlus
                    if abs(re_asym(j1) - re_asym(j2)) < 1e-6
                        delta    = fallback_delta;
                        near_col = rect_zeros( ...
                            abs(real(rect_zeros) - re_asym(j1)) < 1.0 & ...
                            imag(rect_zeros) > 0);
                        if numel(near_col) >= 2
                            [~, sidx] = sort(imag(near_col), 'descend');
                            near_col  = near_col(sidx);
                            for zi = 1:numel(near_col)-1
                                if abs(imag(near_col(zi)) - imag(near_col(zi+1))) < ...
                                        min_step * 0.5
                                    d = abs(real(near_col(zi)) - real(near_col(zi+1)));
                                    if d > 1e-10, delta = d; end
                                    break;
                                end
                            end
                        end
                        fprintf('  Coincident families %d & %d: delta=%.6f\n', j1, j2, delta);
                        re_perturb(j1) =  delta;
                        re_perturb(j2) = -delta;
                    end
                end
            end

            % ---- Pass 1: fsolve from asymptotic seeds ----------------------
            asym_zeros = [];
            start_k    = ceil(B ./ im_step);
            for j = 1:nPlus
                k = start_k(j);
                while true
                    seed = (re_asym(j) + re_perturb(j)) + k * im_step(j) * 1i;
                    if imag(seed) > imag_max, break; end
                    unpert_seed = re_asym(j) + k * im_step(j) * 1i;
                    try
                        obj_fun = @(v) [real(obj.detFunc(v(1)+1i*v(2))); ...
                                        imag(obj.detFunc(v(1)+1i*v(2)))];
                        v_opt = fsolve(obj_fun, [real(seed), imag(seed)], opts);
                        z_opt = v_opt(1) + 1i*v_opt(2);
                        fval  = abs(obj.detFunc(z_opt));
                        deriv = obj.estimateDerivative(z_opt);
                        good  = fval < tol_zero || (deriv > 0 && fval/deriv < 0.1);
                        near  = abs(z_opt - unpert_seed) < 0.5 * im_step(j);
                        if good && near && imag(z_opt) > B && real(z_opt) <= 0.1
                            asym_zeros = [asym_zeros; z_opt]; %#ok<AGROW>
                        elseif good && ~near
                            fprintf('  Wandered: fam=%d row=%d zero=%.4f%+.4fi dist=%.4f\n', ...
                                j, k, real(z_opt), imag(z_opt), abs(z_opt-unpert_seed));
                        end
                    catch
                    end
                    k = k + 1;
                end
            end

            % ---- Coincident family duplication ----------------------------
            % When fsolve converges to the same point from both perturbed seeds,
            % we get one zero where two are expected.  Duplicate it.
            % Guard: only duplicate if still short of expected count.
            % For k>2 coincident families: replace multiplicity 2 with k.
            n_expected_asym = sum(floor((imag_max - B) ./ im_step));
            for j1 = 1:nPlus
                for j2 = j1+1:nPlus
                    if abs(re_asym(j1) - re_asym(j2)) < 1e-6
                        k1 = ceil(B / im_step(j1));
                        while true
                            h1 = k1 * im_step(j1);
                            if h1 > imag_max, break; end
                            k2 = round(h1 / im_step(j2));
                            h2 = k2 * im_step(j2);
                            if abs(h1 - h2) < min_step * 0.1
                                h_mid   = (h1 + h2) / 2;
                                half_w  = min_step * 0.4 * max(1, B / h_mid);
                                in_band = asym_zeros( ...
                                    imag(asym_zeros) >= h_mid - half_w & ...
                                    imag(asym_zeros) <= h_mid + half_w & ...
                                    abs(real(asym_zeros) - re_asym(j1)) < 1.0);
                                if numel(in_band) == 1 && ...
                                        numel(asym_zeros) < n_expected_asym + numel(rect_zeros)
                                    asym_zeros = [asym_zeros; in_band(1)]; %#ok<AGROW>
                                    fprintf('  Duplicated: %.4f%+.4fi (fam %d & %d Im=%.4f)\n', ...
                                        real(in_band(1)), imag(in_band(1)), j1, j2, h_mid);
                                elseif numel(in_band) == 0
                                    fprintf('  WARNING: no zero near Im=%.4f (fam %d & %d)\n', ...
                                        h_mid, j1, j2);
                                end
                            end
                            k1 = k1 + 1;
                        end
                    end
                end
            end

            % ---- Count check ----------------------------------------------
            genuinely_new = obj.removeIfNearExisting(asym_zeros, rect_zeros, tol_dedup);
            n_expected    = sum(floor((imag_max - B) ./ im_step));
            n_found       = numel(genuinely_new);
            fprintf('  Asymptotic zeros above B: found %d genuinely new, expected %d.\n', ...
                n_found, n_expected);
            if n_found < n_expected
                warning('PeriodicFluidModel:asymShortfall', ...
                    ['Asymptotic zero count short: found %d, expected %d above B=%.3f.\n' ...
                     'Consider increasing K_rect or checking for coincident families.'], ...
                    n_found, n_expected, B);
            end
        end

        % -----------------------------------------------------------------
        function candidates = extractLocalMinima(~, RE, IM, svd_grid)
        % EXTRACTLOCALMINIMA  Find grid points that are local minima of svd_grid.
        %
        % A point is a local minimum if it is strictly less than all 8
        % neighbours (or is on the boundary and less than all available
        % neighbours).  Returns a column vector of complex values.
            [nr, nc] = size(svd_grid);
            candidates = [];
            for ii = 1:nr
                for jj = 1:nc
                    val = svd_grid(ii,jj);
                    ii_nb = max(1,ii-1):min(nr,ii+1);
                    jj_nb = max(1,jj-1):min(nc,jj+1);
                    neighbourhood = svd_grid(ii_nb, jj_nb);
                    neighbourhood(ii - max(1,ii-1) + 1, jj - max(1,jj-1) + 1) = Inf;
                    if val < min(neighbourhood(:))
                        candidates = [candidates; RE(ii,jj) + 1i*IM(ii,jj)]; %#ok<AGROW>
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        function [coin_pairs, coin_heights] = findCoincidentFamilies(~, ...
                re_asym, im_step, imag_max, tol_re)
        % FINDCOINCIDENTFAMILIES  Detect family pairs sharing asymptotic real
        % part and compute heights where they produce near-coincident zeros.
        %
        % coin_pairs:   K x 2 matrix of family index pairs (j1, j2)
        % coin_heights: K x 1 cell array; coin_heights{k} lists imaginary
        %               heights in (0, imag_max] where families coin_pairs(k,:)
        %               both have seeds within min(im_step)/2 of each other.
            nPlus      = numel(re_asym);
            min_step   = min(im_step);
            coin_pairs   = zeros(0, 2);
            coin_heights = {};

            for j1 = 1:nPlus
                for j2 = j1+1:nPlus
                    if abs(re_asym(j1) - re_asym(j2)) < tol_re
                        % Coincident real parts: find coincident heights
                        heights = [];
                        k1 = 1;
                        while true
                            h1 = k1 * im_step(j1);
                            if h1 > imag_max, break; end
                            k2 = round(h1 / im_step(j2));
                            h2 = k2 * im_step(j2);
                            if k2 >= 1 && abs(h1 - h2) < min_step/2
                                heights = [heights; (h1+h2)/2]; %#ok<AGROW>
                            end
                            k1 = k1 + 1;
                        end
                        if ~isempty(heights)
                            coin_pairs(end+1,:)   = [j1, j2]; %#ok<AGROW>
                            coin_heights{end+1}   = heights;  %#ok<AGROW>
                            fprintf(['  Coincident families %d & %d ' ...
                                '(re=%.4f, %.4f): %d coincident heights.\n'], ...
                                j1, j2, re_asym(j1), re_asym(j2), numel(heights));
                        end
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        function all_z = deduplicateExcludingCoincident(~, zlist, tol_tight, ...
                coin_pairs, coin_heights, im_step)
        % DEDUPLICATEEXCLUDINGCOINCIDENT  Remove only numerical duplicates,
        % preserving genuinely distinct zeros at coincident heights.
        %
        % Two zeros are candidates for merging only if:
        %   (a) |z1 - z2| < tol_tight   (very close numerically)
        % AND
        %   (b) neither z1 nor z2 is at a known coincident height for any
        %       family pair — at coincident heights we expect two distinct
        %       zeros and never merge them.
            if isempty(zlist)
                all_z = [];
                return
            end

            % Build set of coincident heights (union over all pairs)
            coin_h_all = [];
            for k = 1:numel(coin_heights)
                coin_h_all = [coin_h_all; coin_heights{k}(:)]; %#ok<AGROW>
            end
            min_step = min(im_step);

            n   = numel(zlist);
            keep = true(n, 1);
            for i = 1:n
                if ~keep(i), continue; end
                for j = i+1:n
                    if ~keep(j), continue; end
                    if abs(zlist(i) - zlist(j)) < tol_tight
                        % Check if this imaginary height is a coincident height
                        h_mid = (imag(zlist(i)) + imag(zlist(j))) / 2;
                        at_coin = ~isempty(coin_h_all) && ...
                            any(abs(coin_h_all - h_mid) < min_step/2);
                        if ~at_coin
                            % Genuine duplicate: keep the one with smaller |det|
                            keep(j) = false;
                            fprintf(['  dedup: merged %.4f%+.4fi and ' ...
                                '%.4f%+.4fi (not at coincident height)\n'], ...
                                real(zlist(i)), imag(zlist(i)), ...
                                real(zlist(j)), imag(zlist(j)));
                        end
                    end
                end
            end
            all_z = zlist(keep);
        end

        % -----------------------------------------------------------------
        function checkFamilySpacing(obj, all_zeros, re_asym, im_step, ...
                coin_pairs, coin_heights, B, tol_re)
        % CHECKFAMILYSPACING  Per-family spacing check in asymptotic region.
        %
        % For each S_plus family j, extracts zeros above B whose real part
        % is closest to re_asym(j), checks for gaps > 1.5*im_step(j).
        % At coincident heights, two zeros are expected close together;
        % a single zero there is flagged as a possible miss.
            nPlus    = numel(re_asym);
            min_step = min(im_step);
            asym_z   = all_zeros(imag(all_zeros) > B);

            % Build coincident heights union
            coin_h_all = [];
            for k = 1:numel(coin_heights)
                coin_h_all = [coin_h_all; coin_heights{k}(:)]; %#ok<AGROW>
            end

            for j = 1:nPlus
                % Assign zeros to family j: nearest re_asym
                re_dists = abs(real(asym_z) - re_asym(j));
                % Use threshold: within 1 unit of re_asym(j) and closer to
                % re_asym(j) than to any other family
                in_fam = false(numel(asym_z), 1);
                for k = 1:numel(asym_z)
                    all_dists = abs(real(asym_z(k)) - re_asym);
                    [~, nearest] = min(all_dists);
                    if nearest == j && re_dists(k) < 1.0
                        in_fam(k) = true;
                    end
                end
                fam_z = asym_z(in_fam);
                if isempty(fam_z), continue; end

                fz = sort(imag(fam_z));
                if numel(fz) >= 2
                    gaps = diff(fz);
                    for gi = 1:numel(gaps)
                        if gaps(gi) > 1.5 * im_step(j)
                            gap_h = fz(gi) + gaps(gi)/2;
                            warning('PeriodicFluidModel:spacingGap', ...
                                ['Family %d (re~%.3f): gap of %.3f at ' ...
                                 'Im~%.3f (expected spacing %.3f). ' ...
                                 'Possible missing zero.'], ...
                                j, re_asym(j), gaps(gi), gap_h, im_step(j));
                        end
                    end
                end

                % Check coincident heights: expect 2 zeros within min_step/2
                if ~isempty(coin_h_all)
                    for ci = 1:numel(coin_h_all)
                        h = coin_h_all(ci);
                        if h <= B, continue; end
                        near_h = all_zeros(abs(imag(all_zeros) - h) < min_step/2);
                        if numel(near_h) < 2
                            warning('PeriodicFluidModel:coincidentMiss', ...
                                ['Coincident height Im~%.3f: expected 2 zeros, ' ...
                                 'found %d. Possible miss.'], h, numel(near_h));
                        end
                    end
                end
            end
        end


        function new_z = removeIfNearExisting(~, candidates, existing, tol)
        % REMOVEIFNEAREXISTING  Keep only candidates not within tol of any
        % existing zero.  existing zeros are kept unconditionally.
            if isempty(candidates) || isempty(existing)
                new_z = candidates(:);
                return
            end
            scale = max(abs([candidates(:); existing(:)])) + 1;
            keep  = true(numel(candidates), 1);
            for k = 1:numel(candidates)
                if any(abs(candidates(k) - existing(:)) < tol * scale)
                    keep(k) = false;
                end
            end
            new_z = candidates(keep);
        end


        function out = deduplicateZeros(~, zlist, tol)
        % DEDUPLICATEZEROS  Remove near-duplicate complex zeros.
        %
        % Converts to [real imag] pairs, applies uniquetol row-wise, then
        % reconstructs complex vector.
            if isempty(zlist)
                out = [];
                return
            end
            pairs = [real(zlist(:)), imag(zlist(:))];
            pairs = uniquetol(pairs, tol, 'ByRows', true, ...
                              'DataScale', max(abs(zlist(:)))+1);
            out   = pairs(:,1) + 1i*pairs(:,2);
        end

    end

    %% ----------------------------------------------------------------
    %  PRIVATE - Utilities
    %% ----------------------------------------------------------------
    methods (Access = private)

        function gamma_pm = gammaPhMajor(obj)
        % GAMMAPHMAYOR  gamma(t) in phase-major layout, 1 x (m*n).
            gamma_pm = zeros(1, obj.m * obj.n);
            for q = 1:obj.m
                gamma_pm((q-1)*obj.n + (1:obj.n)) = obj.gamma(:,q).';
            end
        end


        function checkComputed(obj, field)
            if isempty(obj.(field))
                error('PeriodicFluidModel:notComputed', ...
                    'Field "%s" not yet computed. Run the appropriate step first.', field);
            end
        end

        function sortedLam = sortEigenvalues(~, lam)
        % SORTEIGENVALUES  Sort zeros for truncation-safe ordering.
        %
        % Rules:
        %   1. Keep Im >= 0 only (conjugates handled via 2*real elsewhere).
        %   2. Real zeros first, sorted by descending Re (least negative
        %      first, since those decay slowest in x).
        %   3. Complex zeros sorted by ascending Im, ties broken by
        %      descending Re.  This ensures that truncating to any prefix
        %      always retains the least rapidly oscillating zeros across
        %      all families before taking higher rows from any family.
            lam      = lam(imag(lam) >= 0);
            is_real  = (abs(imag(lam)) == 0);

            % Real part: descending (least negative = slowest x-decay first)
            real_lam = sort(real(lam(is_real)), 'descend');

            % Complex part: ascending Im, descending Re for ties
            complex_lam = lam(~is_real);
            [~, idx]    = sortrows([imag(complex_lam(:)), -real(complex_lam(:))], [1, 2]);
            sortedLam   = [real_lam(:); complex_lam(idx)];
        end

        function Ts = buildTarray(obj)
        % Ts(:,:,j) = dt * T((j-0.5)/n).
            Ts = zeros(obj.m, obj.m, obj.n);
            for j = 1:obj.n
                Ts(:,:,j) = obj.dt * obj.T_fun((j-0.5)/obj.n);
            end
        end

        function formatPiAxis(~, imagRange, hax)
        % Label y-axis in multiples of pi.
            if nargin < 3, hax = gca; end
            ticks  = floor(imagRange(1)/pi) : ceil(imagRange(2)/pi);
            labels = cell(size(ticks));
            for i = 1:numel(ticks)
                k = ticks(i);
                if k==0,      labels{i} = '0';
                elseif k==1,  labels{i} = '$\pi$';
                elseif k==-1, labels{i} = '$-\pi$';
                else,         labels{i} = sprintf('$%d\\pi$', k);
                end
            end
            set(hax, 'YTick', ticks*pi, 'YTickLabel', labels, ...
                'TickLabelInterpreter', 'latex');
        end

    end

end
