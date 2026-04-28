function plotAllForModel(example_name, varargin)
% PLOTALLFORMODEL  Generate all diagnostic plots for a named example.
%
% Loads saved model and simulation results, generates all plots, and
% optionally saves figures to a per-model directory.
%
% USAGE
%   plotAllForModel('model2')
%   plotAllForModel('model2', 'save_figs', true, 'close_figs', true)
%
% OPTIONS
%   save_dir   - where to load model from.  Default: 'saved_models'
%   fig_dir    - base directory for figures.  Default: 'figures'
%   save_figs  - save figures to fig_dir/example_name/.  Default: false
%   close_figs - close each figure after saving.  Default: true if
%                save_figs is true, false otherwise.  Set explicitly
%                to keep figures open even when saving, or close without
%                saving (e.g. when running a batch).
%   x_levels   - x thresholds for CDF plots.  Default: [0.2 0.5 1.0 2.0 5.0]
%   x_max      - upper fluid level for density plots.  Default: 3.0
%   n_x        - number of x points for density.  Default: 100
%   n_events   - events per simulation run.  Default: 50000
%   n_runs     - simulation runs.  Default: 20
%   max_fluid  - max fluid for simulation bins.  Default: 5
%   sim_mesh   - simulation time mesh.  Default: 100

    p = inputParser();
    addParameter(p, 'save_dir',   'saved_models',          @ischar);
    addParameter(p, 'fig_dir',    'figures',               @ischar);
    addParameter(p, 'save_figs',  false,                   @islogical);
    addParameter(p, 'close_figs', [],                      @(x) islogical(x) || isempty(x));
    addParameter(p, 'x_levels',   [0.2 0.5 1.0 2.0 5.0], @isnumeric);
    addParameter(p, 'x_max',      3.0,                    @isnumeric);
    addParameter(p, 'n_x',        100,                    @isnumeric);
    addParameter(p, 'n_events',   50000,                  @isnumeric);
    addParameter(p, 'n_runs',     20,                     @isnumeric);
    addParameter(p, 'max_fluid',  5,                      @isnumeric);
    addParameter(p, 'sim_mesh',   100,                    @isnumeric);
    parse(p, varargin{:});
    opt = p.Results;

    % Default: close figures if saving, keep open if just viewing
    if isempty(opt.close_figs)
        opt.close_figs = opt.save_figs;
    end

    tic;

    % ---- Load model -------------------------------------------------------
    fprintf('\n=== plotAllForModel: %s ===\n', example_name);
    model = loadModel(example_name, opt.save_dir);

    % ---- Set up figure directory ------------------------------------------
    if opt.save_figs
        fig_subdir = fullfile(opt.fig_dir, example_name);
        if ~exist(fig_subdir, 'dir')
            mkdir(fig_subdir);
        end
        fprintf('Figures will be saved to %s\n', fig_subdir);
    end

    function save_and_maybe_close(name)
        if opt.save_figs
            fname = fullfile(fig_subdir, sprintf('%s_%s.png', example_name, name));
            exportgraphics(gcf, fname, 'Resolution', 150);
            fprintf('  Saved: %s\n', fname);
        end
        if opt.close_figs
            close(gcf);
        end
    end

    % ---- Gamma ------------------------------------------------------------
    figure;
    model.plotGamma();
    save_and_maybe_close('gamma');

    % ---- Drift ------------------------------------------------------------
    figure;
    model.plotDrift();
    save_and_maybe_close('drift');

    % ---- P0 ---------------------------------------------------------------
    figure;
    model.plotP0();
    save_and_maybe_close('P0');

    % ---- Mean fluid -------------------------------------------------------
    figure;
    model.plotMeanFluid();
    save_and_maybe_close('mean_fluid');

    % ---- Phase portrait with eigenvalues ----------------------------------
    re_min = floor(min(real(model.eigenvalues))) - 1;
    figure;
    model.plotPhasePlot([re_min 0.5], [0 model.imag_max], [300 600]);
    hold on;
    plot(model.eigenvalues, 'w.', 'MarkerSize', 6);
    hold off;
    save_and_maybe_close('phase_portrait');

    % ---- Simulation: load or generate -------------------------------------
    sim_file = fullfile(opt.save_dir, sprintf('%s_sim.mat', example_name));
    if exist(sim_file, 'file')
        fprintf('\nLoading simulation from %s\n', sim_file);
        load(sim_file, 'pdistcum');
    else
        fprintf('\nNo simulation file found. Running simulation...\n');
        pdistcum = model.simulate(opt.n_events, opt.n_runs, ...
                                  opt.max_fluid, opt.sim_mesh);
        save(sim_file, 'pdistcum');
        fprintf('Simulation saved to %s\n', sim_file);
    end

    % ---- CDF with simulation overlay, one figure per phase ---------------
    for j = 1:model.m
        figure;
        model.plotCDFwithSim(pdistcum, j, opt.x_levels, opt.sim_mesh);
        save_and_maybe_close(sprintf('cdf_phase%d', j));
    end

    % ---- Density contourf, one figure per phase --------------------------
    x_vals = linspace(0, opt.x_max, opt.n_x).';
    for j = 1:model.m
        figure;
        model.plotDensity(x_vals, j);
        save_and_maybe_close(sprintf('density_phase%d', j));
    end

    elapsed = toc;
    fprintf('\nplotAllForModel: %s complete in %.1f seconds.\n', ...
        example_name, elapsed);

    if exist('ringBell', 'file') || exist('ringBell', 'builtin')
        ringBell;
    end
end