function model = runExample(example_name, varargin)
% RUNEXAMPLE  Run a named example end-to-end and optionally save results.
%
% Constructs a PeriodicFluidModel from ExampleDefinitions, runs the full
% pipeline (gamma, zeros, P0, eigenfunctions), and optionally saves.
%
% USAGE
%   model = runExample('model2')
%   model = runExample('model2', 'n', 400, 'imag_max', 20*pi)
%   model = runExample('model2', 'save', true, 'save_dir', 'my_results')
%
% NAME-VALUE OPTIONS (in addition to PeriodicFluidModel options)
%   'save'      - save results after computation.  Default: false
%   'save_dir'  - directory for saving.  Default: 'saved_models'
%   'load'      - try to load from file first; compute only if not found.
%                 Default: false

    % Parse our options out before passing rest to constructor
    p = inputParser();
    addParameter(p, 'save',     false, @islogical);
    addParameter(p, 'save_dir', 'saved_models', @ischar);
    addParameter(p, 'load',     false, @islogical);
    % PeriodicFluidModel options
    addParameter(p, 'n',         100,    @(x) isscalar(x) && x>0);
    addParameter(p, 'imag_max',  20*pi,  @(x) isscalar(x) && x>0);
    addParameter(p, 'K_rect',    3,      @(x) isscalar(x) && x>0);
    addParameter(p, 'WfMethod',  'expm', @(x) ismember(x,{'expm','ode'}));
    addParameter(p, 'ResMethod', 'adjugate', @(x) ismember(x,{'adjugate','eigvec'}));
    parse(p, varargin{:});

    do_save  = p.Results.save;
    save_dir = p.Results.save_dir;
    do_load  = p.Results.load;

    % Try loading first if requested
    if do_load
        try
            model = loadModel(example_name, save_dir);
            fprintf('Loaded from file. Skipping computation.\n');
            return
        catch
            fprintf('No saved file found. Computing from scratch.\n');
        end
    end

    % Get example definition
    [T_fun, C_diag, description, notes] = ExampleDefinitions(example_name);
    fprintf('\n=== %s ===\n', description);
    if ~isempty(notes)
        fprintf('Notes: %s\n', notes);
    end

    % Construct model
    model = PeriodicFluidModel(T_fun, C_diag, ...
        'n',         p.Results.n, ...
        'imag_max',  p.Results.imag_max, ...
        'K_rect',    p.Results.K_rect, ...
        'WfMethod',  p.Results.WfMethod, ...
        'ResMethod', p.Results.ResMethod);

    % Run pipeline
    fprintf('\n--- Step 1: computeGamma ---\n');
    model.computeGamma();

    fprintf('\n--- Step 2: findZeros ---\n');
    model.findZeros();

    fprintf('\n--- Step 3: computeP0 ---\n');
    model.computeP0();

    fprintf('\n--- Step 4: computeEigenfunctions ---\n');
    model.computeEigenfunctions();

    fprintf('\n=== Pipeline complete ===\n');
    model.disp();

    % Save if requested
    if do_save
        saveModel(model, example_name, save_dir);
    end
end
