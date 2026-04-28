function model = loadModel(example_name, save_dir, varargin)
% LOADMODEL  Load saved results and reconstruct a PeriodicFluidModel.
%
% INPUTS
%   example_name - string identifier matching ExampleDefinitions
%   save_dir     - directory to load from.  Default: 'saved_models'
%   varargin     - additional name-value pairs passed to PeriodicFluidModel
%
% RETURNS
%   model - PeriodicFluidModel with all saved quantities restored

    if nargin < 2 || isempty(save_dir)
        save_dir = 'saved_models';
    end

    % Find the file
    pattern = fullfile(save_dir, sprintf('%s_n*.mat', example_name));
    files   = dir(pattern);
    if isempty(files)
        error('loadModel: no saved file found for ''%s'' in %s', ...
            example_name, save_dir);
    end
    if numel(files) > 1
        [~, idx] = max([files.datenum]);
        fprintf('Multiple files found; loading most recent: %s\n', files(idx).name);
    else
        idx = 1;
    end
    fname = fullfile(save_dir, files(idx).name);
    fprintf('Loading from %s\n', fname);

    S = load(fname);

    % Reconstruct T_fun from ExampleDefinitions
    [T_fun, ~, description] = ExampleDefinitions(S.example_name);
    fprintf('Example: %s\n', description);

    % Construct model with saved parameters
    model = PeriodicFluidModel(T_fun, S.C_diag, ...
        'n',        S.n, ...
        'imag_max', S.imag_max, ...
        'K_rect',   S.K_rect, ...
        varargin{:});

    % Restore computed quantities via method (respects private SetAccess)
    model.restoreFromStruct(S);

    fprintf('Loaded: %d eigenvalues', numel(model.eigenvalues));
    if ~isempty(model.E)
        fprintf(', %d eigenfunctions', size(model.E,1));
    end
    fprintf('\n');
end
