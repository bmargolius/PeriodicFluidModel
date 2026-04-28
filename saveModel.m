function saveModel(model, example_name, save_dir)
% SAVEMODEL  Save computed results from a PeriodicFluidModel to a .mat file.
%
% The function handle T_fun cannot be stored in a mat file portably.
% Instead, the example_name string is stored so the model can be
% reconstructed via ExampleDefinitions and loadModel.
%
% INPUTS
%   model        - PeriodicFluidModel with at least eigenvalues computed
%   example_name - string identifier matching ExampleDefinitions
%   save_dir     - directory to save to.  Default: 'saved_models'
%
% SAVES
%   example_name, C_diag, n, imag_max, tol_dedup (if available),
%   eigenvalues, E (if computed), P0 (if computed), gamma_data

    if nargin < 3 || isempty(save_dir)
        save_dir = 'saved_models';
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    fname = fullfile(save_dir, sprintf('%s_n%d.mat', example_name, model.n));

    % Always save
    C_diag       = model.C_diag;
    n            = model.n;
    imag_max     = model.imag_max;
    K_rect       = model.K_rect;
    eigenvalues  = model.eigenvalues;

    % Optional computed quantities
    gamma_data   = model.gamma;         % n x m
    P0           = model.P0;            % 1 x (m*n)
    E            = model.E;             % nK x (m*n)

    save(fname, 'example_name', 'C_diag', 'n', 'imag_max', 'K_rect', ...
        'eigenvalues', 'gamma_data', 'P0', 'E');

    fprintf('Model saved to %s\n', fname);
    fprintf('  example:    %s\n', example_name);
    fprintf('  n:          %d\n', n);
    fprintf('  imag_max:   %.4g\n', imag_max);
    fprintf('  eigenvalues: %d\n', numel(eigenvalues));
    if ~isempty(E)
        fprintf('  E:          %d x %d\n', size(E,1), size(E,2));
    end
end
