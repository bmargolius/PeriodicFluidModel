function [T_fun, C_diag, description, notes] = ExampleDefinitions(example_name)
% EXAMPLEDEFINITIONS  Returns T_fun and C_diag for a named example.
%
% USAGE
%   [T_fun, C_diag, description, notes] = ExampleDefinitions('model2')
%   ExampleDefinitions()   % list all available examples
%
% All examples have C_diag in canonical order [+...+ -...- 0...0].
% Unstable examples are caught by computeGamma (mu >= 0 error).
%

    if nargin < 1
        fprintf('Available examples:\n');
        names = {'model1','model2','const2', ...
                 'ex1','ex2','ex3','ex4','ex5','ex6','ex7','ex8','ex9', ...
                 'ex10','ex11','ex12','ex13','ex14','ex15','ex16','ex17','ex18'};
        for k = 1:numel(names)
            try
                [~,~,desc] = ExampleDefinitions(names{k});
                fprintf('  %-12s  %s\n', names{k}, desc);
            catch ME
                fprintf('  %-12s  ERROR: %s\n', names{k}, ME.message);
            end
        end
        T_fun = []; C_diag = []; description = ''; notes = '';
        return
    end

    a  = @(t) (1 + 0.8*sin(2*pi*t));
    b  = @(t) 0.8*(1 - 0.8*sin(2*pi*t));
    z  = @(t) 0*t;
    f  = @(t) 1 + z(t);
    hh = @(t) heaviside(t - floor(t) - 0.5);
    notes = '';

    switch lower(example_name)

        case 'model1'
            description = 'Paper Model 1: piecewise-constant T, C=[4 2 1 -12 0 0]';
            notes = 'Coincident asymptotic families at Re=-6.';
            C_diag = [4 2 1 -12 0 0];
            T1 = @(t) [-10*f(t)  2*f(t)   2*f(t)   6*f(t)     z(t)        z(t);
                          f(t)  -9*f(t)    f(t)     3*f(t)    3*f(t)      f(t);
                          z(t)   z(t)     -7*f(t)   4*f(t)    2*f(t)      f(t);
                         0.8*f(t) 0.8*f(t) 0.8*f(t) -8*f(t)   4.8*f(t)  0.8*f(t);
                         0.8*f(t) 0.8*f(t) 2.6*f(t) 5.2*f(t) -10.2*f(t) 0.8*f(t);
                          f(t)   2*f(t)    f(t)     5.2*f(t)  6*f(t)    -15.2*f(t)];
            T2 = @(t) [-2*f(t)    z(t)    z(t)    2*f(t)    z(t)       z(t);
                         f(t)   -15*f(t)  f(t)    5*f(t)   7*f(t)     f(t);
                        0.4*f(t) 0.3*f(t) -2*f(t) 0.4*f(t) 0.8*f(t)  0.1*f(t);
                        0.8*f(t) 0.8*f(t) 0.8*f(t) -4.8*f(t) 1.6*f(t) 0.8*f(t);
                        0.8*f(t) 0.8*f(t) 0.6*f(t)  1.2*f(t) -4.2*f(t) 0.8*f(t);
                         f(t)    z(t)     f(t)    1.2*f(t)  2*f(t)   -5.2*f(t)];
            T_fun = @(t) (1-hh(t))*T1(t) + hh(t)*T2(t);

        case 'model2'
            description = 'Paper Model 2: sinusoidal T, C=[1 2 3 -4 -3 0]';
            C_diag = [1 2 3 -4 -3 0];
            T_fun = @(t) [-3*a(t)  a(t)    a(t)    a(t)    z(t)    z(t);
                           a(t)   -5*a(t)  a(t)    a(t)    a(t)    a(t);
                           a(t)    a(t)   -5*a(t)  a(t)    a(t)    a(t);
                           b(t)    b(t)    b(t)   -5*b(t)  b(t)    b(t);
                           b(t)    b(t)    b(t)    b(t)   -5*b(t)  b(t);
                           a(t)    a(t)    a(t)    b(t)    5*b(t) -3*a(t)-6*b(t)];

        case 'const2'
            description = 'Constant-rate version of model2, C=[1 2 3 -4 -3 0]';
            notes = 'Constant T: verify against standard fluid model.';
            C_diag = [1 2 3 -4 -3 0];
            T_fun = @(t) [-3    1    1    1    z(t)  z(t);
                           1   -5    1    1    1     1;
                           1    1   -5    1    1     1;
                           0.8  0.8  0.8 -4    0.8   0.8;
                           0.8  0.8  0.8  0.8 -4     0.8;
                           1    1    1    0.8   4    -3-6*0.8];

        case 'ex1'
            description = 'ex1: sinusoidal T, C=[1 2 3 -40 -30 0]';
            C_diag = [1 2 3 -40 -30 0];
            T_fun = @(t) [-6*a(t)   a(t)    a(t)   4*a(t)    z(t)          z(t);
                           a(t)   -12*a(t)  a(t)   4*a(t)   5*a(t)         a(t);
                           a(t)    a(t)   -12*a(t) 4*a(t)   5*a(t)         a(t);
                           b(t)    b(t)    b(t)   -8*b(t)   4*b(t)         b(t);
                           b(t)    b(t)   2*b(t)  4*b(t)   -9*b(t)         b(t);
                           a(t)    a(t)    a(t)   4*b(t)   5*b(t)   -3*a(t)-9*b(t)];

        case 'ex2'
            description = 'ex2: constant T, C=[1 2 3 -40 -30 0]';
            notes = 'Constant T.';
            C_diag = [1 2 3 -40 -30 0];
            T_fun = @(t) f(t)*[-6   1   1   4   0   0;
                                 1 -12   1   4   5   1;
                                 1   1 -12   4   5   1;
                                0.8 0.8 0.8 -6.4 3.2 0.8;
                                0.8 0.8 1.6  3.2 -7.2 0.8;
                                1   1   1   3.2  4  -10.2];

        case 'ex3'
            description = 'ex3: asymmetric sinusoidal, C=[2 3 4 -10 -5 0]';
            C_diag = [2 3 4 -10 -5 0];
            T_fun = @(t) [-6*a(t)             a(t)              a(t)   4*a(t)    z(t)          z(t);
                          0.2*a(t)  -1.2*a(t)-b(t)  0.2*a(t)   b(t)  0.2*a(t)  0.6*a(t);
                          2*a(t)              a(t)            -5*a(t) 0.5*a(t)   a(t)         0.5*a(t);
                          b(t)                b(t)             b(t)  -8*b(t)    4*b(t)          b(t);
                          b(t)                b(t)            2*b(t)  4*b(t)   -9*b(t)          b(t);
                          a(t)                a(t)             a(t)   4*b(t)   5*b(t)  -3*a(t)-9*b(t)];

        case 'ex4'
            description = 'ex4: sinusoidal, C=[4 1 1 -8 -12 0]';
            C_diag = [4 1 1 -8 -12 0];
            T_fun = @(t) [-2*f(t)-4*a(t)  f(t)    a(t)    3*a(t)   f(t)                     z(t);
                           b(t)          -5*b(t)  2*b(t)   z(t)    b(t)                      b(t);
                           5/4*b(t)       z(t)  -5/4*b(t)  z(t)    z(t)                      z(t);
                           f(t)           f(t)    a(t)   -3*f(t)-a(t)  z(t)                  f(t);
                           f(t)           f(t)    a(t)    b(t)   -2*f(t)-2*a(t)-b(t)         a(t);
                           b(t)           a(t)    f(t)    b(t)    f(t)          -2*b(t)-2*f(t)-a(t)];

        case 'ex5'
            description = 'ex5: sinusoidal, C=[2 2 1 -5 -5 0]';
            C_diag = [2 2 1 -5 -5 0];
            T_fun = @(t) [-f(t)-4*a(t)   z(t)    a(t)    3*a(t)   f(t)                     z(t);
                           z(t)         -3*a(t)   z(t)    a(t)    2*a(t)                    z(t);
                           5/4*b(t)      z(t)  -5/4*b(t)  z(t)    z(t)                      z(t);
                           f(t)          f(t)    a(t)   -3*f(t)-a(t)  z(t)                  f(t);
                           f(t)          f(t)    a(t)    b(t)   -2*f(t)-2*a(t)-b(t)         a(t);
                           b(t)          a(t)    f(t)    b(t)    f(t)          -2*b(t)-2*f(t)-a(t)];

        case 'ex6'
            description = 'ex6: sinusoidal T, C=[1 2 3 -4 -3 0]';
            notes = 'Same T structure as model2.';
            C_diag = [1 2 3 -4 -3 0];
            T_fun = @(t) [-3*a(t)  a(t)    a(t)    a(t)    z(t)    z(t);
                           a(t)   -5*a(t)  a(t)    a(t)    a(t)    a(t);
                           a(t)    a(t)   -5*a(t)  a(t)    a(t)    a(t);
                           b(t)    b(t)    b(t)   -5*b(t)  b(t)    b(t);
                           b(t)    b(t)    b(t)    b(t)   -5*b(t)  b(t);
                           a(t)    a(t)    a(t)    b(t)    5*b(t) -3*a(t)-6*b(t)];

        case 'ex7'
            description = 'ex7: sinusoidal variant, C=[1 2 3 -40 -30 0]';
            C_diag = [1 2 3 -40 -30 0];
            T_fun = @(t) [-6*a(t)   a(t)    a(t)   4*a(t)    z(t)          z(t);
                           a(t)   -12*a(t)  a(t)   4*a(t)   5*a(t)         a(t);
                           a(t)    a(t)   -12*a(t) 4*a(t)   5*a(t)         a(t);
                           b(t)    b(t)    b(t)   -9*b(t)   5*b(t)         b(t);
                           b(t)    b(t)    b(t)    4*b(t)   -8*b(t)        b(t);
                           a(t)    a(t)    a(t)   4*b(t)   5*b(t)   -3*a(t)-9*b(t)];

        case 'ex8'
            description = 'ex8: no zero phase, C=[1 2 3 -40 -30]';
            notes = 'No zero-rate phases.';
            C_diag = [1 2 3 -40 -30];
            T_fun = @(t) [-3*a(t)  a(t)    a(t)    a(t)    z(t);
                           a(t)   -4*a(t)  a(t)    a(t)    a(t);
                           a(t)    a(t)   -4*a(t)  a(t)    a(t);
                           b(t)    b(t)    b(t)   -4*b(t)  b(t);
                           b(t)    b(t)    b(t)    b(t)   -4*b(t)];

        case 'ex9'
            description = 'ex9: constant T, no zero phase, C=[1 2 3 -40 -30]';
            notes = 'Constant T, no zero-rate phases.';
            C_diag = [1 2 3 -40 -30];
            T_fun = @(t) f(t)*[-3   1   1   1   0;
                                 1  -4   1   1   1;
                                 1   1  -4   1   1;
                                0.8 0.8 0.8 -3.2 0.8;
                                0.8 0.8 0.8  0.8 -3.2];

        case 'ex10'
            description = 'ex10: small rates, 5 phases, C=[1 1 -1 -1 0]';
            notes = 'Rates scaled by 0.01.';
            C_diag = [1 1 -1 -1 0];
            aa = 1.1;
            T_fun = @(t) 0.01*[...
               -2*f(t)-2.2*f(t)          f(t)    aa*f(t)                  aa*f(t)                    f(t);
                f(t)  -4.2*f(t)-0.5*sin(2*pi*t)  aa*f(t)+0.5*sin(2*pi*t)  aa*f(t)                   f(t);
                f(t)   f(t)            -4.1*f(t)-0.5*sin(2*pi*t)  aa*f(t)+0.5*sin(2*pi*t)            f(t);
                f(t)   f(t)             aa*f(t)                  -4.1*f(t)                            f(t);
                f(t)   f(t)             aa*f(t)           aa*f(t)+0.5*sin(2*pi*t)  -4.2*f(t)-0.5*sin(2*pi*t)];

        case 'ex11'
            description = 'ex11: sinusoidal, C=[1 2 -3 -4 -3 0]';
            C_diag = [1 2 -3 -4 -3 0];
            T_fun = @(t) [-6*a(t)   a(t)    a(t)   4*a(t)    z(t)          z(t);
                           a(t)   -12*a(t)  a(t)   4*a(t)   5*a(t)         a(t);
                           a(t)    a(t)   -12*a(t) 4*a(t)   5*a(t)         a(t);
                           b(t)    b(t)    b(t)   -8*b(t)   4*b(t)         b(t);
                           b(t)    b(t)   2*b(t)  4*b(t)   -9*b(t)         b(t);
                           a(t)    a(t)    a(t)   4*b(t)   5*b(t)   -3*a(t)-9*b(t)];

        case 'ex12'
            description = 'ex12: sinusoidal, C=[1 2 3 -4 0 0]';
            notes = 'Two zero-rate phases.';
            C_diag = [1 2 3 -4 0 0];
            T_fun = @(t) [-6*a(t)   a(t)    a(t)   4*a(t)    z(t)          z(t);
                           a(t)   -12*a(t)  a(t)   4*a(t)   5*a(t)         a(t);
                           a(t)    a(t)   -12*a(t) 4*a(t)   5*a(t)         a(t);
                           b(t)    b(t)    b(t)   -8*b(t)   4*b(t)         b(t);
                           b(t)    b(t)   2*b(t)  4*b(t)   -9*b(t)         b(t);
                           a(t)    a(t)    a(t)   4*b(t)   5*b(t)   -3*a(t)-9*b(t)];

        case 'ex13'
            description = 'ex13: piecewise-constant T variant, C=[1 2 3 -40 -30 0]';
            notes = 'Piecewise-constant switching at t=0.5.';
            C_diag = [1 2 3 -40 -30 0];
            T1 = @(t) [-10*f(t)  2*f(t)   2*f(t)   6*f(t)     z(t)        z(t);
                          f(t)  -9*f(t)    f(t)     3*f(t)    3*f(t)      f(t);
                          z(t)   z(t)     -7*f(t)   4*f(t)    2*f(t)      f(t);
                         0.8*f(t) 0.8*f(t) 0.8*f(t) -8*f(t)   4.8*f(t)  0.8*f(t);
                         0.8*f(t) 0.8*f(t) 2.6*f(t) 5.2*f(t) -10.2*f(t) 0.8*f(t);
                          f(t)   2*f(t)    f(t)     5.2*f(t)  6*f(t)    -15.2*f(t)];
            T2 = @(t) [-2*f(t)    z(t)    z(t)    2*f(t)    z(t)       z(t);
                         f(t)   -15*f(t)  f(t)    5*f(t)   7*f(t)     f(t);
                        2*f(t)   2*f(t) -17*f(t)  4*f(t)   8*f(t)     f(t);
                        0.8*f(t) 0.8*f(t) 0.8*f(t) -4.8*f(t) 1.6*f(t) 0.8*f(t);
                        0.8*f(t) 0.8*f(t) 0.6*f(t)  1.2*f(t) -4.2*f(t) 0.8*f(t);
                         f(t)    z(t)     f(t)    1.2*f(t)  2*f(t)   -5.2*f(t)];
            T_fun = @(t) (1-hh(t))*T1(t) + hh(t)*T2(t);

        case 'ex14'
            description = 'ex14: sinusoidal, C=[1 2 3 -40 -30 0]';
            C_diag = [1 2 3 -40 -30 0];
            T_fun = @(t) [-6*a(t)   z(t)    z(t)    a(t)    a(t)   4*a(t);
                           4*a(t) -12*a(t)  a(t)    a(t)   5*a(t)  a(t);
                           a(t)    a(t)   -12*a(t)  4*a(t)  a(t)   5*a(t);
                           b(t)   2*b(t)   2*b(t)  -8*b(t)  2*b(t) b(t);
                           b(t)   4*b(t)   2*b(t)   b(t)   -9*b(t) b(t);
                           a(t)    a(t)    a(t)    4*b(t)  5*b(t) -3*a(t)-9*b(t)];

        case 'ex15'
            description = 'ex15: 5 phases, C=[2 1 -5 -4 0]';
            C_diag = [2 1 -5 -4 0];
            a2 = @(t) (1 + 0.8*sin(2*pi*t));
            b2 = @(t) 0.8*(1 - 0.8*sin(2*pi*t));
            T_fun = @(t) [-2*a2(t)    z(t)      z(t)     a2(t)    a2(t);
                           4*a2(t) -11*a2(t)   a2(t)    a2(t)   5*a2(t);
                           a2(t)    a2(t)    -7*a2(t)  4*a2(t)  a2(t);
                           3*b2(t)   z(t)      b2(t)   -6*b2(t) 2*b2(t);
                           2*b2(t)  4*b2(t)   b2(t)    b2(t)   -8*b2(t)];

        case 'ex16'
            description = 'ex16: 5 phases, C=[2 1 -8 -4 0]';
            C_diag = [2 1 -8 -4 0];
            a2 = @(t) (1 + 0.8*sin(2*pi*t));
            b2 = @(t) 0.8*(1 - 0.8*sin(2*pi*t));
            T_fun = @(t) [-2*a2(t)    z(t)      z(t)     a2(t)    a2(t);
                           4*a2(t) -11*a2(t)   a2(t)    a2(t)   5*a2(t);
                           a2(t)    a2(t)    -7*a2(t)  4*a2(t)  a2(t);
                           3*b2(t)   z(t)      b2(t)   -6*b2(t) 2*b2(t);
                           2*b2(t)  4*b2(t)   b2(t)    b2(t)   -8*b2(t)];

        case 'ex17'
            description = 'ex17: sinusoidal, C=[1 2 3 -4 -5 0]';
            C_diag = [1 2 3 -4 -5 0];
            T_fun = @(t) [-6*a(t)   a(t)    a(t)   4*a(t)    z(t)          z(t);
                           a(t)   -12*a(t)  a(t)   4*a(t)   5*a(t)         a(t);
                           a(t)    a(t)   -12*a(t) 4*a(t)   5*a(t)         a(t);
                           b(t)    b(t)    b(t)   -8*b(t)   4*b(t)         b(t);
                           b(t)    b(t)   2*b(t)  4*b(t)   -9*b(t)         b(t);
                           a(t)    a(t)    a(t)   4*b(t)   5*b(t)   -3*a(t)-9*b(t)];

        case 'ex18'
            description = 'ex18: sinusoidal 6-phase, C=[4 2 1 -6 -5 0]';
            C_diag = [4 2 1 -6 -5 0];
            a2 = @(t) (1 + 0.8*sin(2*pi*t));
            b2 = @(t) 0.8*(1 - 0.8*sin(2*pi*t));
            T_fun = @(t) [-2*a2(t)-b2(t)        z(t)          z(t)      a2(t)    a2(t)           b2(t);
                           4*a2(t)   -11*a2(t)-2*b2(t)    a2(t)     a2(t)   2*b2(t)         5*a2(t);
                           a2(t)      a2(t)        -7*a2(t)   4*a2(t)  a2(t)      z(t);
                           3*b2(t)    z(t)           b2(t)    -6*b2(t)  z(t)      2*b2(t);
                           2*b2(t)   4*b2(t)         b2(t)     b2(t)   -8*b2(t)-a2(t)  a2(t);
                           3*a2(t)   5*b2(t)        2*a2(t)   a2(t)+b2(t) 6*b2(t) -6*a2(t)-12*b2(t)];

        otherwise
            error('ExampleDefinitions: unknown example ''%s''. Call with no args to list.', ...
                example_name);
    end
end
