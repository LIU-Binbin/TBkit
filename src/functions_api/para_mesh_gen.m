function mesh_list = para_mesh_gen(AtoB, np_list, opts)
    % Generate a mesh in an arbitrary-dimensional parameter space
    % AtoB: Range of parameters [min, max] for each parameter.
    % np_list: Number of points for each dimension or a single number.
    % opts.shift: Type of shift ('none', 'uni-rand', 'tot-rand').
    
    % Validate input arguments
    arguments
        AtoB (:, 2) double = [0, 1; 0, 1; 0, 1]  % Default range for 3 parameters
        np_list double = [10, 10, 10]  % Default resolution for 3 parameters
        opts.shift {mustBeMember(opts.shift, {'none', 'uni-rand', 'tot-rand'})} = 'none'; % Default shift type
    end

    nparas = length(np_list);  % Number of parameters
    vecs = cell(1, nparas);  % Cell array to store the vectors

    % Generate linearly spaced vectors for each parameter
    for i = 1:nparas
        vecs{i} = linspace(AtoB(i, 1), AtoB(i, 2), np_list(i));
    end

    % Determine the step size for each parameter (used for shifts)
    step = (AtoB(:, 2) - AtoB(:, 1)) ./ np_list';

    % Apply shifts to the vectors if specified
    switch opts.shift
        case 'none'
            % No shift applied
        case 'uni-rand'
            % Uniform random shift for each parameter
            for i = 1:nparas
                vecs{i} = vecs{i} + 1e-2 * step(i) * (1 + rand());
            end
        case 'tot-rand'
            % Total random shift for each parameter
            for i = 1:nparas
                vecs{i} = vecs{i} + 1e-2 * step(i) * (1 + rand(1, np_list(i)));
            end
    end

    % Create meshgrid based on the number of parameters (supports up to 4 parameters)
    switch nparas
        case 2
            [p1, p2] = ndgrid(vecs{1}, vecs{2});
            mesh_list = [p1(:), p2(:)];
        case 3
            [p1, p2, p3] = ndgrid(vecs{1}, vecs{2}, vecs{3});
            mesh_list = [p1(:), p2(:), p3(:)];
        case 4
            [p1, p2, p3, p4] = ndgrid(vecs{1}, vecs{2}, vecs{3}, vecs{4});
            mesh_list = [p1(:), p2(:), p3(:), p4(:)];
        otherwise
            error("Para_dim=%d not supported. Please provide between 2 and 4 parameters.", nparas);
    end
end
