function accumulator = kloop_accumulate(klist, kfun, accumulator, optionsParallel, optionsAdapt)
%KLOOP_ACCUMULATE Adaptively accumulate k-dependent contributions.
%
%   This helper evaluates the per-k-point kernel once for every k-point,
%   ranks their relative magnitudes, and keeps only the dominant ones whose
%   weight is comparable to (or larger than) the average 1/Nk contribution.
%   The surviving subset is used for the final summation, realizing an
%   "adapted-k" integration strategy without duplicating boilerplate code.
%
%   KLIST             : [Nk x 3] matrix of k-points.
%   KFUN              : function handle that receives one k-point row and
%                       returns a contribution with the same shape as the
%                       accumulator.
%   ACCUMULATOR       : initial value of the accumulated quantity.
%   OPTIONSPARALLEL   : structure with field `ncore` describing how many
%                       workers should be used (>=1).
%   OPTIONSADAPT      : (optional) structure controlling the adaptive logic
%       .enable            : set false to recover the legacy behaviour.
%       .importance_factor : keep k-points whose weight exceeds
%                            importance_factor * (1 / Nk) of the total.
%       .min_keep          : minimum number of k-points retained.
%       .reuse_samples     : reuse the coarse-evaluated contributions when
%                            true (default), otherwise re-evaluate them.
%       .mesh_shape        : [Nk1 Nk2 Nk3] describing the original tensor
%                            product grid used to build klist.
%       .refine_factor     : integer refinement factor applied to each
%                            active direction when generating the
%                            dense local meshes (e.g. 6×6×6).
arguments
    klist double
    kfun (1,1) function_handle
    accumulator
    optionsParallel.ncore (1,1) double = 1
    optionsAdapt.enable (1,1) logical = true
    optionsAdapt.importance_factor (1,1) double {mustBePositive} = 0.2
    optionsAdapt.min_keep (1,1) double {mustBePositive} = 32
    optionsAdapt.reuse_samples (1,1) logical = true
    optionsAdapt.mesh_shape double = []
    optionsAdapt.refine_factor (1,1) double {mustBePositive} = 1
end

Nk = size(klist, 1);
if Nk == 0
    return;
end

use_parallel = (optionsParallel.ncore > 1);

if ~optionsAdapt.enable
    accumulator = accumulate_all(klist, kfun, accumulator, use_parallel);
    return;
end

% First pass: evaluate every k-point once and measure its contribution.
[contrib_samples, weight_matrix] = evaluate_weighted_samples( ...
    klist, kfun, use_parallel);

total_weights_per_entry = sum(weight_matrix, 1);
if all(total_weights_per_entry == 0)
    % Degenerate case: every contribution vanished, so just accumulate all
    % coarse samples and return.
    accumulator = accumulate_from_samples(contrib_samples, accumulator);
    return;
end

avg_weights = total_weights_per_entry / Nk;
thresholds = avg_weights * optionsAdapt.importance_factor;

% Compare each matrix entry independently and keep the union of the
% positions whose contribution exceeds the per-entry threshold.
active_entries = total_weights_per_entry > 0;
selected = false(Nk, 1);
if any(active_entries)
    entry_weights = weight_matrix(:, active_entries);
    entry_thresholds = thresholds(active_entries);
    selection_matrix = bsxfun(@ge, entry_weights, entry_thresholds);
    selected = any(selection_matrix, 2);
end

% Fall back to the total Frobenius weight for ordering when enforcing
% min_keep.
aggregate_weights = sum(weight_matrix, 2);

min_keep = min(Nk, max(1, round(optionsAdapt.min_keep)));
if nnz(selected) < min_keep
    [~, order] = sort(aggregate_weights, 'descend'); %#ok<UDIM>
    selected(order(1:min_keep)) = true;
end

selected_idx = find(selected);
unselected_idx = find(~selected);

% Always keep the coarse contributions from the suppressed k-points so that
% the integral remains faithful to the original sum.
if ~isempty(unselected_idx)
    accumulator = accumulate_from_samples(contrib_samples, accumulator, unselected_idx);
end

if isempty(selected_idx)
    return;
end

refine_factor = max(1, round(optionsAdapt.refine_factor));

if refine_factor == 1
    if optionsAdapt.reuse_samples
        accumulator = accumulate_from_samples(contrib_samples, accumulator, selected_idx);
    else
        subset_samples = evaluate_samples(klist(selected_idx,:), kfun, use_parallel);
        accumulator = accumulate_from_samples(subset_samples, accumulator);
    end
    return;
end

Nk = size(klist, 1);
mesh_shape = validate_mesh_shape(optionsAdapt.mesh_shape, Nk);
refine_offsets = build_refine_offsets(klist, mesh_shape, refine_factor);
refined_samples = evaluate_refined_samples(klist(selected_idx,:), kfun, refine_offsets, use_parallel);
accumulator = accumulate_from_samples(refined_samples, accumulator);
end

function accumulator = accumulate_all(klist, kfun, accumulator, use_parallel)
samples = evaluate_samples(klist, kfun, use_parallel);
accumulator = accumulate_from_samples(samples, accumulator);
end

function samples = evaluate_samples(klist, kfun, use_parallel)
Nk = size(klist, 1);
samples = cell(Nk, 1);
if use_parallel
    parfor ki = 1:Nk %#ok<PFBNS>
        samples{ki} = kfun(klist(ki,:));
    end
else
    for ki = 1:Nk
        samples{ki} = kfun(klist(ki,:));
    end
end
end

function [samples, weight_matrix] = evaluate_weighted_samples(klist, kfun, use_parallel)
Nk = size(klist, 1);
samples = cell(Nk, 1);
first_sample = kfun(klist(1,:));
samples{1} = first_sample;
entry_count = numel(first_sample);
weight_matrix = zeros(Nk, entry_count);
weight_matrix(1,:) = abs(first_sample(:)).';
if Nk == 1
    return;
end

if use_parallel
    parfor ki = 2:Nk %#ok<PFBNS>
        sample = kfun(klist(ki,:));
        samples{ki} = sample;
        weight_matrix(ki,:) = flatten_sample(sample, entry_count);
    end
else
    for ki = 2:Nk
        sample = kfun(klist(ki,:));
        samples{ki} = sample;
        weight_matrix(ki,:) = flatten_sample(sample, entry_count);
    end
end
end

function flat_vals = flatten_sample(sample, expected_len)
flat_vals = abs(sample(:)).';
if numel(flat_vals) ~= expected_len
    error("kloop_accumulate:SampleShapeChanged", ...
        "All kfun outputs must share the same number of elements.");
end
end

function accumulator = accumulate_from_samples(samples, accumulator, selection)
if nargin < 3
    selection = 1:numel(samples);
end
for idx = selection(:).'
    accumulator = accumulator + samples{idx};
end
end

function mesh_shape = validate_mesh_shape(mesh_shape, Nk)
if isempty(mesh_shape)
    error("kloop_accumulate:MeshShapeMissing", ...
        "optionsAdapt.mesh_shape must be provided when refinement is enabled.");
end
validateattributes(mesh_shape, {"numeric"}, {"numel", 3, "positive", "finite", "integer"}, mfilename, "mesh_shape");
mesh_shape = double(mesh_shape(:).');
if prod(mesh_shape) ~= Nk
    error("kloop_accumulate:MeshShapeMismatch", ...
        "mesh_shape product (%d) must equal numel(klist) (%d).", prod(mesh_shape), Nk);
end
end

function offsets = build_refine_offsets(klist, mesh_shape, refine_factor)
step_vectors = infer_step_vectors(klist, mesh_shape);
subdivisions = ones(1, 3);
subdivisions(mesh_shape > 1) = refine_factor;
[g1, g2, g3] = ndgrid(0:subdivisions(1)-1, 0:subdivisions(2)-1, 0:subdivisions(3)-1);
count = numel(g1);
offsets = zeros(count, size(klist, 2));
step1 = step_vectors(1,:) / subdivisions(1);
step2 = step_vectors(2,:) / subdivisions(2);
step3 = step_vectors(3,:) / subdivisions(3);
offsets = g1(:) .* step1 + g2(:) .* step2 + g3(:) .* step3;
end

function step_vectors = infer_step_vectors(klist, mesh_shape)
Nk1 = mesh_shape(1);
Nk2 = mesh_shape(2);
Nk3 = mesh_shape(3);
dims = size(klist, 2);
step_vectors = zeros(3, dims);
if Nk1 > 1
    step_vectors(1,:) = klist(2,:) - klist(1,:);
end
if Nk2 > 1
    idx = Nk1 + 1;
    if idx <= size(klist, 1)
        step_vectors(2,:) = klist(idx,:) - klist(1,:);
    end
end
if Nk3 > 1
    idx = Nk1 * Nk2 + 1;
    if idx <= size(klist, 1)
        step_vectors(3,:) = klist(idx,:) - klist(1,:);
    end
end
end

function samples = evaluate_refined_samples(base_klist, kfun, offsets, use_parallel)
nsel = size(base_klist, 1);
samples = cell(nsel, 1);
if use_parallel
    parfor idx = 1:nsel %#ok<PFBNS>
        samples{idx} = average_refined_kernel(base_klist(idx,:), kfun, offsets);
    end
else
    for idx = 1:nsel
        samples{idx} = average_refined_kernel(base_klist(idx,:), kfun, offsets);
    end
end
end

function avg_sample = average_refined_kernel(base_k, kfun, offsets)
nref = size(offsets, 1);
avg_sample = kfun(base_k + offsets(1,:));
for s = 2:nref
    avg_sample = avg_sample + kfun(base_k + offsets(s,:));
end
avg_sample = avg_sample ./ nref;
end
end
