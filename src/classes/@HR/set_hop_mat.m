function H_hr = set_hop_mat(H_hr, amp, vector, mode)
%SET_HOP_MAT Sets or modifies the hopping matrix in a HR object.
%   H_hr = SET_HOP_MAT(H_hr, amp, vector, mode) updates the hopping matrix
%   component of the tight-binding model object H_hr based on the input
%   amplitude AMP, lattice translation VECTOR, and operation MODE.
%
%   Inputs:
%       H_hr   - TBmodel object with hopping matrix fields.
%       amp    - Hopping amplitude(s), numeric or symbolic, matrix or scalar.
%       vector - Lattice translation vector [Rx, Ry, Rz] (or higher dimension).
%       mode   - String specifying operation mode:
%                  'set'     : Overwrites existing hopping matrix.
%                  'add'     : Adds to existing numeric hopping matrix.
%                  'sym'     : Overwrites symbolic hopping matrix (not implemented for sparse).
%                  'symadd'  : Adds to symbolic hopping matrix (not implemented for sparse).
%
%   Output:
%       H_hr   - Updated HR object with modified hopping matrices.

V = H_hr.vectorL;

switch H_hr.Type
    case 'list'
        sizeamp = size(amp);
        for n = 1:numel(amp)
            [hi, hj] = ind2sub(sizeamp, n);
            if ~isequal(sym(amp(n)), sym(0))
                H_hr = set_hop_single(H_hr, amp(n), hi, hj, vector, mode);
            end
        end

    case 'sparse'
        [~, seq] = ismember(vector, V, 'rows');
        if seq == 0
            seq = H_hr.NRPTS + 1;
            H_hr = H_hr.add_empty_one(vector);
        end

        switch mode
            case 'set'
                H_hr.HnumL{seq} = amp;

            case 'add'
                H_hr.HnumL{seq} = H_hr.HnumL{seq} + amp;

            case 'sym'
                error('Symbolic operations not implemented for sparse format.');

            case 'symadd'
                error('Symbolic operations not implemented for sparse format.');
        end

    otherwise
        [~, seq] = ismember(vector, V, 'rows');
        if seq == 0
            seq = H_hr.NRPTS + 1;
            H_hr = H_hr.add_empty_one(vector);
        end

        switch mode
            case 'set'
                zeros_num_mat = zeros(H_hr.WAN_NUM);
                if any(H_hr.HnumL(:,:,seq) ~= zeros_num_mat, 'all')
                    warning('Hopping already exists. Consider using ''add'' mode.');
                    Format = fold(@strcat, repmat("%d ", [1, H_hr.Dim]));
                    Input = num2cell(vector(1,:));
                    fprintf(strcat(Format, "\n"), Input{:});
                end
                H_hr.HnumL(:,:,seq) = amp;

            case 'add'
                H_hr.HnumL(:,:,seq) = H_hr.HnumL(:,:,seq) + amp;

            case 'sym'
                zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
                if any(H_hr.HcoeL(:,:,seq) ~= zeros_coe_mat, 'all')
                    warning('Symbolic hopping already exists. Consider using ''symadd'' mode.');
                    Format = fold(@strcat, repmat("%d ", [1, H_hr.Dim]));
                    Input = num2cell(vector(1,:));
                    fprintf(strcat(Format, "\n"), Input{:});
                end
                H_hr.HcoeL(:,:,seq) = amp;

            case 'symadd'
                H_hr.HcoeL(:,:,seq) = H_hr.HcoeL(:,:,seq) + amp;
        end
end
end
