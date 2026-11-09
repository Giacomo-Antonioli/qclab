function circuit = mottonen_state_prep(psi, offset)
%MOTTONEN_STATE_PREP  Möttönen et al. state preparation as a QCLAB circuit.
%
%   circuit = MOTTONEN_STATE_PREP(psi) returns a qclab.QCircuit on
%   n = log2(length(psi)) qubits such that applying the circuit to
%   |0...0> reproduces psi exactly (up to the exact global phase).
%
%   circuit = MOTTONEN_STATE_PREP(psi, offset) uses an offset so that
%   the prepared state is encoded on qubits offset:offset+n-1,
%   with the rest of the (n+offset)-qubit register untouched. All gate
%   targets and controls are shifted by offset.
%
%   QUBIT CONVENTION (QCLAB): qubit 0 is the most significant qubit
%   (drawn on top). This matches the "qubit 1 = MSB" convention used to
%   derive the rotation angles, shifted by one (math qubit k -> QCLAB
%   qubit k-1).
%
%   Reference: Möttönen, Vartiainen, Bergholm, Salomaa, "Transformation
%   of quantum states using uniformly controlled rotations",
%   arXiv:quant-ph/0407010.
%
%   Each uniformly-controlled (multiplexed) rotation with target qubit t
%   and s controls is built recursively: split the 2^s angles into the
%   top-control-qubit's two halves, take the elementwise average
%   (a_plus) and half-difference (a_minus), recurse on the remaining
%   s-1 controls for each half, and sandwich as:
%
%       [recurse on a_minus]  CNOT  [recurse on a_plus]  CNOT
%
%   (applied left-to-right in time: a_plus block first).
%
%   GLOBAL PHASE
%   RY, RZ, and CNOT all have determinant +-1 when embedded in the full
%   N x N space, so a circuit built purely from them cannot realise a
%   generic global phase e^{i*phi}*I_N. The Phase gate P(phi) =
%   diag(1, e^{i*phi}) breaks this: sandwiching it with PauliX on the
%   same qubit gives e^{i*phi}*I_2 on that qubit. This sequence is
%   appended on qubit `offset` at the end.

    if nargin < 2
        offset = 0;
    end
    % Validate offset
    if ~isscalar(offset) || offset < 0 || offset ~= floor(offset)
        error('mottonen_state_prep:badOffset', ...
            'offset must be a non‑negative integer.');
    end

    psi = psi(:);
    N = numel(psi);
    n = round(log2(N));

    if 2^n ~= N
        error('mottonen_state_prep:badLength', ...
            'Length of psi must be a power of two (got %d).', N);
    end

    nrm = norm(psi);
    if abs(nrm - 1) > 1e-8
        warning('mottonen_state_prep:notNormalized', ...
            '||psi|| = %.6g != 1; renormalizing.', nrm);
    end
    psi = psi / nrm;

    %% ------------------------------------------------------------
    % Magnitude/phase tree (identical to the dense-matrix version):
    % level k (k = n downto 1) yields the 2^(k-1) angles for the
    % multiplexed rotation that targets qubit k (math, 1-indexed).
    %% ------------------------------------------------------------
    r_k   = abs(psi);
    phi_k = angle(psi);

    ry_angles = cell(1, n);
    rz_angles = cell(1, n);

    for k = n:-1:1
        R   = reshape(r_k,   2, []);
        Phi = reshape(phi_k, 2, []);

        ry_angles{k} = (2 * atan2(R(2, :), R(1, :))).';
        rz_angles{k} = (Phi(2, :) - Phi(1, :)).';

        r_k   = hypot(R(1, :), R(2, :)).';
        phi_k = ((Phi(1, :) + Phi(2, :)) / 2).';
    end

    global_phase = phi_k(1);

    %% ------------------------------------------------------------
    % Build the circuit. Math qubit k -> QCLAB qubit (k-1)+offset.
    % Controls for level-k are math qubits 1..k-1, i.e.
    % QCLAB qubits offset .. offset+k-2, in that order.
    %% ------------------------------------------------------------
    circuit = qclab.QCircuit(n + offset);

    % Add Y-rotations (first pass)
    for k = 1:n
        target   = (k - 1) + offset;
        controls = (0:(k - 2)) + offset;   % empty when k == 1
        add_multiplexed_rotation(circuit, 'ry', target, controls, ry_angles{k});
    end

    % Add Z-rotations (second pass)
    for k = 1:n
        target   = (k - 1) + offset;
        controls = (0:(k - 2)) + offset;
        add_multiplexed_rotation(circuit, 'rz', target, controls, rz_angles{k});
    end

    % --- mandatory global-phase correction (see header) ---
    circuit.push_back(qclab.qgates.Phase(offset, global_phase));
    circuit.push_back(qclab.qgates.PauliX(offset));
    circuit.push_back(qclab.qgates.Phase(offset, global_phase));
    circuit.push_back(qclab.qgates.PauliX(offset));

    %% ================================================================
    %  Nested helper function (recursive decomposition of UCR)
    %  (required to handle the recursive sum/difference + CNOT sandwich)
    %% ================================================================
    function add_multiplexed_rotation(circuit, gate_type, target, controls, angles)
        % Appends the uniformly-controlled rotation to `circuit`.
        % `controls(1)` is the most-significant control qubit,
        % matching angles' bit-pattern indexing.
        s = numel(controls);

        % Base case: no controls -> single rotation
        if s == 0
            if strcmp(gate_type, 'ry')
                circuit.push_back(qclab.qgates.RotationY(target, angles(1)));
            else  % 'rz'
                circuit.push_back(qclab.qgates.RotationZ(target, angles(1)));
            end
            return;
        end

        half = numel(angles) / 2;
        a = angles(1:half);           % top control (controls(1)) == 0
        b = angles(half+1:end);       % top control (controls(1)) == 1

        a_plus  = (a + b) / 2;
        a_minus = (a - b) / 2;

        c1   = controls(1);
        rest = controls(2:end);

        % a_plus branch (applied first in time)
        add_multiplexed_rotation(circuit, gate_type, target, rest, a_plus);
        circuit.push_back(qclab.qgates.CNOT(c1, target));
        % a_minus branch
        add_multiplexed_rotation(circuit, gate_type, target, rest, a_minus);
        circuit.push_back(qclab.qgates.CNOT(c1, target));
    end
end