%> @file QControlledGate3.m
%> @brief Implements QControlledGate3 class
% ==============================================================================
%> @class QControlledGate3
%> @brief Base class for 3-qubit gates of controlled 2-qubit gates.
%>
%> This class models a 3-qubit quantum gate consisting of one control qubit
%> and a 2-qubit gate applied to two target qubits, conditioned on the
%> control qubit being in state `controlState` (0 or 1).
%>
% (C) Copyright Daan Camps, Sophia Keip and Roel Van Beeumen 2025.
% ==============================================================================
classdef QControlledGate3 < qclab.qgates.QGate3

  properties (Access = protected)
    %> Control qubit of this 3-qubit gate.
    control_(1,1)      int64
    %> Control state of this 3-qubit gate.
    controlState_(1,1) int64
  end

  methods

    % Class constructor  =======================================================
    %> @brief Constructor for 3-qubit gates of controlled 2-qubit gates.
    %>
    %> Constructs a 3-qubit gate with the given control qubit `control` on the
    %> control state `controlState`. The default control state is 1.
    %>
    %> @param control      control qubit index
    %> @param controlState (optional) control state for controlled gate: 0 or 1.
    %>        Default: 1.
    % ==========================================================================
    function obj = QControlledGate3(control, controlState)
      if nargin < 2, controlState = 1; end

      obj.control_      = control;
      obj.controlState_ = controlState;

      assert( qclab.isNonNegInteger(control) );
      assert( (controlState == 0) || (controlState == 1) );
    end

    % qubit - returns the lowest qubit index spanned by this gate
    function [qubit] = qubit(obj)
      qubit = min([obj.control_, obj.target1, obj.target2]);
    end

    % qubits - returns sorted qubit indices [min, ..., max] spanned by this gate
    function [qubits] = qubits(obj)
      qubits = sort([obj.control_, obj.target1, obj.target2]);
    end

    % setQubits
    function setQubits(obj, qubits)
      assert( qclab.isNonNegIntegerArray(qubits) );
      assert( numel(unique(qubits)) == 3 );
      obj.control_ = qubits(1);
      obj.setTarget1( qubits(2) );
      obj.setTarget2( qubits(3) );
    end

    % matrix
    %> @brief Compute the spanning matrix of size 2^span x 2^span, where
    %> span = maxQ - minQ + 1 and maxQ/minQ are the highest/lowest qubit indices
    %> among control, target1, target2.
    %>
    %> Gap qubits between minQ and maxQ act as identity. This is consistent with
    %> QControlledGate2.matrix() which uses Imid to pad gap qubits.
    function [mat] = matrix(obj)
      E0 = [1 0; 0 0];
      E1 = [0 0; 0 1];
      I4 = eye(4);
      CG = obj.gate.matrix;  % 4x4 matrix of the controlled 2-qubit gate

      % Sorted qubit positions and control position index (0=MSB, 2=LSB)
      allQ    = sort([obj.control_, obj.target1, obj.target2]);
      ctrlPos = find(allQ == obj.control_) - 1;  % 0, 1, or 2

      % Build 8x8 local matrix for the 3 active qubits
      switch ctrlPos
        case 0
          if obj.controlState_ == 0
            mat8 = kron(E0, CG) + kron(E1, I4);
          else
            mat8 = kron(E0, I4) + kron(E1, CG);
          end
        case 2
          if obj.controlState_ == 0
            mat8 = kron(CG, E0) + kron(I4, E1);
          else
            mat8 = kron(I4, E0) + kron(CG, E1);
          end
        case 1
          % Control is sandwiched: partition by middle bit value
          idx0 = [1, 2, 5, 6];
          idx1 = [3, 4, 7, 8];
          mat8 = zeros(8, 8);
          if obj.controlState_ == 0
            mat8(idx0, idx0) = CG;
            mat8(idx1, idx1) = I4;
          else
            mat8(idx0, idx0) = I4;
            mat8(idx1, idx1) = CG;
          end
        otherwise
          error('QControlledGate3:matrix', 'Unexpected ctrlPos value %d.', ctrlPos);
      end

      % Embed into 2^span x 2^span, placing identity on any gap qubits
      allQ_d = double(allQ);
      minQ   = allQ_d(1);
      maxQ   = allQ_d(3);
      span   = maxQ - minQ + 1;   % pure double arithmetic
      if span == 3
        mat = mat8;
      else
        mat = qclab.qgates.QControlledGate3.buildSpanningMatrix_( ...
          mat8, minQ, allQ_d, span);
      end
    end

    % ==========================================================================
    %> @brief Apply the QControlledGate3 to a matrix or a struct of state
    %> vectors.
    %>
    %> @param obj      instance of QControlledGate3 class.
    %> @param side     'L' or 'R' for left or right application.
    %> @param op       'N', 'T' or 'C' for normal, transpose, or conjugate
    %>                 transpose application.
    %> @param nbQubits qubit size of `current`.
    %> @param current  matrix or struct of state vectors to which the gate is
    %>                 applied.
    %> @param offset   offset applied to qubits. Default: 0.
    % ==========================================================================
    function [current] = apply(obj, side, op, nbQubits, current, offset)
      if nargin == 5, offset = 0; end
      isSparse = qclab.isSparse(nbQubits);
      assert( nbQubits >= 3 );
      if isa(current, 'double')
        if strcmp(side, 'L')
          assert( size(current, 2) == 2^nbQubits );
        else
          assert( size(current, 1) == 2^nbQubits );
        end
      else
        assert( length(current.states{1}) == 2^nbQubits );
      end

      qubits = obj.qubits + offset;  % sorted [minQ, midQ, maxQ]
      assert( qubits(1) < nbQubits );
      assert( qubits(3) < nbQubits );

      % Get the spanning matrix (2^span x 2^span) with correct op variant.
      % obj.matrix already handles gap qubits and returns the right size.
      switch op
        case 'N',  mats = obj.matrix;
        case 'T',  mats = obj.matrix.';
        otherwise, mats = obj.matrix';
      end

      minQ = qubits(1);
      maxQ = qubits(3);
      span = maxQ - minQ + 1;

      % Embed in the full n-qubit space: kron(I_left, mats, I_right)
      if ( minQ == 0 && maxQ == nbQubits - 1 )
        matn = mats;
      elseif ( minQ == 0 )
        matn = kron(mats, qclab.qId(nbQubits - span, isSparse));
      elseif ( maxQ == nbQubits - 1 )
        matn = kron(qclab.qId(nbQubits - span, isSparse), mats);
      else
        matn = kron(kron(qclab.qId(minQ, isSparse), mats), ...
                    qclab.qId(nbQubits - maxQ - 1, isSparse));
      end

      current = qclab.applyGateTo(current, matn, side);
    end

    %> @brief Returns the control qubit of this 3-qubit gate.
    function [control] = control(obj)
      control = obj.control_;
    end

    %> @brief Returns the control state of this 3-qubit gate.
    function [controlState] = controlState(obj)
      controlState = obj.controlState_;
    end

    %> @brief Sets the control qubit of this controlled gate.
    function setControl(obj, control)
      assert( qclab.isNonNegInteger(control) );
      assert( control ~= obj.target1() );
      assert( control ~= obj.target2() );
      obj.control_ = control;
    end

    %> @brief Sets the control state of this controlled gate.
    function setControlState(obj, controlState)
      assert( (controlState == 0) || (controlState == 1) );
      obj.controlState_ = controlState;
    end

    % ==========================================================================
    %> @brief Draw a 3-qubit gate of a controlled 2-qubit gate.
    %>
    %> @param obj       3-qubit gate of a controlled 2-qubit gate.
    %> @param fid       file id to draw to:
    %>                   - 0  : return cell array with ascii characters as `out`
    %>                   - 1  : draw to command window (default)
    %>                   - >1 : draw to (open) file id
    %> @param parameter 'N' don't print parameter (default), 'S' print short
    %>                  parameter, 'L' print long parameter.
    %> @param offset    qubit offset. Default is 0.
    %>
    %> @retval out if fid > 0 then out == 0 on successful completion, otherwise
    %>             out contains a cell array with the drawing info.
    % ==========================================================================
    function [varargout] = draw(obj, fid, parameter, offset)
      if nargin < 2, fid = 1; end
      if nargin < 3, parameter = 'N'; end
      if nargin < 4, offset = 0; end
      qclab.drawCommands;  % load draw commands

      allQubits = obj.qubits;
      minQ = allQubits(1);
      maxQ = allQubits(end);
      span  = maxQ - minQ + 1;
      gateCell = cell(3 * span, 1);

      label = obj.gate.label(parameter);
      width = length(label);
      if mod(width, 2) == 0
        width = width + 1;
        label = [label space];
      end
      hwidth = floor(width / 2);

      % Default: pass-through wire for every row
      for i = 1:3:length(gateCell)
        gateCell{i}   = [ repmat(space, 1, hwidth + 2), ...
                          v, repmat(space, 1, width - hwidth) ];
        gateCell{i+1} = [ repmat(h, 1, hwidth + 2), ...
                          vh, repmat(h, 1, width - hwidth) ];
        gateCell{i+2} = [ repmat(space, 1, hwidth + 2), ...
                          v, repmat(space, 1, width - hwidth) ];
      end

      ctrl = obj.control_;
      t1   = obj.target1;
      t2   = obj.target2;
      rowOf = @(q) 3 * (q - minQ) + 1;  % 1-based row for a given qubit

      % --- Control qubit ---
      r = rowOf(ctrl);
      gateCell{r}   = repmat(space, 1, width + 3);
      if obj.controlState_ == 0
        gateCell{r+1} = [ repmat(h, 1, hwidth+2), ctrlo, ...
                          repmat(h, 1, width - hwidth) ];
      else
        gateCell{r+1} = [ repmat(h, 1, hwidth+2), ctrl, ...
                          repmat(h, 1, width - hwidth) ];
      end
      % vertical connector leaving the control toward the targets
      if ctrl < min(t1, t2)
        gateCell{r+2} = [ repmat(space, 1, hwidth+2), v, ...
                          repmat(space, 1, width - hwidth) ];
      else
        gateCell{r} = [ repmat(space, 1, hwidth+2), v, ...
                        repmat(space, 1, width - hwidth) ];
        gateCell{r+2} = repmat(space, 1, width + 3);
      end

      % --- First target qubit ---
      r = rowOf(t1);
      gateCell{r}   = [ space, ulc, repmat(h, 1, hwidth), hu, ...
                        repmat(h, 1, width - hwidth - 1), urc ];
      gateCell{r+1} = [ h, vl, label, vr ];
      gateCell{r+2} = [ space, blc, repmat(h, 1, hwidth), hb, ...
                        repmat(h, 1, width - hwidth - 1), brc ];

      % --- Second target qubit ---
      r = rowOf(t2);
      gateCell{r}   = [ space, ulc, repmat(h, 1, hwidth), hu, ...
                        repmat(h, 1, width - hwidth - 1), urc ];
      gateCell{r+1} = [ h, vl, label, vr ];
      gateCell{r+2} = [ space, blc, repmat(h, 1, width), brc ];

      % Top of first target box should have no upward connector if t1 is first
      if min(t1, t2) < ctrl
        % targets are above the control — remove upward stub on upper target
        r = rowOf(min(t1, t2));
        gateCell{r} = [ space, ulc, repmat(h, 1, width), urc ];
      end

      if fid > 0
        drawQubits = (minQ:maxQ) + offset;
        qclab.drawCellArray(fid, gateCell, drawQubits);
        out = 0;
      else
        out = gateCell;
      end

      if nargout > 0
        varargout = {out};
      else
        varargout = {};
      end
    end

    % ==========================================================================
    %> @brief Save a 3-qubit gate of a controlled 2-qubit gate to a TeX file.
    %>
    %> @param obj       3-qubit gate of a controlled 2-qubit gate.
    %> @param fid       file id to draw to:
    %>                   - 0  : return cell array with ascii characters as `out`
    %>                   - 1  : draw to command window (default)
    %>                   - >1 : draw to (open) file id
    %> @param parameter 'N' don't print parameter (default), 'S' print short
    %>                  parameter, 'L' print long parameter.
    %> @param offset    qubit offset. Default is 0.
    %>
    %> @retval out if fid > 0 then out == 0 on successful completion, otherwise
    %>             out contains a cell array with the drawing info.
    % ==========================================================================
    function [out] = toTex(obj, fid, parameter, offset)
      if nargin < 2, fid = 1; end
      if nargin < 3, parameter = 'N'; end
      if nargin < 4, offset = 0; end

      allQubits = obj.qubits;
      minQ = allQubits(1);
      maxQ = allQubits(end);
      span  = maxQ - minQ + 1;
      gateCell = cell(span, 1);
      label = obj.gate.label(parameter, true);

      % Default: pass-through wires
      for i = 1:span
        gateCell{i} = '&\t\\qw\t';
      end

      ctrl    = obj.control_;
      t1      = obj.target1;
      t2      = obj.target2;
      ctrlRow = ctrl - minQ + 1;
      t1Row   = t1   - minQ + 1;
      t2Row   = t2   - minQ + 1;

      % Control — offset to nearest target
      if ctrl < min(t1, t2)
        diffCtrl = min(t1Row, t2Row) - ctrlRow;
      else
        diffCtrl = -(ctrlRow - max(t1Row, t2Row));
      end

      if obj.controlState_ == 0
        gateCell{ctrlRow} = ['&\t\\ctrlo{', num2str(diffCtrl), '}\t'];
      else
        gateCell{ctrlRow} = ['&\t\\ctrl{',  num2str(diffCtrl), '}\t'];
      end

      % Two target qubits: use multigate spanning from t1Row to t2Row
      topRow = min(t1Row, t2Row);
      botRow = max(t1Row, t2Row);
      nTargetRows = botRow - topRow + 1;  % number of rows the gate spans

      gateCell{topRow} = ['&\t\\multigate{', num2str(nTargetRows - 1), ...
                          '}{', label, '}\t'];
      for r = topRow+1 : botRow
        gateCell{r} = '&\t\\ghost{}\t';
      end

      if fid > 0
        drawQubits = (minQ:maxQ) + offset;
        qclab.toTexCellArray(fid, gateCell, drawQubits);
        out = 0;
      else
        out = gateCell;
      end
    end

  end

  methods (Static)
    %> @brief Returns true: this is a controlled gate.
    function [bool] = controlled
      bool = true;
    end

    % ==========================================================================
    %> @brief Embed an 8x8 local gate matrix into a 2^span x 2^span spanning
    %> subspace matrix, placing identity on all "gap" qubits.
    %>
    %> The local matrix mat8 acts on three qubits whose absolute indices are
    %> given in allQ (sorted ascending).  The spanning subspace covers qubit
    %> indices minQ..minQ+span-1; any index not in allQ is a gap qubit on which
    %> the operation acts as identity.
    %>
    %> Algorithm (O(2^span * 8)):
    %>   For each basis state inSpan in the spanning subspace:
    %>     1. Extract the three gate-qubit bits -> in3  (3-bit index into mat8)
    %>     2. Compute base = inSpan with those three bits cleared
    %>     3. For each out3 in 0..7:
    %>          outSpan = base | (out3 bits placed at gate positions)
    %>          mats(outSpan+1, inSpan+1) += mat8(out3+1, in3+1)
    %>
    %> @param mat8   8x8 local matrix (basis ordered by allQ ascending)
    %> @param minQ   lowest qubit index in the spanning subspace
    %> @param allQ   1x3 sorted array of the three gate qubit indices
    %> @param span   number of qubits in the spanning subspace
    %> @retval mats  2^span x 2^span spanning subspace matrix
    % ==========================================================================
    function mats = buildSpanningMatrix_(mat8, minQ, allQ, span)
      % Convert everything to double to avoid int64 division / mixed-type issues
      minQ = double(minQ);
      allQ = double(allQ(:)');   % 1x3 row vector
      span = double(span);

      N    = 2^span;
      mats = zeros(N, N);

      % Bit weight of each active qubit within the span-bit index
      % (MSB = position 0 = weight 2^(span-1))
      wA = 2 .^ (span - 1 - (allQ - minQ));   % [wA1, wA2, wA3]

      % Gap qubit positions (0-indexed from MSB) and their weights
      gap_pos  = setdiff(0 : span-1, allQ - minQ);
      num_gaps = numel(gap_pos);
      wG       = 2 .^ (span - 1 - gap_pos);   % one entry per gap qubit

      % Iterate over all 2^num_gaps gap-qubit configurations
      for g = 0 : 2^num_gaps - 1

        % Span-index contribution from gap qubits for this configuration
        gap_contrib = 0;
        for j = 1:num_gaps
          bit_j       = mod(floor(g / 2^(num_gaps - j)), 2);
          gap_contrib = gap_contrib + bit_j * wG(j);
        end

        % For every (in3, out3) pair, map to (in_full, out_full)
        for in3 = 0:7
          in_active = 0;
          for k = 1:3
            in_active = in_active + mod(floor(in3 / 2^(3-k)), 2) * wA(k);
          end
          in_full = gap_contrib + in_active;

          for out3 = 0:7
            val = mat8(out3+1, in3+1);
            if val == 0, continue; end

            out_active = 0;
            for k = 1:3
              out_active = out_active + mod(floor(out3 / 2^(3-k)), 2) * wA(k);
            end
            out_full = gap_contrib + out_active;

            mats(out_full+1, in_full+1) = val;
          end
        end
      end
    end

  end

  methods (Abstract)
    %> @brief Returns a copy of the controlled 2-qubit gate of this 3-qubit gate.
    [gate] = gate(obj)
    %> @brief Returns the first target qubit of this 3-qubit gate.
    [target] = target1(obj)
    %> @brief Returns the second target qubit of this 3-qubit gate.
    [target] = target2(obj)
    %> @brief Sets the first target qubit of this controlled gate.
    setTarget1(obj, target)
    %> @brief Sets the second target qubit of this controlled gate.
    setTarget2(obj, target)
  end

  methods (Access = protected)

    %> Property display groups
    function groups = getPropertyGroups(obj)
      import matlab.mixin.util.PropertyGroup
      props = struct();
      props.nbQubits  = obj.nbQubits;
      props.Control   = obj.control;
      props.Target1   = obj.target1;
      props.Target2   = obj.target2;
      groups = PropertyGroup(props);
    end

  end
end