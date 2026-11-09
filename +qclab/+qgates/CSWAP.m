% CSWAP - Controlled SWAP gate (Fredkin gate)
% The CSWAP class implements a 3-qubit controlled SWAP gate with a control
% qubit, two target qubits, and an optional control state (0 or 1).
%
% The matrix representation is
%     [ 1  0  0  0  0  0  0  0;
%       0  1  0  0  0  0  0  0;
%       0  0  1  0  0  0  0  0;
%       0  0  0  1  0  0  0  0;
%       0  0  0  0  1  0  0  0;
%       0  0  0  0  0  0  1  0;
%       0  0  0  0  0  1  0  0;
%       0  0  0  0  0  0  0  1 ]
% if control < target1 < target2 and controlState = 1.
%
% Creation
%   Syntax
%     G = qclab.qgates.CSWAP()
%       - Creates a CSWAP gate with control = 0, target1 = 1, target2 = 2,
%         and controlState = 1.
%
%     G = qclab.qgates.CSWAP(control, target1, target2)
%       - Creates a CSWAP gate with the given control and target qubits,
%         using controlState = 1 by default.
%
%     G = qclab.qgates.CSWAP(control, target1, target2, controlState)
%       - Creates a CSWAP gate with the given control, target qubits,
%         and the specified control state (0 or 1).
%
% Input Arguments
%     control        - Integer index of control qubit.
%     target1        - Integer index of first target qubit.
%     target2        - Integer index of second target qubit.
%     controlState   - (optional) Control state (0 or 1). Default: 1.
%
% Output
%     G - A quantum object of type `CSWAP`, representing the controlled SWAP
%         (Fredkin) gate.
%
% Examples:
%   CSWAP gate with default parameters:
%     G = qclab.qgates.CSWAP();
%
%   CSWAP gate with control=0, target1=2, target2=4:
%     G = qclab.qgates.CSWAP(0, 2, 4);
%
%   CSWAP gate triggered on controlState = 0:
%     G = qclab.qgates.CSWAP(1, 2, 3, 0);

%> @file CSWAP.m
%> @brief Implements CSWAP class (Fredkin gate)
% ==============================================================================
%> @class CSWAP
%> @brief Controlled SWAP gate (CSWAP / Fredkin gate).
%>
%> 3-qubit CSWAP gate with default matrix representation:
%>
%> \f[\begin{bmatrix}
%>   1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
%>   0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
%>   0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
%>   0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
%>   0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
%>   0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
%>   0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
%>   0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
%> \end{bmatrix}, \f]
%> if `control` < `target1` < `target2` and `controlState` = 1.
%>
% (C) Copyright 2021.
% ==============================================================================
classdef CSWAP < qclab.qgates.QControlledGate3

  properties (Access = protected)
    %> Property storing the 2-qubit SWAP gate of this 3-qubit controlled gate.
    gate_(1,1)  qclab.qgates.SWAP
  end

  methods

    % Class constructor  =======================================================
    %> @brief Constructor for CSWAP gates.
    %>
    %> Constructs a CSWAP gate with given control qubit `control`,
    %> first target qubit `target1`, second target qubit `target2`,
    %> and control state `controlState`.
    %> The default control state is 1.
    %>
    %> @param control      (optional) control qubit index, default = 0
    %> @param target1      (optional) first target qubit index, default = 1
    %> @param target2      (optional) second target qubit index, default = 2
    %> @param controlState (optional) control state for controlled gate: 0 or 1,
    %>        default = 1
    % ==========================================================================
    function obj = CSWAP(control, target1, target2, controlState)
      % set defaults
      if nargin < 1, control = 0; end
      if nargin < 2, target1 = 1; end
      if nargin < 3, target2 = 2; end
      if nargin < 4, controlState = 1; end

      obj@qclab.qgates.QControlledGate3(control, controlState);
      obj.gate_ = qclab.qgates.SWAP(target1, target2);

      assert( qclab.isNonNegInteger(target1) );
      assert( qclab.isNonNegInteger(target2) );
      assert( target1 ~= target2 );
      assert( control ~= target1 );
      assert( control ~= target2 );
    end

    % toQASM
    function [out] = toQASM(obj, fid, offset)
      if nargin == 2, offset = 0; end
      qclab.IO.qasmCSWAP( fid, obj.control + offset, ...
                          obj.target1 + offset, obj.target2 + offset, ...
                          obj.controlState );
      out = 0;
    end

    % equals
    function [bool] = equals(obj, other)
      bool = false;
      if isa(other, 'qclab.qgates.CSWAP') && ...
          (obj.controlState == other.controlState)
        % Gates are equal if control is consistently before or after both
        % targets (same relative orientation)
        obj_ctrl_before  = (obj.control  < obj.target1)  && ...
                           (obj.control  < obj.target2);
        other_ctrl_before = (other.control < other.target1) && ...
                            (other.control < other.target2);
        obj_ctrl_after   = (obj.control  > obj.target1)  && ...
                           (obj.control  > obj.target2);
        other_ctrl_after  = (other.control > other.target1) && ...
                            (other.control > other.target2);
        bool = (obj_ctrl_before && other_ctrl_before) || ...
               (obj_ctrl_after  && other_ctrl_after);
      end
    end

    % target1 - index of first SWAP target qubit (lower index, since SWAP sorts)
    function [t] = target1(obj)
      q = obj.gate_.qubits();
      t = q(1);
    end

    % target2 - index of second SWAP target qubit (higher index)
    function [t] = target2(obj)
      q = obj.gate_.qubits();
      t = q(2);
    end

    %> Copy of 2-qubit SWAP gate of controlled-SWAP gate
    function [gate] = gate(obj)
      gate = copy(obj.gate_);
    end

    % setTarget1 - update the lower target qubit of the inner SWAP gate
    function setTarget1(obj, target1)
      assert( qclab.isNonNegInteger(target1) );
      assert( target1 ~= obj.control() );
      assert( target1 ~= obj.target2() );
      % SWAP.setQubits sorts internally, so pass new pair
      obj.gate_.setQubits([target1, obj.target2()]);
    end

    % setTarget2 - update the higher target qubit of the inner SWAP gate
    function setTarget2(obj, target2)
      assert( qclab.isNonNegInteger(target2) );
      assert( target2 ~= obj.control() );
      assert( target2 ~= obj.target1() );
      obj.gate_.setQubits([obj.target1(), target2]);
    end

    % qubits - returns sorted qubit indices spanned by this gate
    function [qubits] = qubits(obj)
      qubits = sort([obj.control, obj.target1, obj.target2]);
    end

    % ==========================================================================
    %> @brief Draw a 3-qubit CSWAP gate.
    %>
    %> @param obj       3-qubit CSWAP gate
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
      qclab.drawCommands; % load draw commands

      allQubits = obj.qubits;
      minQ = allQubits(1);
      maxQ = allQubits(end);
      span  = maxQ - minQ + 1;
      gateCell = cell(3 * span, 1);

      % Helper: 1-based row index for a given qubit
      rowOf = @(q) 3 * (q - minQ) + 1;

      % Fill all wire rows with a plain pass-through by default
      for i = 1:3:length(gateCell)
        gateCell{i}   = [ space, space ];
        gateCell{i+1} = [ h, h ];
        gateCell{i+2} = [ space, space ];
      end

      ctrlQ = obj.control;   % qubit index (avoid shadowing draw command 'ctrl')
      t1    = obj.target1;
      t2    = obj.target2;

      % --- Control qubit ---
      r = rowOf(ctrlQ);
      if ctrlQ == minQ
        gateCell{r}   = [ space, space ];
      else
        gateCell{r}   = [ space, v ];
      end
      if obj.controlState == 0
        gateCell{r+1} = [ h, ctrlo ];
      else
        gateCell{r+1} = [ h, ctrl ];
      end
      if ctrlQ == maxQ
        gateCell{r+2} = [ space, space ];
      else
        gateCell{r+2} = [ space, v ];
      end

      % --- First swap target ---
      r = rowOf(t1);
      if t1 == minQ
        gateCell{r}   = [ space, space ];
      else
        gateCell{r}   = [ space, v ];
      end
      gateCell{r+1} = [ h, sw ];
      if t1 == maxQ
        gateCell{r+2} = [ space, space ];
      else
        gateCell{r+2} = [ space, v ];
      end

      % --- Second swap target ---
      r = rowOf(t2);
      if t2 == minQ
        gateCell{r}   = [ space, space ];
      else
        gateCell{r}   = [ space, v ];
      end
      gateCell{r+1} = [ h, sw ];
      if t2 == maxQ
        gateCell{r+2} = [ space, space ];
      else
        gateCell{r+2} = [ space, v ];
      end

      % --- Gap qubits: vertical pass-through ---
      for q = minQ+1 : maxQ-1
        if q == ctrlQ || q == t1 || q == t2, continue; end
        r = rowOf(q);
        gateCell{r}   = [ space, v ];
        gateCell{r+1} = [ h, vh ];
        gateCell{r+2} = [ space, v ];
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
    %> @brief Save a 3-qubit CSWAP gate to a Tex file.
    %>
    %> @param obj       3-qubit CSWAP gate
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

      % Default: pass-through wires
      for i = 1:span
        gateCell{i} = '&\t\\qw\t';
      end

      ctrlQ   = obj.control;
      t1      = obj.target1;
      t2      = obj.target2;

      % 1-based row indices within the cell array
      ctrlRow = ctrlQ - minQ + 1;
      topSwap = min(t1, t2) - minQ + 1;   % upper swap row (1-based)
      botSwap = max(t1, t2) - minQ + 1;   % lower swap row (1-based)

      % Distance between the two swap qubits (used by \qwx)
      diffSwap = botSwap - topSwap;

      % Top swap qubit: \qswap with a wire connecting down to the bottom
      gateCell{topSwap} = ['&\t\\qswap\\qwx[', num2str(diffSwap), ']\t'];
      % Bottom swap qubit: plain \qswap
      gateCell{botSwap} = '&\t\\qswap\t';

      % Control: point to the nearest swap target
      nearSwap = [topSwap, botSwap];
      [~, nearIdx] = min(abs(nearSwap - ctrlRow));
      diffCtrl = nearSwap(nearIdx) - ctrlRow;   % +ve = down, -ve = up

      if obj.controlState == 0
        gateCell{ctrlRow} = ['&\t\\ctrlo{', num2str(diffCtrl), '}\t'];
      else
        gateCell{ctrlRow} = ['&\t\\ctrl{',  num2str(diffCtrl), '}\t'];
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

    function [bool] = fixed
      bool = true;
    end

  end

  methods (Access = protected)

    %> @brief Override copyElement to allow correct deep copy of handle property.
    function cp = copyElement(obj)
      cp = copyElement@matlab.mixin.Copyable(obj);
      cp.gate_ = obj.gate();
    end

  end

end % CSWAP