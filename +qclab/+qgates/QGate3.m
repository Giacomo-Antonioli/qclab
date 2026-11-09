%> @file QGate3.m
%> @brief Implements QGate3 class
% ==============================================================================
%> @class QGate3
%> @brief Base class for 3-qubit gates.
%>
%>  This class implements functionalities that are shared between different
%>  types of 3-qubit gates.
%>
%>  Abstract base class for 3-qubit gates.
%>  Concrete subclasses have to implement following methods:
%>  - equals      % defined in qclab.QObject
%>  - toQASM      % defined in qclab.QObject
%>  - matrix      % defined in qclab.QObject
%>  - fixed       % defined in qclab.QObject
%>
% (C) Copyright Daan Camps, Sophia Keip and Roel Van Beeumen 2025.
% ==============================================================================
classdef QGate3 < qclab.QObject

  methods

    % ==========================================================================
    %> @brief Apply the QGate3 to a matrix or a struct of state vectors.
    %>
    %> @param obj      instance of QGate3 class.
    %> @param side     'L' or 'R' for respectively left or right side of
    %>                 application (in quantum circuit ordering).
    %> @param op       'N', 'T' or 'C' for respectively normal, transpose or
    %>                 conjugate transpose application of QGate3.
    %> @param nbQubits qubit size of `current`.
    %> @param current  matrix or struct of state vectors to which QGate3 is
    %>                 applied.
    %> @param offset   offset applied to qubits. Default: 0.
    % ==========================================================================
    function [current] = apply(obj, side, op, nbQubits, current, offset)
      if nargin == 5, offset = 0; end
      assert( nbQubits >= 3 );
      isSparse = qclab.isSparse(nbQubits);
      if isa(current, 'double')
        if strcmp(side, 'L') % left
          assert( size(current, 2) == 2^nbQubits );
        else                 % right
          assert( size(current, 1) == 2^nbQubits );
        end
      else
        assert( length(current.states{1}) == 2^nbQubits );
      end
      qubits = obj.qubits + offset;
      assert( qubits(1) < nbQubits );
      assert( qubits(3) < nbQubits );
      assert( qubits(1) + 1 == qubits(2) ); % nearest neighbor qubits
      assert( qubits(2) + 1 == qubits(3) ); % nearest neighbor qubits
      % operation
      if strcmp(op, 'N')      % normal
        mat3 = obj.matrix;
      elseif strcmp(op, 'T')  % transpose
        mat3 = obj.matrix.';
      else                    % conjugate transpose
        mat3 = obj.matrix';
      end
      % kron( Ileft, mat3, Iright )
      if ( nbQubits == 3 )
        matn = mat3;
      elseif ( qubits(1) == 0 )
        matn = kron(mat3, qclab.qId(nbQubits - 3, isSparse));
      elseif ( qubits(3) == nbQubits - 1 )
        matn = kron(qclab.qId(nbQubits - 3, isSparse), mat3);
      else
        matn = kron(kron(qclab.qId(qubits(1), isSparse), mat3), ...
                    qclab.qId(nbQubits - qubits(3) - 1, isSparse));
      end
      current = qclab.applyGateTo(current, matn, side);
    end

    % ctranspose
    function objprime = ctranspose(obj)
      objprime = copy(obj);
    end

    % ==========================================================================
    %> @brief Draw a 3-qubit gate acting on nearest-neighbor qubits.
    %>
    %> @param obj       3-qubit gate
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
      qubits = obj.qubits + offset;
      assert( qubits(1) + 1 == qubits(2) ); % nearest neighbor qubits
      assert( qubits(2) + 1 == qubits(3) ); % nearest neighbor qubits
      qclab.drawCommands; % load draw commands
      gateCell = cell(9, 1);
      label = obj.label(parameter);
      width = length(label);
      % Top qubit row
      gateCell{1} = [ space, ulc, repmat(h, 1, width), urc ];
      gateCell{2} = [ h, vl, repmat(space, 1, width), vr ];
      gateCell{3} = [ space, v, label, v ];
      % Middle qubit row (interior of the gate box)
      gateCell{4} = [ space, v, repmat(space, 1, width), v ];
      gateCell{5} = [ h, vl, repmat(space, 1, width), vr ];
      gateCell{6} = [ space, v, repmat(space, 1, width), v ];
      % Bottom qubit row
      gateCell{7} = [ space, v, repmat(space, 1, width), v ];
      gateCell{8} = [ h, vl, repmat(space, 1, width), vr ];
      gateCell{9} = [ space, blc, repmat(h, 1, width), brc ];
      if fid > 0
        qclab.drawCellArray(fid, gateCell, qubits);
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
    %> @brief Save a 3-qubit gate acting on nearest-neighbor qubits to TeX file.
    %>
    %> @param obj       3-qubit gate
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
      qubits = obj.qubits + offset;
      assert( qubits(1) + 1 == qubits(2) ); % nearest neighbor qubits
      assert( qubits(2) + 1 == qubits(3) ); % nearest neighbor qubits
      gateCell = cell(3, 1);
      label = obj.label(parameter, true);
      % Top qubit: multigate spanning 3 rows (index 2 = span of 2 additional rows)
      gateCell{1} = ['&\t\\multigate{2}{', label, '}\t'];
      % Middle and bottom qubits: ghost entries
      gateCell{2} = ['&\t\\ghost{', label, '}\t'];
      gateCell{3} = ['&\t\\ghost{', label, '}\t'];
      if fid > 0
        qclab.toTexCellArray(fid, gateCell, qubits);
        out = 0;
      else
        out = gateCell;
      end
    end

  end

  methods (Static)

    % nbQubits
    function [nbQubits] = nbQubits(~)
      nbQubits = int64(3);
    end

    % setQubit — not supported for multi-qubit gates
    function setQubit(~, ~)
      assert( false );
    end

  end

end % class QGate3