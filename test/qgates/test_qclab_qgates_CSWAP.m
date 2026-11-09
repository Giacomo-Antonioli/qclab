classdef test_qclab_qgates_CSWAP < matlab.unittest.TestCase
  methods (Test)
    function test_CSWAP(test)
      g = qclab.qgates.CSWAP();
      test.verifyEqual( g.nbQubits, int64(3) );
      test.verifyTrue( g.fixed );
      test.verifyTrue( g.controlled );
      test.verifyEqual( g.control, int64(0) );
      test.verifyEqual( g.target1, int64(1) );
      test.verifyEqual( g.target2, int64(2) );
      test.verifyEqual( g.controlState, int64(1) );

      % matrix for control < target1 < target2, controlState = 1
      M = eye(8);
      M([6 7],:) = M([7 6],:);
      test.verifyEqual( g.matrix, M );

      % controlState = 0 swaps when the control is |0>
      g0 = qclab.qgates.CSWAP(0, 1, 2, 0);
      M0 = eye(8);
      M0([2 3],:) = M0([3 2],:);
      test.verifyEqual( g0.matrix, M0 );

      % equals and copy
      test.verifyTrue( g.equals( qclab.qgates.CSWAP(0, 1, 2) ) );
      test.verifyFalse( g.equals( g0 ) );
      c = copy(g);
      c.setTarget2(3);
      test.verifyEqual( g.target2, int64(2) );

      % qubits are reported sorted
      test.verifyEqual( qclab.qgates.CSWAP(2, 0, 1).qubits, int64([0 1 2]) );

      % toQASM
      [T, out] = evalc('g.toQASM(1)');
      test.verifyEqual( out, 0 );
      test.verifyEqual( strtrim(T), 'cswap q[0], q[1], q[2];' );
      [T, ~] = evalc('g0.toQASM(1)');
      test.verifyEqual( strsplit(strtrim(T), newline), ...
        {'x q[0];', 'cswap q[0], q[1], q[2];', 'x q[0];'} );
    end

    function test_CSWAP_in_circuit(test)
      % control not adjacent to the targets
      circuit = qclab.QCircuit(4);
      circuit.push_back( qclab.qgates.CSWAP(0, 2, 3) );
      U = circuit.matrix;
      for k = 0:15
        bits = bitget(k, 4:-1:1);       % qubit 0 is the most significant
        if bits(1) == 1, bits([3 4]) = bits([4 3]); end
        kk = sum(bits .* 2.^(3:-1:0));
        test.verifyEqual( U(kk+1, k+1), 1 );
      end

      % control between the targets
      circuit = qclab.QCircuit(3);
      circuit.push_back( qclab.qgates.CSWAP(1, 0, 2) );
      U = circuit.matrix;
      for k = 0:7
        bits = bitget(k, 3:-1:1);
        if bits(2) == 1, bits([1 3]) = bits([3 1]); end
        kk = sum(bits .* 2.^(2:-1:0));
        test.verifyEqual( U(kk+1, k+1), 1 );
      end
    end
  end
end
