classdef test_qclab_QCircuit_stats < matlab.unittest.TestCase
  methods (Test)
    function test_counts_and_depth(test)
      c = qclab.QCircuit(4);
      test.verifyEqual( c.getGateCount, 0 );
      test.verifyEqual( c.getDepth, 0 );
      test.verifyEqual( c.getQubitCount, int64(4) );

      c.push_back( qclab.qgates.Hadamard(0) );
      c.push_back( qclab.qgates.Hadamard(1) );     % parallel with the first
      c.push_back( qclab.qgates.CNOT(0, 1) );
      c.push_back( qclab.qgates.MCX([0 1], 2, [1 1]) );
      c.push_back( qclab.qgates.CSWAP(0, 2, 3) );
      c.push_back( qclab.qgates.RotationY(3, 0.3) );
      test.verifyEqual( c.getGateCount, 6 );
      test.verifyEqual( c.getQubitCount, int64(4) );
      % layers: H H | CNOT | MCX | CSWAP | RY(3) -> CSWAP and RY share qubit 3
      test.verifyEqual( c.getDepth, 5 );
      s = c.getControlledStats;
      test.verifyEqual( s.single, 2 );             % CNOT, CSWAP
      test.verifyEqual( s.multi, 1 );              % MCX
      test.verifyEqual( s.total, 3 );
    end

    function test_barrier_is_ignored(test)
      c = qclab.QCircuit(2);
      c.push_back( qclab.qgates.Hadamard(0) );
      c.barrier();
      c.push_back( qclab.qgates.Hadamard(1) );
      test.verifyEqual( c.getGateCount, 2 );
      test.verifyEqual( c.getDepth, 1 );
    end

    function test_nested_circuit(test)
      inner = qclab.QCircuit(2, 1);                % acts on qubits 1 and 2
      inner.push_back( qclab.qgates.CNOT(0, 1) );
      outer = qclab.QCircuit(3);
      outer.push_back( qclab.qgates.Hadamard(0) );
      outer.push_back( inner );
      test.verifyEqual( outer.getGateCount, 2 );
      test.verifyEqual( outer.getDepth, 1 );
      test.verifyEqual( outer.getControlledStats.single, 1 );
    end
  end
end
