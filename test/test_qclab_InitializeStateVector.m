classdef test_qclab_InitializeStateVector < matlab.unittest.TestCase
  methods (Test)
    function test_mottonen_real_and_complex(test)
      rng(7);
      for n = 1:4
        psi = randn(2^n, 1) + 1i * randn(2^n, 1);
        psi = psi / norm(psi);
        c = qclab.mottonen_state_prep(psi, 0);
        e0 = zeros(2^n, 1); e0(1) = 1;
        out = c.matrix * e0;
        test.verifyEqual( out, psi, 'AbsTol', 1e-10 );
      end
    end

    function test_initialize_in_circuit(test)
      psi = [1; 2; 3; 4] / sqrt(30);
      c = qclab.QCircuit(2);
      c.InitializeStateVector(psi);
      out = c.matrix * [1; 0; 0; 0];
      test.verifyEqual( out, psi, 'AbsTol', 1e-10 );
    end

    function test_initialize_with_offset_and_label(test)
      psi = [1; 1i] / sqrt(2);
      c = qclab.QCircuit(3);
      c.InitializeStateVector(psi, 1, 'prep');
      out = c.matrix * [1; zeros(7, 1)];
      expected = kron(kron([1; 0], psi), [1; 0]);
      test.verifyEqual( out, expected, 'AbsTol', 1e-10 );
    end

    function test_bad_length(test)
      test.verifyError( @() qclab.mottonen_state_prep([1; 2; 3], 0), ...
        'mottonen_state_prep:badLength' );
    end
  end
end
