@testset "Frame #01" begin
    # Define an empty model:
    model = Model()

    # Define the nodes and DOF supports:
    node!(model,  1, 0.0000e+00, 0.0000e+00, 0.0000e+00, u_x = true, u_y = true, θ_z = true)
    node!(model,  2, 3.6000e+02, 0.0000e+00, 0.0000e+00, u_x = true, u_y = true, θ_z = true)
    node!(model,  3, 0.0000e+00, 4.5000e+01, 0.0000e+00)
    node!(model,  4, 3.6000e+02, 4.5000e+01, 0.0000e+00)
    node!(model,  5, 0.0000e+00, 9.0000e+01, 0.0000e+00)
    node!(model,  6, 3.6000e+02, 9.0000e+01, 0.0000e+00)
    node!(model,  7, 0.0000e+00, 1.3500e+02, 0.0000e+00)
    node!(model,  8, 3.6000e+02, 1.3500e+02, 0.0000e+00)
    node!(model,  9, 0.0000e+00, 1.8000e+02, 0.0000e+00)
    node!(model, 10, 9.0000e+01, 1.8000e+02, 0.0000e+00)
    node!(model, 11, 1.8000e+02, 1.8000e+02, 0.0000e+00)
    node!(model, 12, 2.7000e+02, 1.8000e+02, 0.0000e+00)
    node!(model, 13, 3.6000e+02, 1.8000e+02, 0.0000e+00)

    # Define the sections:
    section!(model, 1, 2.0000e+01, 7.2200e+02, 1.2100e+02, 3.0100e+00)
    section!(model, 2, 1.0300e+01, 5.1000e+02, 1.5300e+01, 5.0600e-01)

    # Define the materials:
    material!(model, 1, 29000, 0.3, 0.290)

    # Define the elements:
    element!(model,  1,  1,  3, 1, 1)
    element!(model,  2,  3,  5, 1, 1)
    element!(model,  3,  5,  7, 1, 1)
    element!(model,  4,  7,  9, 1, 1)
    element!(model,  5,  9, 10, 2, 1)
    element!(model,  6, 10, 11, 2, 1)
    element!(model,  7, 11, 12, 2, 1)
    element!(model,  8, 12, 13, 2, 1)
    element!(model,  9,  2,  4, 1, 1)
    element!(model, 10,  4,  6, 1, 1)
    element!(model, 11,  6,  8, 1, 1)
    element!(model, 12,  8, 13, 1, 1)

    # Define the loads:
    concload!(model,  9, 3.4849e+01, -3.4849e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 13, 0.0000e+00, -3.4849e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)

    # Solve the model using a linear elastic analysis:
    solve!(model, LinearElasticAnalysis())

    # Check the nodal displacements:
    @test getnodaldisplacements(model,  1) ≈ [0.0000e+00,  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  0.0000e+00] rtol = 1E-3
    @test getnodaldisplacements(model,  2) ≈ [0.0000e+00,  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  0.0000e+00] rtol = 1E-3
    @test getnodaldisplacements(model,  3) ≈ [8.8550e-02, -2.6579e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -3.6521e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  4) ≈ [8.6543e-02, -2.7497e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -3.5681e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  5) ≈ [3.0317e-01, -5.3158e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -5.6032e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  6) ≈ [2.9609e-01, -5.4994e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -5.4667e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  7) ≈ [5.6733e-01, -7.9737e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -5.8535e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  8) ≈ [5.5350e-01, -8.2491e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -5.6958e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  9) ≈ [8.0447e-01, -1.0632e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00, -4.4028e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 10) ≈ [7.9927e-01, -2.5797e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00,  4.9284e-04] rtol = 1E-3
    @test getnodaldisplacements(model, 11) ≈ [7.9407e-01, -1.1478e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00,  2.1493e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 12) ≈ [7.8887e-01,  3.1722e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  5.6651e-04] rtol = 1E-3
    @test getnodaldisplacements(model, 13) ≈ [7.8367e-01, -1.0999e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00, -4.2555e-03] rtol = 1E-3

    # Check the nodal reactions:
    # @test getnodalreactions(model, solution, 1) ≈ [-1.7587e+01, 3.4257e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 2.0950e+03] rtol = 1E-3
    # @test getnodalreactions(model, solution, 2) ≈ [-1.7262e+01, 3.5440e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 2.0486e+03] rtol = 1E-3
end

@testset "Frame #05" begin
    # Define an empty model:
    model = Model()

    # Define the nodes and DOF supports:
    node!(model,  1, 0.0000e+00, 0.0000e+00, 0.0000e+00, u_x = true, u_y = true)
    node!(model,  2, 2.8800e+02, 0.0000e+00, 0.0000e+00, u_x = true, u_y = true)
    node!(model,  3, 5.0400e+02, 0.0000e+00, 0.0000e+00, u_x = true, u_y = true)
    node!(model,  4, 0.0000e+00, 3.6000e+01, 0.0000e+00)
    node!(model,  5, 2.8800e+02, 3.6000e+01, 0.0000e+00)
    node!(model,  6, 5.0400e+02, 3.6000e+01, 0.0000e+00)
    node!(model,  7, 0.0000e+00, 7.2000e+01, 0.0000e+00)
    node!(model,  8, 2.8800e+02, 7.2000e+01, 0.0000e+00)
    node!(model,  9, 5.0400e+02, 7.2000e+01, 0.0000e+00)
    node!(model, 10, 0.0000e+00, 1.0800e+02, 0.0000e+00)
    node!(model, 11, 2.8800e+02, 1.0800e+02, 0.0000e+00)
    node!(model, 12, 5.0400e+02, 1.0800e+02, 0.0000e+00)
    node!(model, 13, 0.0000e+00, 1.4400e+02, 0.0000e+00)
    node!(model, 14, 2.8800e+02, 1.4400e+02, 0.0000e+00)
    node!(model, 15, 3.4200e+02, 1.4400e+02, 0.0000e+00)
    node!(model, 16, 3.9600e+02, 1.4400e+02, 0.0000e+00)
    node!(model, 17, 4.5000e+02, 1.4400e+02, 0.0000e+00)
    node!(model, 18, 5.0400e+02, 1.4400e+02, 0.0000e+00)
    node!(model, 19, 0.0000e+00, 1.6800e+02, 0.0000e+00)
    node!(model, 20, 2.8800e+02, 1.6800e+02, 0.0000e+00)
    node!(model, 21, 0.0000e+00, 1.9200e+02, 0.0000e+00)
    node!(model, 22, 2.8800e+02, 1.9200e+02, 0.0000e+00)
    node!(model, 23, 0.0000e+00, 2.1600e+02, 0.0000e+00)
    node!(model, 24, 2.8800e+02, 2.1600e+02, 0.0000e+00)
    node!(model, 25, 0.0000e+00, 2.4000e+02, 0.0000e+00)
    node!(model, 26, 7.2000e+01, 2.4000e+02, 0.0000e+00)
    node!(model, 27, 1.4400e+02, 2.4000e+02, 0.0000e+00)
    node!(model, 28, 2.1600e+02, 2.4000e+02, 0.0000e+00)
    node!(model, 29, 2.8800e+02, 2.4000e+02, 0.0000e+00)

    # Define the sections:
    section!(model, 1, 1.4700e+01, 8.0000e+02, 4.0100e+01, 1.2400e+00)
    section!(model, 2, 1.7100e+01, 2.2800e+02, 7.5100e+01, 3.3300e+00)
    section!(model, 3, 8.8500e+00, 2.9100e+02, 1.9600e+01, 3.8000e-01)
    section!(model, 4, 1.1700e+01, 1.4600e+02, 4.9100e+01, 1.1200e+00)

    # Define the materials:
    material!(model, 1, 29000, 0.3, 0.290)

    # Define the elements:
    element!(model,  1,  1,  4, 4, 1)
    element!(model,  2,  4,  7, 4, 1)
    element!(model,  3,  7, 10, 4, 1)
    element!(model,  4, 10, 13, 4, 1)
    element!(model,  5,  2,  5, 2, 1)
    element!(model,  6,  5,  8, 2, 1)
    element!(model,  7,  8, 11, 2, 1)
    element!(model,  8, 11, 14, 2, 1)
    element!(model,  9,  3,  6, 4, 1)
    element!(model, 10,  6,  9, 4, 1)
    element!(model, 11,  9, 12, 4, 1)
    element!(model, 12, 12, 18, 4, 1)
    element!(model, 13, 14, 15, 3, 1)
    element!(model, 14, 15, 16, 3, 1)
    element!(model, 15, 16, 17, 3, 1)
    element!(model, 16, 17, 18, 3, 1)
    element!(model, 17, 13, 19, 4, 1)
    element!(model, 18, 19, 21, 4, 1)
    element!(model, 19, 21, 23, 4, 1)
    element!(model, 20, 23, 25, 4, 1)
    element!(model, 21, 14, 20, 2, 1)
    element!(model, 22, 20, 22, 2, 1)
    element!(model, 23, 22, 24, 2, 1)
    element!(model, 24, 24, 29, 2, 1)
    element!(model, 25, 25, 26, 1, 1)
    element!(model, 26, 26, 27, 1, 1)
    element!(model, 27, 27, 28, 1, 1)
    element!(model, 28, 28, 29, 1, 1)

    # Define the loads:
    concload!(model, 14, -2.3362e+00, -1.4773e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 15,  0.0000e+00, -2.9546e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 16,  0.0000e+00, -2.9546e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 17,  0.0000e+00, -2.9546e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 18, -3.5042e+00, -1.4773e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 25,  0.0000e+00, -1.9697e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 26,  0.0000e+00, -3.9394e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 27,  0.0000e+00, -3.9394e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 28,  0.0000e+00, -3.9394e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)
    concload!(model, 29, -2.3362e+00, -1.9697e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00)

    # Solve the model using a linear elastic analysis:
    solve!(model, LinearElasticAnalysis())

    # Check the nodal displacements:
    @test getnodaldisplacements(model,  1) ≈ [ 0.0000e+00,  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.4588e-02] rtol = 1E-3
    @test getnodaldisplacements(model,  2) ≈ [ 0.0000e+00,  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  8.6070e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  3) ≈ [ 0.0000e+00,  0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.3185e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  4) ≈ [-5.1719e-01, -8.2521e-03, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.3923e-02] rtol = 1E-3
    @test getnodaldisplacements(model,  5) ≈ [-3.0066e-01, -1.0652e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  7.8412e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  6) ≈ [-5.4774e-02, -5.4379e-03, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.9275e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  7) ≈ [-9.8654e-01, -1.6504e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.1930e-02] rtol = 1E-3
    @test getnodaldisplacements(model,  8) ≈ [-5.4619e-01, -2.1304e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  5.5439e-03] rtol = 1E-3
    @test getnodaldisplacements(model,  9) ≈ [-1.5340e-01, -1.0876e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  3.7545e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 10) ≈ [-1.3602e+00, -2.4756e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  8.6073e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 11) ≈ [-6.8144e-01, -3.1956e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.7149e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 12) ≈ [-3.3972e-01, -1.6314e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  6.7996e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 13) ≈ [-1.5903e+00, -3.3008e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  3.9558e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 14) ≈ [-6.5128e-01, -4.2609e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -3.6457e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 15) ≈ [-6.5286e-01, -4.6883e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00, -9.1366e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 16) ≈ [-6.5443e-01, -7.9684e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00, -1.7094e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 17) ≈ [-6.5601e-01, -6.0469e-01, 0.0000e+00, 0.0000e+00, 0.0000e+00,  8.4267e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 18) ≈ [-6.5758e-01, -2.1751e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.1063e-02] rtol = 1E-3
    @test getnodaldisplacements(model, 19) ≈ [-1.6404e+00, -3.8510e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.1642e-04] rtol = 1E-3
    @test getnodaldisplacements(model, 20) ≈ [-6.1419e-01, -4.6471e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  5.8395e-04] rtol = 1E-3
    @test getnodaldisplacements(model, 21) ≈ [-1.5912e+00, -4.4011e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -4.3136e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 22) ≈ [-6.8071e-01, -5.0333e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  4.9883e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 23) ≈ [-1.4286e+00, -4.9513e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -9.3344e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 24) ≈ [-8.5503e-01, -5.4195e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  9.5673e-03] rtol = 1E-3
    @test getnodaldisplacements(model, 25) ≈ [-1.1384e+00, -5.5014e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00, -1.4946e-02] rtol = 1E-3
    @test getnodaldisplacements(model, 26) ≈ [-1.1391e+00, -1.0918e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, -1.1691e-02] rtol = 1E-3
    @test getnodaldisplacements(model, 27) ≈ [-1.1399e+00, -1.5327e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.4031e-04] rtol = 1E-3
    @test getnodaldisplacements(model, 28) ≈ [-1.1406e+00, -1.0770e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.1745e-02] rtol = 1E-3
    @test getnodaldisplacements(model, 29) ≈ [-1.1413e+00, -5.8057e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,  1.4321e-02] rtol = 1E-3

    # Check the nodal reactions:
    # @test getnodalreactions(model, solution, 1) ≈ [ 4.3419e+00, 7.7776e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00] rtol = 1E-3
    # @test getnodalreactions(model, solution, 2) ≈ [ 7.8139e+00, 1.4673e+02, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00] rtol = 1E-3
    # @test getnodalreactions(model, solution, 3) ≈ [-3.9792e+00, 5.1252e+01, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00] rtol = 1E-3
end
