using Hephaestus

M = Model()

node!(M, 1, 0, 0, 0)
node!(M, 2, 1, 0, 0)

section!(M, 1, 10, 100, 100, 5)

material!(M, 1, 29000, 0.3, 0.290)

element!(M, 1, 1, 2, 1, 1)

concload!(M, 1, 0, 0, -10, 0, 0, 0)
distload!(M, 1, 0, 0, -10)

generatereport(M)