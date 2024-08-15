using Hephaestus

M = Model()

addnode!(M, 1,  0.0, 0.0, 0.0)
addnode!(M, 2, 10.0, 0.0, 0.0)

addmaterial!(M, 1, 29000, 0.3)

addsection!(M, 1, 100, 1000, 1000, 500)

addelement!(M, 1, 1, 2, 1, 1)