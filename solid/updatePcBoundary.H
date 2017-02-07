forAll(pc.boundaryField(), patchi)
{ 
	pc.boundaryField()[patchi] = pc.boundaryField()[patchi].patchInternalField();
}
