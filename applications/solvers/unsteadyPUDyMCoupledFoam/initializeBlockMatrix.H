
// Create block matrix
BlockLduMatrix<vector4> A( mesh );

// Block matrix - create x and b
vector4Field & x = Up.internalField();
vector4Field b( x.size(), vector4::zero );

// Set block interfaces properly
A.interfaces() = Up.boundaryField().blockInterfaces();
