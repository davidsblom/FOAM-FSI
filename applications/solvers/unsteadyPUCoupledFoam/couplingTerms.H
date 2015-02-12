
{
  // Calculate grad p coupling matrix. Needs to be here if one uses
  // gradient schemes with limiters.  VV, 9/June/2014
  BlockLduSystem<vector, vector> pInU( fvm::grad( p ) );

  // Calculate div U coupling.  Could be calculated only once since
  // it is only geometry dependent.  VV, 9/June/2014
  BlockLduSystem<vector, scalar> UInp( fvm::UDiv( U ) );

  // Last argument in insertBlockCoupling says if the column direction
  // should be incremented. This is needed for arbitrary positioning
  // of U and p in the system. This could be better. VV, 30/April/2014
  UpEqn.insertBlockCoupling( 0, 3, pInU, true );
  UpEqn.insertBlockCoupling( 3, 0, UInp, false );
}
