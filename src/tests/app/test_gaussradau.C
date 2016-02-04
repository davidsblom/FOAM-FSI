
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "GaussRadau.H"

TEST( GaussRadauTest, object )
{
    const size_t num_nodes = 2;
    fsi::quadrature::GaussRadau<double> gaussRadau( num_nodes );
}

TEST( GaussRadauTest, nodes )
{
    const size_t num_nodes = 5;
    fsi::quadrature::GaussRadau<double> gaussRadau( num_nodes );
    const vector<double> nodes = gaussRadau.get_nodes();

    double tol = 1.0e-13;
    ASSERT_NEAR( nodes[0], 0.05710419611451768219312119255411562123507794559875016432780829293093467820207316458611381681984273686, tol );
    ASSERT_NEAR( nodes[1], 0.2768430136381238276800459976856251411108891696950304683494420488311213396837080367725415285640511309, tol );
    ASSERT_NEAR( nodes[2], 0.5835904323689168200566976686629172486934326398967716401762938418317475019618310120056322774674562993, tol );
    ASSERT_NEAR( nodes[3], 0.8602401356562194478479129188751197667383780225872255049242335941839742579301655644134901549264276107, tol );
    ASSERT_NEAR( nodes[4], 1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000, tol );
}

TEST( GaussRadauTest, api )
{
    const size_t num_nodes = 5;
    fsi::quadrature::GaussRadau<double> gaussRadau( num_nodes );
    gaussRadau.get_nodes();
    gaussRadau.get_q_mat();
    gaussRadau.get_s_mat();
    gaussRadau.get_b_mat();
    gaussRadau.get_q_vec();
    ASSERT_EQ( gaussRadau.get_num_nodes(), num_nodes );
}
