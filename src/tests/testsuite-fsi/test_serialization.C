
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <fstream>
#include "gtest/gtest.h"

#ifndef EIGEN_CONFIG_H_
 #define EIGEN_CONFIG_H_

 #include <boost/serialization/array.hpp>
 #define EIGEN_DENSEBASE_PLUGIN "EigenDenseBaseAddons.h"

 #include <Eigen/Core>

#endif // EIGEN_CONFIG_H_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

/// Boost Serialization Helper

template <typename T>
bool serialize(
    const T & data,
    const std::string & filename
    )
{
    std::ofstream ofs( filename.c_str(), std::ios::out );

    if ( !ofs.is_open() )
        return false;

    {
        boost::archive::binary_oarchive oa( ofs );
        oa << data;
    }
    ofs.close();
    return true;
}

template <typename T>
bool deSerialize(
    T & data,
    const std::string & filename
    )
{
    std::ifstream ifs( filename.c_str(), std::ios::in );

    if ( !ifs.is_open() )
        return false;

    {
        boost::archive::binary_iarchive ia( ifs );
        ia >> data;
    }
    ifs.close();
    return true;
}

TEST( serialization, file )
{
    serialize( Eigen::Vector2f( 0.5f, 0.5f ), "a.bin" );
    Eigen::Vector2f a2f;
    deSerialize( a2f, "a.bin" );
    assert( a2f.isApprox( Eigen::Vector2f( 0.5f, 0.5f ) ) );
    Eigen::VectorXf axf;
    deSerialize( axf, "a.bin" );
    assert( axf.isApprox( Eigen::Vector2f( 0.5f, 0.5f ) ) );

    Eigen::Matrix4f m( Eigen::Matrix4f::Random() ), m_in;
    serialize( m.topRows<2>(), "m.bin" );

    Eigen::MatrixXd mat_out( 2343, 3 ), mat_in;
    mat_out.setRandom();
    serialize( mat_out, "mat.bin" );

    deSerialize( mat_in, "mat.bin" );
    assert( mat_in.isApprox( mat_out ) );

    Eigen::MatrixXd mat_out2( 2343, 2 ), mat_in2;
    mat_out.setRandom();
    serialize( mat_out2, "mat2.bin" );

    deSerialize( mat_in2, "mat2.bin" );
    assert( mat_in2.isApprox( mat_out2 ) );
}
