// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "DistVector.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "distvector.h"

#include "matrix.h"

#include <QRegExp>

#include <QString>

#include <cmath>

#include "distvector.h"

SireMaths::DistVector __copy__(const SireMaths::DistVector &other){ return SireMaths::DistVector(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_DistVector_class(){

    { //::SireMaths::DistVector
        typedef bp::class_< SireMaths::DistVector > DistVector_exposer_t;
        DistVector_exposer_t DistVector_exposer = DistVector_exposer_t( "DistVector", "This is a vector that stores the vector as a unit vector giving\nthe direction, and a scalar giving the magnitude\n\nAuthor: Christopher Woods\n", bp::init< >("Create the zero vector") );
        bp::scope DistVector_scope( DistVector_exposer );
        DistVector_exposer.def( bp::init< SireMaths::Vector const & >(( bp::arg("vec") ), "Create from the passed vector") );
        DistVector_exposer.def( bp::init< SireMaths::DistVector const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMaths::DistVector::angle
        
            typedef ::SireUnits::Dimension::Angle ( *angle_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            angle_function_type angle_function_value( &::SireMaths::DistVector::angle );
            
            DistVector_exposer.def( 
                "angle"
                , angle_function_value
                , ( bp::arg("v0"), bp::arg("v1") )
                , bp::release_gil_policy()
                , "Return the angle between vectors v0 and v1 - this is the smallest\nangle, and will always lie between 0 and 180 degrees" );
        
        }
        { //::SireMaths::DistVector::angle
        
            typedef ::SireUnits::Dimension::Angle ( *angle_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            angle_function_type angle_function_value( &::SireMaths::DistVector::angle );
            
            DistVector_exposer.def( 
                "angle"
                , angle_function_value
                , ( bp::arg("v0"), bp::arg("v1"), bp::arg("v2") )
                , bp::release_gil_policy()
                , "Return the angle between v0-v1-v2 (treating the vectors as points in space)" );
        
        }
        { //::SireMaths::DistVector::at
        
            typedef double ( ::SireMaths::DistVector::*at_function_type)( unsigned int ) const;
            at_function_type at_function_value( &::SireMaths::DistVector::at );
            
            DistVector_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Access elements by index" );
        
        }
        { //::SireMaths::DistVector::b
        
            typedef double ( ::SireMaths::DistVector::*b_function_type)(  ) const;
            b_function_type b_function_value( &::SireMaths::DistVector::b );
            
            DistVector_exposer.def( 
                "b"
                , b_function_value
                , bp::release_gil_policy()
                , "Return the components via rgb (limited between 0 and 1)" );
        
        }
        { //::SireMaths::DistVector::bearing
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::DistVector::*bearing_function_type)(  ) const;
            bearing_function_type bearing_function_value( &::SireMaths::DistVector::bearing );
            
            DistVector_exposer.def( 
                "bearing"
                , bearing_function_value
                , bp::release_gil_policy()
                , "Return the bearing of this vector against (0,1,0) (north) on the xy plane" );
        
        }
        { //::SireMaths::DistVector::bearingXY
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::DistVector::*bearingXY_function_type)( ::SireMaths::DistVector const & ) const;
            bearingXY_function_type bearingXY_function_value( &::SireMaths::DistVector::bearingXY );
            
            DistVector_exposer.def( 
                "bearingXY"
                , bearingXY_function_value
                , ( bp::arg("v") )
                , bp::release_gil_policy()
                , "Return the bearing of this vector against v on the xy plane" );
        
        }
        { //::SireMaths::DistVector::bearingXZ
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::DistVector::*bearingXZ_function_type)( ::SireMaths::DistVector const & ) const;
            bearingXZ_function_type bearingXZ_function_value( &::SireMaths::DistVector::bearingXZ );
            
            DistVector_exposer.def( 
                "bearingXZ"
                , bearingXZ_function_value
                , ( bp::arg("v") )
                , bp::release_gil_policy()
                , "Return the bearing of this vector against v on the xz plane" );
        
        }
        { //::SireMaths::DistVector::bearingYZ
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMaths::DistVector::*bearingYZ_function_type)( ::SireMaths::DistVector const & ) const;
            bearingYZ_function_type bearingYZ_function_value( &::SireMaths::DistVector::bearingYZ );
            
            DistVector_exposer.def( 
                "bearingYZ"
                , bearingYZ_function_value
                , ( bp::arg("v") )
                , bp::release_gil_policy()
                , "Return the bearing of this vector against v on the yz plane" );
        
        }
        { //::SireMaths::DistVector::count
        
            typedef unsigned int ( ::SireMaths::DistVector::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMaths::DistVector::count );
            
            DistVector_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the size of the Vector (always 3 - unless you disagree\nwith me that we should be living in a 3-dimensional space)" );
        
        }
        { //::SireMaths::DistVector::cross
        
            typedef ::SireMaths::DistVector ( *cross_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            cross_function_type cross_function_value( &::SireMaths::DistVector::cross );
            
            DistVector_exposer.def( 
                "cross"
                , cross_function_value
                , ( bp::arg("v0"), bp::arg("v1") )
                , bp::release_gil_policy()
                , "Return the cross product of v0 and v1" );
        
        }
        { //::SireMaths::DistVector::dihedral
        
            typedef ::SireUnits::Dimension::Angle ( *dihedral_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const &,::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            dihedral_function_type dihedral_function_value( &::SireMaths::DistVector::dihedral );
            
            DistVector_exposer.def( 
                "dihedral"
                , dihedral_function_value
                , ( bp::arg("v0"), bp::arg("v1"), bp::arg("v2"), bp::arg("v3") )
                , bp::release_gil_policy()
                , "Return the dihedral angle between v0-v1-v2-v3 (treating the vectors as points)" );
        
        }
        { //::SireMaths::DistVector::direction
        
            typedef ::SireMaths::Vector const & ( ::SireMaths::DistVector::*direction_function_type)(  ) const;
            direction_function_type direction_function_value( &::SireMaths::DistVector::direction );
            
            DistVector_exposer.def( 
                "direction"
                , direction_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the direction of this vector" );
        
        }
        { //::SireMaths::DistVector::distance
        
            typedef double ( *distance_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            distance_function_type distance_function_value( &::SireMaths::DistVector::distance );
            
            DistVector_exposer.def( 
                "distance"
                , distance_function_value
                , ( bp::arg("v1"), bp::arg("v2") )
                , bp::release_gil_policy()
                , "Return the distance between two vectors" );
        
        }
        { //::SireMaths::DistVector::distance2
        
            typedef double ( *distance2_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            distance2_function_type distance2_function_value( &::SireMaths::DistVector::distance2 );
            
            DistVector_exposer.def( 
                "distance2"
                , distance2_function_value
                , ( bp::arg("v1"), bp::arg("v2") )
                , bp::release_gil_policy()
                , "Return the distance squared between two vectors" );
        
        }
        { //::SireMaths::DistVector::dot
        
            typedef double ( *dot_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            dot_function_type dot_function_value( &::SireMaths::DistVector::dot );
            
            DistVector_exposer.def( 
                "dot"
                , dot_function_value
                , ( bp::arg("v0"), bp::arg("v1") )
                , bp::release_gil_policy()
                , "Return the dot product of v0 and v1" );
        
        }
        { //::SireMaths::DistVector::fromString
        
            typedef ::SireMaths::DistVector ( *fromString_function_type )( ::QString const & );
            fromString_function_type fromString_function_value( &::SireMaths::DistVector::fromString );
            
            DistVector_exposer.def( 
                "fromString"
                , fromString_function_value
                , ( bp::arg("str") )
                , bp::release_gil_policy()
                , "Construct a DistVector from the QString representation returned by toString()\nThrow: SireError::invalid_arg\n" );
        
        }
        { //::SireMaths::DistVector::g
        
            typedef double ( ::SireMaths::DistVector::*g_function_type)(  ) const;
            g_function_type g_function_value( &::SireMaths::DistVector::g );
            
            DistVector_exposer.def( 
                "g"
                , g_function_value
                , bp::release_gil_policy()
                , "Return the components via rgb (limited between 0 and 1)" );
        
        }
        { //::SireMaths::DistVector::generate
        
            typedef ::SireMaths::DistVector ( *generate_function_type )( double,::SireMaths::DistVector const &,::SireUnits::Dimension::Angle const &,::SireMaths::DistVector const &,::SireUnits::Dimension::Angle const &,::SireMaths::DistVector const & );
            generate_function_type generate_function_value( &::SireMaths::DistVector::generate );
            
            DistVector_exposer.def( 
                "generate"
                , generate_function_value
                , ( bp::arg("dst"), bp::arg("v1"), bp::arg("ang"), bp::arg("v2"), bp::arg("dih"), bp::arg("v3") )
                , bp::release_gil_policy()
                , "Generate a vector, v0, that has distance dst v0-v1, angle ang v0-v1-v2,\nand dihedral dih v0-v1-v2-v3" );
        
        }
        { //::SireMaths::DistVector::getitem
        
            typedef double ( ::SireMaths::DistVector::*getitem_function_type)( int ) const;
            getitem_function_type getitem_function_value( &::SireMaths::DistVector::getitem );
            
            DistVector_exposer.def( 
                "getitem"
                , getitem_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::DistVector::invDistance
        
            typedef double ( *invDistance_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            invDistance_function_type invDistance_function_value( &::SireMaths::DistVector::invDistance );
            
            DistVector_exposer.def( 
                "invDistance"
                , invDistance_function_value
                , ( bp::arg("v1"), bp::arg("v2") )
                , bp::release_gil_policy()
                , "Return the 1  distance between two vectors" );
        
        }
        { //::SireMaths::DistVector::invDistance2
        
            typedef double ( *invDistance2_function_type )( ::SireMaths::DistVector const &,::SireMaths::DistVector const & );
            invDistance2_function_type invDistance2_function_value( &::SireMaths::DistVector::invDistance2 );
            
            DistVector_exposer.def( 
                "invDistance2"
                , invDistance2_function_value
                , ( bp::arg("v1"), bp::arg("v2") )
                , bp::release_gil_policy()
                , "Return 1  distance2 between two vectors" );
        
        }
        { //::SireMaths::DistVector::invLength
        
            typedef double ( ::SireMaths::DistVector::*invLength_function_type)(  ) const;
            invLength_function_type invLength_function_value( &::SireMaths::DistVector::invLength );
            
            DistVector_exposer.def( 
                "invLength"
                , invLength_function_value
                , bp::release_gil_policy()
                , "Return the inverse of the length of the vector" );
        
        }
        { //::SireMaths::DistVector::invLength2
        
            typedef double ( ::SireMaths::DistVector::*invLength2_function_type)(  ) const;
            invLength2_function_type invLength2_function_value( &::SireMaths::DistVector::invLength2 );
            
            DistVector_exposer.def( 
                "invLength2"
                , invLength2_function_value
                , bp::release_gil_policy()
                , "Return the inverse length squared" );
        
        }
        { //::SireMaths::DistVector::isZero
        
            typedef bool ( ::SireMaths::DistVector::*isZero_function_type)(  ) const;
            isZero_function_type isZero_function_value( &::SireMaths::DistVector::isZero );
            
            DistVector_exposer.def( 
                "isZero"
                , isZero_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a zero length vector" );
        
        }
        { //::SireMaths::DistVector::length
        
            typedef double ( ::SireMaths::DistVector::*length_function_type)(  ) const;
            length_function_type length_function_value( &::SireMaths::DistVector::length );
            
            DistVector_exposer.def( 
                "length"
                , length_function_value
                , bp::release_gil_policy()
                , "Return the length of the vector" );
        
        }
        { //::SireMaths::DistVector::length2
        
            typedef double ( ::SireMaths::DistVector::*length2_function_type)(  ) const;
            length2_function_type length2_function_value( &::SireMaths::DistVector::length2 );
            
            DistVector_exposer.def( 
                "length2"
                , length2_function_value
                , bp::release_gil_policy()
                , "Return the length^2 of the vector" );
        
        }
        { //::SireMaths::DistVector::magnitude
        
            typedef double ( ::SireMaths::DistVector::*magnitude_function_type)(  ) const;
            magnitude_function_type magnitude_function_value( &::SireMaths::DistVector::magnitude );
            
            DistVector_exposer.def( 
                "magnitude"
                , magnitude_function_value
                , bp::release_gil_policy()
                , "Return the magnitude of this vector" );
        
        }
        { //::SireMaths::DistVector::manhattanLength
        
            typedef double ( ::SireMaths::DistVector::*manhattanLength_function_type)(  ) const;
            manhattanLength_function_type manhattanLength_function_value( &::SireMaths::DistVector::manhattanLength );
            
            DistVector_exposer.def( 
                "manhattanLength"
                , manhattanLength_function_value
                , bp::release_gil_policy()
                , "Return the manhattan length of the vector" );
        
        }
        { //::SireMaths::DistVector::max
        
            typedef ::SireMaths::DistVector ( ::SireMaths::DistVector::*max_function_type)( ::SireMaths::DistVector const & ) const;
            max_function_type max_function_value( &::SireMaths::DistVector::max );
            
            DistVector_exposer.def( 
                "max"
                , max_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return a vector that has the maximum xyz components out of this\nand other" );
        
        }
        { //::SireMaths::DistVector::metricTensor
        
            typedef ::SireMaths::Matrix ( ::SireMaths::DistVector::*metricTensor_function_type)(  ) const;
            metricTensor_function_type metricTensor_function_value( &::SireMaths::DistVector::metricTensor );
            
            DistVector_exposer.def( 
                "metricTensor"
                , metricTensor_function_value
                , bp::release_gil_policy()
                , "Return the metric tensor of a vector, i.e.\n| yy + zz,    -xy    -xz      |\n|    -yx,   xx + zz  -yz      |\n|    -zx       -zy    xx + yy |\n" );
        
        }
        { //::SireMaths::DistVector::min
        
            typedef ::SireMaths::DistVector ( ::SireMaths::DistVector::*min_function_type)( ::SireMaths::DistVector const & ) const;
            min_function_type min_function_value( &::SireMaths::DistVector::min );
            
            DistVector_exposer.def( 
                "min"
                , min_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return a vector that has the minimum components" );
        
        }
        { //::SireMaths::DistVector::normalise
        
            typedef ::SireMaths::DistVector ( ::SireMaths::DistVector::*normalise_function_type)(  ) const;
            normalise_function_type normalise_function_value( &::SireMaths::DistVector::normalise );
            
            DistVector_exposer.def( 
                "normalise"
                , normalise_function_value
                , bp::release_gil_policy()
                , "Return a normalised form of the vector" );
        
        }
        DistVector_exposer.def( bp::self != bp::self );
        DistVector_exposer.def( bp::self *= bp::other< double >() );
        DistVector_exposer.def( bp::self += bp::self );
        DistVector_exposer.def( -bp::self );
        DistVector_exposer.def( bp::self -= bp::self );
        DistVector_exposer.def( bp::self /= bp::other< double >() );
        { //::SireMaths::DistVector::operator=
        
            typedef ::SireMaths::DistVector const & ( ::SireMaths::DistVector::*assign_function_type)( ::SireMaths::DistVector const & ) ;
            assign_function_type assign_function_value( &::SireMaths::DistVector::operator= );
            
            DistVector_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        DistVector_exposer.def( bp::self == bp::self );
        { //::SireMaths::DistVector::operator[]
        
            typedef double ( ::SireMaths::DistVector::*__getitem___function_type)( unsigned int ) const;
            __getitem___function_type __getitem___function_value( &::SireMaths::DistVector::operator[] );
            
            DistVector_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMaths::DistVector::r
        
            typedef double ( ::SireMaths::DistVector::*r_function_type)(  ) const;
            r_function_type r_function_value( &::SireMaths::DistVector::r );
            
            DistVector_exposer.def( 
                "r"
                , r_function_value
                , bp::release_gil_policy()
                , "Return the components via rgb (limited between 0 and 1)" );
        
        }
        { //::SireMaths::DistVector::setMax
        
            typedef void ( ::SireMaths::DistVector::*setMax_function_type)( ::SireMaths::DistVector const & ) ;
            setMax_function_type setMax_function_value( &::SireMaths::DistVector::setMax );
            
            DistVector_exposer.def( 
                "setMax"
                , setMax_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Set this Vector so that it has the maximum xyz components out of\nthis and other (e.g. this->x = max(this->x(),other.x() etc.)" );
        
        }
        { //::SireMaths::DistVector::setMin
        
            typedef void ( ::SireMaths::DistVector::*setMin_function_type)( ::SireMaths::DistVector const & ) ;
            setMin_function_type setMin_function_value( &::SireMaths::DistVector::setMin );
            
            DistVector_exposer.def( 
                "setMin"
                , setMin_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Set this Vector so that it has the minimum xyz components" );
        
        }
        { //::SireMaths::DistVector::toString
        
            typedef ::QString ( ::SireMaths::DistVector::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMaths::DistVector::toString );
            
            DistVector_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a QString representation of the vector" );
        
        }
        { //::SireMaths::DistVector::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::DistVector::typeName );
            
            DistVector_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::DistVector::what
        
            typedef char const * ( ::SireMaths::DistVector::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::DistVector::what );
            
            DistVector_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMaths::DistVector::x
        
            typedef double ( ::SireMaths::DistVector::*x_function_type)(  ) const;
            x_function_type x_function_value( &::SireMaths::DistVector::x );
            
            DistVector_exposer.def( 
                "x"
                , x_function_value
                , bp::release_gil_policy()
                , "Return the x component of the vector" );
        
        }
        { //::SireMaths::DistVector::y
        
            typedef double ( ::SireMaths::DistVector::*y_function_type)(  ) const;
            y_function_type y_function_value( &::SireMaths::DistVector::y );
            
            DistVector_exposer.def( 
                "y"
                , y_function_value
                , bp::release_gil_policy()
                , "Return the y component of the vector" );
        
        }
        { //::SireMaths::DistVector::z
        
            typedef double ( ::SireMaths::DistVector::*z_function_type)(  ) const;
            z_function_type z_function_value( &::SireMaths::DistVector::z );
            
            DistVector_exposer.def( 
                "z"
                , z_function_value
                , bp::release_gil_policy()
                , "Return the z component of the vector" );
        
        }
        DistVector_exposer.staticmethod( "angle" );
        DistVector_exposer.staticmethod( "cross" );
        DistVector_exposer.staticmethod( "dihedral" );
        DistVector_exposer.staticmethod( "distance" );
        DistVector_exposer.staticmethod( "distance2" );
        DistVector_exposer.staticmethod( "dot" );
        DistVector_exposer.staticmethod( "fromString" );
        DistVector_exposer.staticmethod( "generate" );
        DistVector_exposer.staticmethod( "invDistance" );
        DistVector_exposer.staticmethod( "invDistance2" );
        DistVector_exposer.staticmethod( "typeName" );
        DistVector_exposer.def( "__copy__", &__copy__);
        DistVector_exposer.def( "__deepcopy__", &__copy__);
        DistVector_exposer.def( "clone", &__copy__);
        DistVector_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::DistVector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DistVector_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::DistVector >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        DistVector_exposer.def_pickle(sire_pickle_suite< ::SireMaths::DistVector >());
        DistVector_exposer.def( "__str__", &__str__< ::SireMaths::DistVector > );
        DistVector_exposer.def( "__repr__", &__str__< ::SireMaths::DistVector > );
        DistVector_exposer.def( "__len__", &__len_count< ::SireMaths::DistVector > );
        DistVector_exposer.def( "__getitem__", &::SireMaths::DistVector::getitem );
    }

}
