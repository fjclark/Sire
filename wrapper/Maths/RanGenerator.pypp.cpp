// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "RanGenerator.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "rangenerator.h"

#include "vector.h"

#include <QDebug>

#include <QMutex>

#include <QUuid>

#include <QVector>

#include <boost/noncopyable.hpp>

#include <boost/scoped_array.hpp>

#include <limits>

#include "rangenerator.h"

SireMaths::RanGenerator __copy__(const SireMaths::RanGenerator &other){ return SireMaths::RanGenerator(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireMaths::RanGenerator&){ return "SireMaths::RanGenerator";}

void register_RanGenerator_class(){

    { //::SireMaths::RanGenerator
        typedef bp::class_< SireMaths::RanGenerator > RanGenerator_exposer_t;
        RanGenerator_exposer_t RanGenerator_exposer = RanGenerator_exposer_t( "RanGenerator", bp::init< >() );
        bp::scope RanGenerator_scope( RanGenerator_exposer );
        RanGenerator_exposer.def( bp::init< quint32 >(( bp::arg("seed") )) );
        RanGenerator_exposer.def( bp::init< QVector< unsigned int > const & >(( bp::arg("seed") )) );
        RanGenerator_exposer.def( bp::init< SireMaths::RanGenerator const & >(( bp::arg("other") )) );
        { //::SireMaths::RanGenerator::getState
        
            typedef ::QVector< unsigned int > ( ::SireMaths::RanGenerator::*getState_function_type)(  ) const;
            getState_function_type getState_function_value( &::SireMaths::RanGenerator::getState );
            
            RanGenerator_exposer.def( 
                "getState"
                , getState_function_value );
        
        }
        { //::SireMaths::RanGenerator::global
        
            typedef ::SireMaths::RanGenerator const & ( *global_function_type )(  );
            global_function_type global_function_value( &::SireMaths::RanGenerator::global );
            
            RanGenerator_exposer.def( 
                "global"
                , global_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        RanGenerator_exposer.def( bp::self != bp::self );
        { //::SireMaths::RanGenerator::operator=
        
            typedef ::SireMaths::RanGenerator & ( ::SireMaths::RanGenerator::*assign_function_type)( ::SireMaths::RanGenerator const & ) ;
            assign_function_type assign_function_value( &::SireMaths::RanGenerator::operator= );
            
            RanGenerator_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        RanGenerator_exposer.def( bp::self == bp::self );
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)(  ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value );
        
        }
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)( double ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value
                , ( bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::rand
        
            typedef double ( ::SireMaths::RanGenerator::*rand_function_type)( double,double ) const;
            rand_function_type rand_function_value( &::SireMaths::RanGenerator::rand );
            
            RanGenerator_exposer.def( 
                "rand"
                , rand_function_value
                , ( bp::arg("minval"), bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)(  ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)( double ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value
                , ( bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::rand53
        
            typedef double ( ::SireMaths::RanGenerator::*rand53_function_type)( double,double ) const;
            rand53_function_type rand53_function_value( &::SireMaths::RanGenerator::rand53 );
            
            RanGenerator_exposer.def( 
                "rand53"
                , rand53_function_value
                , ( bp::arg("minval"), bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::randBool
        
            typedef bool ( ::SireMaths::RanGenerator::*randBool_function_type)(  ) const;
            randBool_function_type randBool_function_value( &::SireMaths::RanGenerator::randBool );
            
            RanGenerator_exposer.def( 
                "randBool"
                , randBool_function_value );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::quint32 ( ::SireMaths::RanGenerator::*randInt_function_type)(  ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::quint32 ( ::SireMaths::RanGenerator::*randInt_function_type)( ::quint32 ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value
                , ( bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::randInt
        
            typedef ::qint32 ( ::SireMaths::RanGenerator::*randInt_function_type)( ::qint32,::qint32 ) const;
            randInt_function_type randInt_function_value( &::SireMaths::RanGenerator::randInt );
            
            RanGenerator_exposer.def( 
                "randInt"
                , randInt_function_value
                , ( bp::arg("minval"), bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::quint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)(  ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::quint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)( ::quint64 ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value
                , ( bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::randInt64
        
            typedef ::qint64 ( ::SireMaths::RanGenerator::*randInt64_function_type)( ::qint64,::qint64 ) const;
            randInt64_function_type randInt64_function_value( &::SireMaths::RanGenerator::randInt64 );
            
            RanGenerator_exposer.def( 
                "randInt64"
                , randInt64_function_value
                , ( bp::arg("minval"), bp::arg("maxval") ) );
        
        }
        { //::SireMaths::RanGenerator::randNorm
        
            typedef double ( ::SireMaths::RanGenerator::*randNorm_function_type)( double,double ) const;
            randNorm_function_type randNorm_function_value( &::SireMaths::RanGenerator::randNorm );
            
            RanGenerator_exposer.def( 
                "randNorm"
                , randNorm_function_value
                , ( bp::arg("mean"), bp::arg("variance") ) );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)(  ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::quint32 ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("seed") ) );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::QVector< unsigned int > const & ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("seed") ) );
        
        }
        { //::SireMaths::RanGenerator::seed
        
            typedef void ( ::SireMaths::RanGenerator::*seed_function_type)( ::SireMaths::RanGenerator const & ) ;
            seed_function_type seed_function_value( &::SireMaths::RanGenerator::seed );
            
            RanGenerator_exposer.def( 
                "seed"
                , seed_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )(  );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::quint32 );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("seed") ) );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::QVector< unsigned int > const & );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("seed") ) );
        
        }
        { //::SireMaths::RanGenerator::seedGlobal
        
            typedef void ( *seedGlobal_function_type )( ::SireMaths::RanGenerator const & );
            seedGlobal_function_type seedGlobal_function_value( &::SireMaths::RanGenerator::seedGlobal );
            
            RanGenerator_exposer.def( 
                "seedGlobal"
                , seedGlobal_function_value
                , ( bp::arg("other") ) );
        
        }
        { //::SireMaths::RanGenerator::setState
        
            typedef void ( ::SireMaths::RanGenerator::*setState_function_type)( ::QVector< unsigned int > const & ) ;
            setState_function_type setState_function_value( &::SireMaths::RanGenerator::setState );
            
            RanGenerator_exposer.def( 
                "setState"
                , setState_function_value
                , ( bp::arg("state") ) );
        
        }
        { //::SireMaths::RanGenerator::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMaths::RanGenerator::typeName );
            
            RanGenerator_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMaths::RanGenerator::vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*vectorOnSphere_function_type)(  ) const;
            vectorOnSphere_function_type vectorOnSphere_function_value( &::SireMaths::RanGenerator::vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "vectorOnSphere"
                , vectorOnSphere_function_value );
        
        }
        { //::SireMaths::RanGenerator::vectorOnSphere
        
            typedef ::SireMaths::Vector ( ::SireMaths::RanGenerator::*vectorOnSphere_function_type)( double ) const;
            vectorOnSphere_function_type vectorOnSphere_function_value( &::SireMaths::RanGenerator::vectorOnSphere );
            
            RanGenerator_exposer.def( 
                "vectorOnSphere"
                , vectorOnSphere_function_value
                , ( bp::arg("radius") ) );
        
        }
        { //::SireMaths::RanGenerator::what
        
            typedef char const * ( ::SireMaths::RanGenerator::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMaths::RanGenerator::what );
            
            RanGenerator_exposer.def( 
                "what"
                , what_function_value );
        
        }
        RanGenerator_exposer.staticmethod( "global" );
        RanGenerator_exposer.staticmethod( "seedGlobal" );
        RanGenerator_exposer.staticmethod( "typeName" );
        RanGenerator_exposer.def( "__copy__", &__copy__);
        RanGenerator_exposer.def( "__deepcopy__", &__copy__);
        RanGenerator_exposer.def( "clone", &__copy__);
        RanGenerator_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMaths::RanGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RanGenerator_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMaths::RanGenerator >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        RanGenerator_exposer.def( "__str__", &pvt_get_name);
        RanGenerator_exposer.def( "__repr__", &pvt_get_name);
    }

}
