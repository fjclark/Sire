// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChainsWithAtoms.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "withatoms.h"

#include "withatoms.h"

SireMol::ChainsWithAtoms __copy__(const SireMol::ChainsWithAtoms &other){ return SireMol::ChainsWithAtoms(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChainsWithAtoms_class(){

    { //::SireMol::ChainsWithAtoms
        typedef bp::class_< SireMol::ChainsWithAtoms, bp::bases< SireMol::ChainID, SireID::ID > > ChainsWithAtoms_exposer_t;
        ChainsWithAtoms_exposer_t ChainsWithAtoms_exposer = ChainsWithAtoms_exposer_t( "ChainsWithAtoms", bp::init< >() );
        bp::scope ChainsWithAtoms_scope( ChainsWithAtoms_exposer );
        ChainsWithAtoms_exposer.def( bp::init< SireMol::AtomID const & >(( bp::arg("atomid") )) );
        ChainsWithAtoms_exposer.def( bp::init< SireMol::ChainsWithAtoms const & >(( bp::arg("other") )) );
        { //::SireMol::ChainsWithAtoms::atomID
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ChainsWithAtoms::*atomID_function_type)(  ) const;
            atomID_function_type atomID_function_value( &::SireMol::ChainsWithAtoms::atomID );
            
            ChainsWithAtoms_exposer.def( 
                "atomID"
                , atomID_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::ChainsWithAtoms::hash
        
            typedef ::uint ( ::SireMol::ChainsWithAtoms::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ChainsWithAtoms::hash );
            
            ChainsWithAtoms_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::ChainsWithAtoms::isNull
        
            typedef bool ( ::SireMol::ChainsWithAtoms::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ChainsWithAtoms::isNull );
            
            ChainsWithAtoms_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::ChainsWithAtoms::map
        
            typedef ::QList< SireMol::ChainIdx > ( ::SireMol::ChainsWithAtoms::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ChainsWithAtoms::map );
            
            ChainsWithAtoms_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        ChainsWithAtoms_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainsWithAtoms::operator=
        
            typedef ::SireMol::ChainsWithAtoms & ( ::SireMol::ChainsWithAtoms::*assign_function_type)( ::SireMol::ChainsWithAtoms const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainsWithAtoms::operator= );
            
            ChainsWithAtoms_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        ChainsWithAtoms_exposer.def( bp::self == bp::other< SireID::ID >() );
        ChainsWithAtoms_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainsWithAtoms::toString
        
            typedef ::QString ( ::SireMol::ChainsWithAtoms::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainsWithAtoms::toString );
            
            ChainsWithAtoms_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::ChainsWithAtoms::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainsWithAtoms::typeName );
            
            ChainsWithAtoms_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::ChainsWithAtoms::what
        
            typedef char const * ( ::SireMol::ChainsWithAtoms::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ChainsWithAtoms::what );
            
            ChainsWithAtoms_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ChainsWithAtoms_exposer.staticmethod( "typeName" );
        ChainsWithAtoms_exposer.def( "__copy__", &__copy__);
        ChainsWithAtoms_exposer.def( "__deepcopy__", &__copy__);
        ChainsWithAtoms_exposer.def( "clone", &__copy__);
        ChainsWithAtoms_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainsWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainsWithAtoms_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainsWithAtoms >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainsWithAtoms_exposer.def( "__str__", &__str__< ::SireMol::ChainsWithAtoms > );
        ChainsWithAtoms_exposer.def( "__repr__", &__str__< ::SireMol::ChainsWithAtoms > );
    }

}
