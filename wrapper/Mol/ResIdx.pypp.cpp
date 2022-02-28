// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResIdx.pypp.hpp"

namespace bp = boost::python;

#include "residx.h"

#include "residx.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "chain.h"

#include "chainresid.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "groupgroupids.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "resid.h"

#include "residentifier.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "withres.h"

#include "resid.h"

SireMol::ResIdx __copy__(const SireMol::ResIdx &other){ return SireMol::ResIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResIdx_class(){

    { //::SireMol::ResIdx
        typedef bp::class_< SireMol::ResIdx, bp::bases< SireMol::ResID, SireID::ID, SireID::IndexBase > > ResIdx_exposer_t;
        ResIdx_exposer_t ResIdx_exposer = ResIdx_exposer_t( "ResIdx", "This is an ID object that is used to index CutGroups\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope ResIdx_scope( ResIdx_exposer );
        ResIdx_exposer.def( bp::init< quint32 >(( bp::arg("idx") ), "") );
        ResIdx_exposer.def( bp::init< SireMol::ResIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::ResIdx::hash
        
            typedef ::uint ( ::SireMol::ResIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ResIdx::hash );
            
            ResIdx_exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireMol::ResIdx::isNull
        
            typedef bool ( ::SireMol::ResIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ResIdx::isNull );
            
            ResIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireMol::ResIdx::map
        
            typedef ::QList< SireMol::ResIdx > ( ::SireMol::ResIdx::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ResIdx::map );
            
            ResIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::ResIdx::null
        
            typedef ::SireMol::ResIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::ResIdx::null );
            
            ResIdx_exposer.def( 
                "null"
                , null_function_value
                , "" );
        
        }
        { //::SireMol::ResIdx::operator=
        
            typedef ::SireMol::ResIdx & ( ::SireMol::ResIdx::*assign_function_type)( ::SireMol::ResIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResIdx::operator= );
            
            ResIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::ResIdx::toString
        
            typedef ::QString ( ::SireMol::ResIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResIdx::toString );
            
            ResIdx_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::ResIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResIdx::typeName );
            
            ResIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ResIdx::what
        
            typedef char const * ( ::SireMol::ResIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ResIdx::what );
            
            ResIdx_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ResIdx_exposer.staticmethod( "null" );
        ResIdx_exposer.staticmethod( "typeName" );
        ResIdx_exposer.def( "__copy__", &__copy__);
        ResIdx_exposer.def( "__deepcopy__", &__copy__);
        ResIdx_exposer.def( "clone", &__copy__);
        ResIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdx_exposer.def_pickle(sire_pickle_suite< ::SireMol::ResIdx >());
        ResIdx_exposer.def( "__str__", &__str__< ::SireMol::ResIdx > );
        ResIdx_exposer.def( "__repr__", &__str__< ::SireMol::ResIdx > );
        ResIdx_exposer.def( "__hash__", &::SireMol::ResIdx::hash );
    }

}
