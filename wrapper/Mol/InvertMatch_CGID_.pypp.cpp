// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "InvertMatch_CGID_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "cgid.h"

#include "cgidentifier.h"

#include "chain.h"

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

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "cgid.h"

SireID::InvertMatch<SireMol::CGID> __copy__(const SireID::InvertMatch<SireMol::CGID> &other){ return SireID::InvertMatch<SireMol::CGID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_InvertMatch_CGID__class(){

    { //::SireID::InvertMatch< SireMol::CGID >
        typedef bp::class_< SireID::InvertMatch< SireMol::CGID >, bp::bases< SireMol::CGID, SireID::ID > > InvertMatch_CGID__exposer_t;
        InvertMatch_CGID__exposer_t InvertMatch_CGID__exposer = InvertMatch_CGID__exposer_t( "InvertMatch_CGID_", "", bp::init< >("") );
        bp::scope InvertMatch_CGID__scope( InvertMatch_CGID__exposer );
        InvertMatch_CGID__exposer.def( bp::init< SireMol::CGID const & >(( bp::arg("id") ), "") );
        InvertMatch_CGID__exposer.def( bp::init< SireID::InvertMatch< SireMol::CGID > const & >(( bp::arg("other") ), "") );
        { //::SireID::InvertMatch< SireMol::CGID >::hash
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef ::uint ( ::SireID::InvertMatch< SireMol::CGID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::InvertMatch< SireMol::CGID >::hash );
            
            InvertMatch_CGID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::InvertMatch< SireMol::CGID >::isNull
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef bool ( ::SireID::InvertMatch< SireMol::CGID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::InvertMatch< SireMol::CGID >::isNull );
            
            InvertMatch_CGID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::InvertMatch< SireMol::CGID >::map
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef ::QList< SireMol::CGIdx > ( ::SireID::InvertMatch< SireMol::CGID >::*map_function_type)( ::SireMol::CGID::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::InvertMatch< SireMol::CGID >::map );
            
            InvertMatch_CGID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , bp::release_gil_policy()
                , "" );
        
        }
        InvertMatch_CGID__exposer.def( bp::self != bp::self );
        InvertMatch_CGID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireID::InvertMatch< SireMol::CGID >::operator=
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef ::SireID::InvertMatch< SireMol::CGID > & ( ::SireID::InvertMatch< SireMol::CGID >::*assign_function_type)( ::SireID::InvertMatch< SireMol::CGID > const & ) ;
            assign_function_type assign_function_value( &::SireID::InvertMatch< SireMol::CGID >::operator= );
            
            InvertMatch_CGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        InvertMatch_CGID__exposer.def( bp::self == bp::self );
        InvertMatch_CGID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireID::InvertMatch< SireMol::CGID >::toString
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef ::QString ( ::SireID::InvertMatch< SireMol::CGID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::InvertMatch< SireMol::CGID >::toString );
            
            InvertMatch_CGID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::InvertMatch< SireMol::CGID >::typeName
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::InvertMatch< SireMol::CGID >::typeName );
            
            InvertMatch_CGID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::InvertMatch< SireMol::CGID >::what
        
            typedef SireID::InvertMatch< SireMol::CGID > exported_class_t;
            typedef char const * ( ::SireID::InvertMatch< SireMol::CGID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::InvertMatch< SireMol::CGID >::what );
            
            InvertMatch_CGID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        InvertMatch_CGID__exposer.staticmethod( "typeName" );
        InvertMatch_CGID__exposer.def( "__copy__", &__copy__);
        InvertMatch_CGID__exposer.def( "__deepcopy__", &__copy__);
        InvertMatch_CGID__exposer.def( "clone", &__copy__);
        InvertMatch_CGID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::InvertMatch<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InvertMatch_CGID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::InvertMatch<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        InvertMatch_CGID__exposer.def_pickle(sire_pickle_suite< ::SireID::InvertMatch<SireMol::CGID> >());
        InvertMatch_CGID__exposer.def( "__str__", &__str__< ::SireID::InvertMatch<SireMol::CGID> > );
        InvertMatch_CGID__exposer.def( "__repr__", &__str__< ::SireID::InvertMatch<SireMol::CGID> > );
        InvertMatch_CGID__exposer.def( "__hash__", &::SireID::InvertMatch<SireMol::CGID>::hash );
    }

}
