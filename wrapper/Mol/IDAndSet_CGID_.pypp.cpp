// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IDAndSet_CGID_.pypp.hpp"

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

SireID::IDAndSet<SireMol::CGID> __copy__(const SireID::IDAndSet<SireMol::CGID> &other){ return SireID::IDAndSet<SireMol::CGID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_IDAndSet_CGID__class(){

    { //::SireID::IDAndSet< SireMol::CGID >
        typedef bp::class_< SireID::IDAndSet< SireMol::CGID >, bp::bases< SireMol::CGID, SireID::ID > > IDAndSet_CGID__exposer_t;
        IDAndSet_CGID__exposer_t IDAndSet_CGID__exposer = IDAndSet_CGID__exposer_t( "IDAndSet_CGID_", "", bp::init< >("") );
        bp::scope IDAndSet_CGID__scope( IDAndSet_CGID__exposer );
        IDAndSet_CGID__exposer.def( bp::init< SireMol::CGID const & >(( bp::arg("id") ), "") );
        IDAndSet_CGID__exposer.def( bp::init< SireMol::CGID const &, SireMol::CGID const & >(( bp::arg("id0"), bp::arg("id1") ), "") );
        IDAndSet_CGID__exposer.def( bp::init< QList< SireMol::CGIdentifier > const & >(( bp::arg("ids") ), "") );
        IDAndSet_CGID__exposer.def( bp::init< SireID::IDAndSet< SireMol::CGID > const & >(( bp::arg("ids") ), "") );
        IDAndSet_CGID__exposer.def( bp::init< SireID::IDAndSet< SireMol::CGID > const & >(( bp::arg("other") ), "") );
        { //::SireID::IDAndSet< SireMol::CGID >::IDs
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::QSet< SireMol::CGIdentifier > const & ( ::SireID::IDAndSet< SireMol::CGID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDAndSet< SireMol::CGID >::IDs );
            
            IDAndSet_CGID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::hash
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::uint ( ::SireID::IDAndSet< SireMol::CGID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDAndSet< SireMol::CGID >::hash );
            
            IDAndSet_CGID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::isNull
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef bool ( ::SireID::IDAndSet< SireMol::CGID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDAndSet< SireMol::CGID >::isNull );
            
            IDAndSet_CGID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::map
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::QList< SireMol::CGIdx > ( ::SireID::IDAndSet< SireMol::CGID >::*map_function_type)( ::SireID::IDAndSet< SireMol::CGID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireMol::CGID >::map );
            
            IDAndSet_CGID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , bp::release_gil_policy()
                , "" );
        
        }
        IDAndSet_CGID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDAndSet_CGID__exposer.def( bp::self != bp::self );
        IDAndSet_CGID__exposer.def( bp::self != bp::other< SireMol::CGID >() );
        { //::SireID::IDAndSet< SireMol::CGID >::operator=
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::CGID > & ( ::SireID::IDAndSet< SireMol::CGID >::*assign_function_type)( ::SireID::IDAndSet< SireMol::CGID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::CGID >::operator= );
            
            IDAndSet_CGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::operator=
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::CGID > & ( ::SireID::IDAndSet< SireMol::CGID >::*assign_function_type)( ::SireMol::CGID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::CGID >::operator= );
            
            IDAndSet_CGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IDAndSet_CGID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDAndSet_CGID__exposer.def( bp::self == bp::self );
        IDAndSet_CGID__exposer.def( bp::self == bp::other< SireMol::CGID >() );
        { //::SireID::IDAndSet< SireMol::CGID >::toString
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef ::QString ( ::SireID::IDAndSet< SireMol::CGID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDAndSet< SireMol::CGID >::toString );
            
            IDAndSet_CGID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::typeName
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDAndSet< SireMol::CGID >::typeName );
            
            IDAndSet_CGID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireMol::CGID >::what
        
            typedef SireID::IDAndSet< SireMol::CGID > exported_class_t;
            typedef char const * ( ::SireID::IDAndSet< SireMol::CGID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDAndSet< SireMol::CGID >::what );
            
            IDAndSet_CGID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IDAndSet_CGID__exposer.staticmethod( "typeName" );
        IDAndSet_CGID__exposer.def( "__copy__", &__copy__);
        IDAndSet_CGID__exposer.def( "__deepcopy__", &__copy__);
        IDAndSet_CGID__exposer.def( "clone", &__copy__);
        IDAndSet_CGID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDAndSet<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_CGID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDAndSet<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_CGID__exposer.def_pickle(sire_pickle_suite< ::SireID::IDAndSet<SireMol::CGID> >());
        IDAndSet_CGID__exposer.def( "__str__", &__str__< ::SireID::IDAndSet<SireMol::CGID> > );
        IDAndSet_CGID__exposer.def( "__repr__", &__str__< ::SireID::IDAndSet<SireMol::CGID> > );
        IDAndSet_CGID__exposer.def( "__hash__", &::SireID::IDAndSet<SireMol::CGID>::hash );
    }

}
