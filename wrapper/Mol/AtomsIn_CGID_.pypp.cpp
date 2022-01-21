// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomsIn_CGID_.pypp.hpp"

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

SireMol::AtomsIn<SireMol::CGID> __copy__(const SireMol::AtomsIn<SireMol::CGID> &other){ return SireMol::AtomsIn<SireMol::CGID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomsIn_CGID__class(){

    { //::SireMol::AtomsIn< SireMol::CGID >
        typedef bp::class_< SireMol::AtomsIn< SireMol::CGID >, bp::bases< SireMol::AtomID, SireID::ID > > AtomsIn_CGID__exposer_t;
        AtomsIn_CGID__exposer_t AtomsIn_CGID__exposer = AtomsIn_CGID__exposer_t( "AtomsIn_CGID_", "", bp::init< >("") );
        bp::scope AtomsIn_CGID__scope( AtomsIn_CGID__exposer );
        AtomsIn_CGID__exposer.def( bp::init< SireMol::CGID const & >(( bp::arg("id") ), "") );
        AtomsIn_CGID__exposer.def( bp::init< SireMol::CGID const &, qint32 >(( bp::arg("id"), bp::arg("i") ), "") );
        AtomsIn_CGID__exposer.def( bp::init< SireMol::CGID const &, qint32, qint32 >(( bp::arg("id"), bp::arg("i"), bp::arg("j") ), "") );
        AtomsIn_CGID__exposer.def( bp::init< SireMol::AtomsIn< SireMol::CGID > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomsIn< SireMol::CGID >::hash
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef ::uint ( ::SireMol::AtomsIn< SireMol::CGID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::AtomsIn< SireMol::CGID >::hash );
            
            AtomsIn_CGID__exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::CGID >::isNull
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef bool ( ::SireMol::AtomsIn< SireMol::CGID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomsIn< SireMol::CGID >::isNull );
            
            AtomsIn_CGID__exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::CGID >::map
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::AtomsIn< SireMol::CGID >::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::AtomsIn< SireMol::CGID >::map );
            
            AtomsIn_CGID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        AtomsIn_CGID__exposer.def( bp::self != bp::self );
        AtomsIn_CGID__exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::CGID >::operator=
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef ::SireMol::AtomsIn< SireMol::CGID > & ( ::SireMol::AtomsIn< SireMol::CGID >::*assign_function_type)( ::SireMol::AtomsIn< SireMol::CGID > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomsIn< SireMol::CGID >::operator= );
            
            AtomsIn_CGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomsIn_CGID__exposer.def( bp::self == bp::self );
        AtomsIn_CGID__exposer.def( bp::self == bp::other< SireID::ID >() );
        { //::SireMol::AtomsIn< SireMol::CGID >::toString
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef ::QString ( ::SireMol::AtomsIn< SireMol::CGID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomsIn< SireMol::CGID >::toString );
            
            AtomsIn_CGID__exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::CGID >::typeName
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomsIn< SireMol::CGID >::typeName );
            
            AtomsIn_CGID__exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::AtomsIn< SireMol::CGID >::what
        
            typedef SireMol::AtomsIn< SireMol::CGID > exported_class_t;
            typedef char const * ( ::SireMol::AtomsIn< SireMol::CGID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::AtomsIn< SireMol::CGID >::what );
            
            AtomsIn_CGID__exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        AtomsIn_CGID__exposer.staticmethod( "typeName" );
        AtomsIn_CGID__exposer.def( "__copy__", &__copy__);
        AtomsIn_CGID__exposer.def( "__deepcopy__", &__copy__);
        AtomsIn_CGID__exposer.def( "clone", &__copy__);
        AtomsIn_CGID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomsIn<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_CGID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomsIn<SireMol::CGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomsIn_CGID__exposer.def( "__getstate_manages_dict__", true);
        AtomsIn_CGID__exposer.def( "__safe_for_unpickling__", true);
        AtomsIn_CGID__exposer.def( "__setstate__", &__setstate__base64< ::SireMol::AtomsIn<SireMol::CGID> > );
        AtomsIn_CGID__exposer.def( "__getstate__", &__getstate__base64< ::SireMol::AtomsIn<SireMol::CGID> > );
        AtomsIn_CGID__exposer.def( "__str__", &__str__< ::SireMol::AtomsIn<SireMol::CGID> > );
        AtomsIn_CGID__exposer.def( "__repr__", &__str__< ::SireMol::AtomsIn<SireMol::CGID> > );
        AtomsIn_CGID__exposer.def( "__hash__", &::SireMol::AtomsIn<SireMol::CGID>::hash );
    }

}
