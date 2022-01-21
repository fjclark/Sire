// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IDAndSet_MonitorID_.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireSystem/errors.h"

#include "monitorid.h"

#include "monitoridx.h"

#include "monitorname.h"

#include "systemmonitors.h"

#include "monitorid.h"

SireID::IDAndSet<SireSystem::MonitorID> __copy__(const SireID::IDAndSet<SireSystem::MonitorID> &other){ return SireID::IDAndSet<SireSystem::MonitorID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_IDAndSet_MonitorID__class(){

    { //::SireID::IDAndSet< SireSystem::MonitorID >
        typedef bp::class_< SireID::IDAndSet< SireSystem::MonitorID >, bp::bases< SireSystem::MonitorID, SireID::ID > > IDAndSet_MonitorID__exposer_t;
        IDAndSet_MonitorID__exposer_t IDAndSet_MonitorID__exposer = IDAndSet_MonitorID__exposer_t( "IDAndSet_MonitorID_", "", bp::init< >("") );
        bp::scope IDAndSet_MonitorID__scope( IDAndSet_MonitorID__exposer );
        IDAndSet_MonitorID__exposer.def( bp::init< SireSystem::MonitorID const & >(( bp::arg("id") ), "") );
        IDAndSet_MonitorID__exposer.def( bp::init< SireSystem::MonitorID const &, SireSystem::MonitorID const & >(( bp::arg("id0"), bp::arg("id1") ), "") );
        IDAndSet_MonitorID__exposer.def( bp::init< QList< SireSystem::MonitorIdentifier > const & >(( bp::arg("ids") ), "") );
        IDAndSet_MonitorID__exposer.def( bp::init< SireID::IDAndSet< SireSystem::MonitorID > const & >(( bp::arg("ids") ), "") );
        IDAndSet_MonitorID__exposer.def( bp::init< SireID::IDAndSet< SireSystem::MonitorID > const & >(( bp::arg("other") ), "") );
        { //::SireID::IDAndSet< SireSystem::MonitorID >::IDs
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::QSet< SireSystem::MonitorIdentifier > const & ( ::SireID::IDAndSet< SireSystem::MonitorID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::IDs );
            
            IDAndSet_MonitorID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::hash
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::uint ( ::SireID::IDAndSet< SireSystem::MonitorID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::hash );
            
            IDAndSet_MonitorID__exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::isNull
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef bool ( ::SireID::IDAndSet< SireSystem::MonitorID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::isNull );
            
            IDAndSet_MonitorID__exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::map
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::QList< SireSystem::MonitorName > ( ::SireID::IDAndSet< SireSystem::MonitorID >::*map_function_type)( ::SireID::IDAndSet< SireSystem::MonitorID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::map );
            
            IDAndSet_MonitorID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , "" );
        
        }
        IDAndSet_MonitorID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDAndSet_MonitorID__exposer.def( bp::self != bp::self );
        IDAndSet_MonitorID__exposer.def( bp::self != bp::other< SireSystem::MonitorID >() );
        { //::SireID::IDAndSet< SireSystem::MonitorID >::operator=
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::SireID::IDAndSet< SireSystem::MonitorID > & ( ::SireID::IDAndSet< SireSystem::MonitorID >::*assign_function_type)( ::SireID::IDAndSet< SireSystem::MonitorID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::operator= );
            
            IDAndSet_MonitorID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::operator=
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::SireID::IDAndSet< SireSystem::MonitorID > & ( ::SireID::IDAndSet< SireSystem::MonitorID >::*assign_function_type)( ::SireSystem::MonitorID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::operator= );
            
            IDAndSet_MonitorID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IDAndSet_MonitorID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDAndSet_MonitorID__exposer.def( bp::self == bp::self );
        IDAndSet_MonitorID__exposer.def( bp::self == bp::other< SireSystem::MonitorID >() );
        { //::SireID::IDAndSet< SireSystem::MonitorID >::toString
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef ::QString ( ::SireID::IDAndSet< SireSystem::MonitorID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::toString );
            
            IDAndSet_MonitorID__exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::typeName
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::typeName );
            
            IDAndSet_MonitorID__exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireID::IDAndSet< SireSystem::MonitorID >::what
        
            typedef SireID::IDAndSet< SireSystem::MonitorID > exported_class_t;
            typedef char const * ( ::SireID::IDAndSet< SireSystem::MonitorID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDAndSet< SireSystem::MonitorID >::what );
            
            IDAndSet_MonitorID__exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        IDAndSet_MonitorID__exposer.staticmethod( "typeName" );
        IDAndSet_MonitorID__exposer.def( "__copy__", &__copy__);
        IDAndSet_MonitorID__exposer.def( "__deepcopy__", &__copy__);
        IDAndSet_MonitorID__exposer.def( "clone", &__copy__);
        IDAndSet_MonitorID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDAndSet<SireSystem::MonitorID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_MonitorID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDAndSet<SireSystem::MonitorID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_MonitorID__exposer.def( "__getstate_manages_dict__", true);
        IDAndSet_MonitorID__exposer.def( "__safe_for_unpickling__", true);
        IDAndSet_MonitorID__exposer.def( "__setstate__", &__setstate__base64< ::SireID::IDAndSet<SireSystem::MonitorID> > );
        IDAndSet_MonitorID__exposer.def( "__getstate__", &__getstate__base64< ::SireID::IDAndSet<SireSystem::MonitorID> > );
        IDAndSet_MonitorID__exposer.def( "__str__", &__str__< ::SireID::IDAndSet<SireSystem::MonitorID> > );
        IDAndSet_MonitorID__exposer.def( "__repr__", &__str__< ::SireID::IDAndSet<SireSystem::MonitorID> > );
        IDAndSet_MonitorID__exposer.def( "__hash__", &::SireID::IDAndSet<SireSystem::MonitorID>::hash );
    }

}
