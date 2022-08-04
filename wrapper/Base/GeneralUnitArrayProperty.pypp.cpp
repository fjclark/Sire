// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GeneralUnitArrayProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/generalunitproperty.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "generalunitproperty.h"

#include "generalunitproperty.h"

SireBase::GeneralUnitArrayProperty __copy__(const SireBase::GeneralUnitArrayProperty &other){ return SireBase::GeneralUnitArrayProperty(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_GeneralUnitArrayProperty_class(){

    { //::SireBase::GeneralUnitArrayProperty
        typedef bp::class_< SireBase::GeneralUnitArrayProperty, bp::bases< SireBase::Property > > GeneralUnitArrayProperty_exposer_t;
        GeneralUnitArrayProperty_exposer_t GeneralUnitArrayProperty_exposer = GeneralUnitArrayProperty_exposer_t( "GeneralUnitArrayProperty", "", bp::init< >("") );
        bp::scope GeneralUnitArrayProperty_scope( GeneralUnitArrayProperty_exposer );
        GeneralUnitArrayProperty_exposer.def( bp::init< QVector< SireUnits::Dimension::GeneralUnit > const & >(( bp::arg("units") ), "") );
        GeneralUnitArrayProperty_exposer.def( bp::init< QList< SireUnits::Dimension::GeneralUnit > const & >(( bp::arg("units") ), "") );
        GeneralUnitArrayProperty_exposer.def( bp::init< SireBase::GeneralUnitArrayProperty const & >(( bp::arg("other") ), "") );
        GeneralUnitArrayProperty_exposer.def( bp::self != bp::self );
        { //::SireBase::GeneralUnitArrayProperty::operator=
        
            typedef ::SireBase::GeneralUnitArrayProperty & ( ::SireBase::GeneralUnitArrayProperty::*assign_function_type)( ::SireBase::GeneralUnitArrayProperty const & ) ;
            assign_function_type assign_function_value( &::SireBase::GeneralUnitArrayProperty::operator= );
            
            GeneralUnitArrayProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GeneralUnitArrayProperty_exposer.def( bp::self == bp::self );
        { //::SireBase::GeneralUnitArrayProperty::toString
        
            typedef ::QString ( ::SireBase::GeneralUnitArrayProperty::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::GeneralUnitArrayProperty::toString );
            
            GeneralUnitArrayProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::GeneralUnitArrayProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::GeneralUnitArrayProperty::typeName );
            
            GeneralUnitArrayProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::GeneralUnitArrayProperty::what
        
            typedef char const * ( ::SireBase::GeneralUnitArrayProperty::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireBase::GeneralUnitArrayProperty::what );
            
            GeneralUnitArrayProperty_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        GeneralUnitArrayProperty_exposer.staticmethod( "typeName" );
        GeneralUnitArrayProperty_exposer.def( "__copy__", &__copy__);
        GeneralUnitArrayProperty_exposer.def( "__deepcopy__", &__copy__);
        GeneralUnitArrayProperty_exposer.def( "clone", &__copy__);
        GeneralUnitArrayProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::GeneralUnitArrayProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GeneralUnitArrayProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::GeneralUnitArrayProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GeneralUnitArrayProperty_exposer.def_pickle(sire_pickle_suite< ::SireBase::GeneralUnitArrayProperty >());
        GeneralUnitArrayProperty_exposer.def( "__str__", &__str__< ::SireBase::GeneralUnitArrayProperty > );
        GeneralUnitArrayProperty_exposer.def( "__repr__", &__str__< ::SireBase::GeneralUnitArrayProperty > );
        GeneralUnitArrayProperty_exposer.def( "__len__", &__len_size< ::SireBase::GeneralUnitArrayProperty > );
        GeneralUnitArrayProperty_exposer.def( "__getitem__", &::SireBase::GeneralUnitArrayProperty::getitem );
    }

}
