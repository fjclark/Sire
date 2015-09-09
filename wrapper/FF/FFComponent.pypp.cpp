// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FFComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "ff.h"

#include "ffcomponent.h"

#include <QRegExp>

#include "ffcomponent.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_FFComponent_class(){

    { //::SireFF::FFComponent
        typedef bp::class_< SireFF::FFComponent, bp::bases< SireCAS::Symbol, SireCAS::ExBase >, boost::noncopyable > FFComponent_exposer_t;
        FFComponent_exposer_t FFComponent_exposer = FFComponent_exposer_t( "FFComponent", bp::no_init );
        bp::scope FFComponent_scope( FFComponent_exposer );
        { //::SireFF::FFComponent::componentName
        
            typedef ::QString ( ::SireFF::FFComponent::*componentName_function_type)(  ) const;
            componentName_function_type componentName_function_value( &::SireFF::FFComponent::componentName );
            
            FFComponent_exposer.def( 
                "componentName"
                , componentName_function_value );
        
        }
        { //::SireFF::FFComponent::forceFieldName
        
            typedef ::SireFF::FFName ( ::SireFF::FFComponent::*forceFieldName_function_type)(  ) const;
            forceFieldName_function_type forceFieldName_function_value( &::SireFF::FFComponent::forceFieldName );
            
            FFComponent_exposer.def( 
                "forceFieldName"
                , forceFieldName_function_value );
        
        }
        { //::SireFF::FFComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireFF::FFComponent::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireFF::FFComponent::symbols );
            
            FFComponent_exposer.def( 
                "symbols"
                , symbols_function_value );
        
        }
        { //::SireFF::FFComponent::total
        
            typedef ::SireFF::FFComponent const & ( ::SireFF::FFComponent::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireFF::FFComponent::total );
            
            FFComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireFF::FFComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::FFComponent::typeName );
            
            FFComponent_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireFF::FFComponent::what
        
            typedef char const * ( ::SireFF::FFComponent::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::FFComponent::what );
            
            FFComponent_exposer.def( 
                "what"
                , what_function_value );
        
        }
        FFComponent_exposer.staticmethod( "typeName" );
        FFComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::FFComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::FFComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFComponent_exposer.def( "__str__", &__str__< ::SireFF::FFComponent > );
        FFComponent_exposer.def( "__repr__", &__str__< ::SireFF::FFComponent > );
    }

}
