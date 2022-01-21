// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ExpressionProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireCAS/expressionproperty.h"

#include "SireCAS/values.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "expressionproperty.h"

#include "expressionproperty.h"

SireCAS::ExpressionProperty __copy__(const SireCAS::ExpressionProperty &other){ return SireCAS::ExpressionProperty(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ExpressionProperty_class(){

    { //::SireCAS::ExpressionProperty
        typedef bp::class_< SireCAS::ExpressionProperty, bp::bases< SireBase::Property > > ExpressionProperty_exposer_t;
        ExpressionProperty_exposer_t ExpressionProperty_exposer = ExpressionProperty_exposer_t( "ExpressionProperty", "This class provides a thin Property wrapper around SireCAS objects\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor - this constructs an empty expression") );
        bp::scope ExpressionProperty_scope( ExpressionProperty_exposer );
        ExpressionProperty_exposer.def( bp::init< SireCAS::ExBase const & >(( bp::arg("exbase") ), "Construct from the passed expression") );
        ExpressionProperty_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("expression") ), "Construct from the passed expression") );
        ExpressionProperty_exposer.def( bp::init< SireCAS::ExpressionProperty const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::ExpressionProperty::asABoolean
        
            typedef bool ( ::SireCAS::ExpressionProperty::*asABoolean_function_type)(  ) const;
            asABoolean_function_type asABoolean_function_value( &::SireCAS::ExpressionProperty::asABoolean );
            
            ExpressionProperty_exposer.def( 
                "asABoolean"
                , asABoolean_function_value
                , "Return this number cast as a boolean" );
        
        }
        { //::SireCAS::ExpressionProperty::asADouble
        
            typedef double ( ::SireCAS::ExpressionProperty::*asADouble_function_type)(  ) const;
            asADouble_function_type asADouble_function_value( &::SireCAS::ExpressionProperty::asADouble );
            
            ExpressionProperty_exposer.def( 
                "asADouble"
                , asADouble_function_value
                , "Return this number cast as a double" );
        
        }
        { //::SireCAS::ExpressionProperty::asAnInteger
        
            typedef int ( ::SireCAS::ExpressionProperty::*asAnInteger_function_type)(  ) const;
            asAnInteger_function_type asAnInteger_function_value( &::SireCAS::ExpressionProperty::asAnInteger );
            
            ExpressionProperty_exposer.def( 
                "asAnInteger"
                , asAnInteger_function_value
                , "Return this number cast as an integer" );
        
        }
        { //::SireCAS::ExpressionProperty::isABoolean
        
            typedef bool ( ::SireCAS::ExpressionProperty::*isABoolean_function_type)(  ) const;
            isABoolean_function_type isABoolean_function_value( &::SireCAS::ExpressionProperty::isABoolean );
            
            ExpressionProperty_exposer.def( 
                "isABoolean"
                , isABoolean_function_value
                , "" );
        
        }
        { //::SireCAS::ExpressionProperty::isADouble
        
            typedef bool ( ::SireCAS::ExpressionProperty::*isADouble_function_type)(  ) const;
            isADouble_function_type isADouble_function_value( &::SireCAS::ExpressionProperty::isADouble );
            
            ExpressionProperty_exposer.def( 
                "isADouble"
                , isADouble_function_value
                , "" );
        
        }
        { //::SireCAS::ExpressionProperty::isAnInteger
        
            typedef bool ( ::SireCAS::ExpressionProperty::*isAnInteger_function_type)(  ) const;
            isAnInteger_function_type isAnInteger_function_value( &::SireCAS::ExpressionProperty::isAnInteger );
            
            ExpressionProperty_exposer.def( 
                "isAnInteger"
                , isAnInteger_function_value
                , "" );
        
        }
        ExpressionProperty_exposer.def( bp::self != bp::self );
        { //::SireCAS::ExpressionProperty::operator=
        
            typedef ::SireCAS::ExpressionProperty & ( ::SireCAS::ExpressionProperty::*assign_function_type)( ::SireCAS::ExpressionProperty const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ExpressionProperty::operator= );
            
            ExpressionProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ExpressionProperty_exposer.def( bp::self == bp::self );
        { //::SireCAS::ExpressionProperty::toString
        
            typedef ::QString ( ::SireCAS::ExpressionProperty::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::ExpressionProperty::toString );
            
            ExpressionProperty_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireCAS::ExpressionProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::ExpressionProperty::typeName );
            
            ExpressionProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireCAS::ExpressionProperty::value
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExpressionProperty::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireCAS::ExpressionProperty::value );
            
            ExpressionProperty_exposer.def( 
                "value"
                , value_function_value
                , "Return this number cast as a double" );
        
        }
        ExpressionProperty_exposer.staticmethod( "typeName" );
        ExpressionProperty_exposer.def( "__copy__", &__copy__);
        ExpressionProperty_exposer.def( "__deepcopy__", &__copy__);
        ExpressionProperty_exposer.def( "clone", &__copy__);
        ExpressionProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::ExpressionProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExpressionProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::ExpressionProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExpressionProperty_exposer.def_pickle(sire_pickle_suite< ::SireCAS::ExpressionProperty >());
        ExpressionProperty_exposer.def( "__str__", &__str__< ::SireCAS::ExpressionProperty > );
        ExpressionProperty_exposer.def( "__repr__", &__str__< ::SireCAS::ExpressionProperty > );
    }

}
