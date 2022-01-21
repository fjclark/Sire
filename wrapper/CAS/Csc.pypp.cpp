// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Csc.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "identities.h"

#include "trigfuncs.h"

#include "trigfuncs.h"

SireCAS::Csc __copy__(const SireCAS::Csc &other){ return SireCAS::Csc(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Csc_class(){

    { //::SireCAS::Csc
        typedef bp::class_< SireCAS::Csc, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > Csc_exposer_t;
        Csc_exposer_t Csc_exposer = Csc_exposer_t( "Csc", "Cosecant", bp::init< >("Null constructor") );
        bp::scope Csc_scope( Csc_exposer );
        Csc_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") ), "Construct cos(expression)") );
        Csc_exposer.def( bp::init< SireCAS::Csc const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Csc::evaluate
        
            typedef double ( ::SireCAS::Csc::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Csc::evaluate );
            
            Csc_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this function" );
        
        }
        { //::SireCAS::Csc::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Csc::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Csc::evaluate );
            
            Csc_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Complex evaluation" );
        
        }
        Csc_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Csc::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Csc::typeName );
            
            Csc_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireCAS::Csc::what
        
            typedef char const * ( ::SireCAS::Csc::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Csc::what );
            
            Csc_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        Csc_exposer.staticmethod( "typeName" );
        Csc_exposer.def( "__copy__", &__copy__);
        Csc_exposer.def( "__deepcopy__", &__copy__);
        Csc_exposer.def( "clone", &__copy__);
        Csc_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Csc >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Csc_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Csc >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Csc_exposer.def( "__setstate__", &__setstate__base64< ::SireCAS::Csc > );
        Csc_exposer.def( "__getstate__", &__getstate__base64< ::SireCAS::Csc > );
        Csc_exposer.def( "__str__", &__str__< ::SireCAS::Csc > );
        Csc_exposer.def( "__repr__", &__str__< ::SireCAS::Csc > );
        Csc_exposer.def( "__hash__", &::SireCAS::Csc::hash );
    }

}
