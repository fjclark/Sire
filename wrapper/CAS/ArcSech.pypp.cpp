// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ArcSech.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "invhyperbolicfuncs.h"

#include "invtrigfuncs.h"

#include "trigfuncs.h"

#include "invhyperbolicfuncs.h"

SireCAS::ArcSech __copy__(const SireCAS::ArcSech &other){ return SireCAS::ArcSech(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ArcSech_class(){

    { //::SireCAS::ArcSech
        typedef bp::class_< SireCAS::ArcSech, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > ArcSech_exposer_t;
        ArcSech_exposer_t ArcSech_exposer = ArcSech_exposer_t( "ArcSech", "Inverse-hyperbolic-secant", bp::init< >("Null constructor") );
        bp::scope ArcSech_scope( ArcSech_exposer );
        ArcSech_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") ), "Construct cos(expression)") );
        ArcSech_exposer.def( bp::init< SireCAS::ArcSech const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::ArcSech::evaluate
        
            typedef double ( ::SireCAS::ArcSech::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcSech::evaluate );
            
            ArcSech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this function" );
        
        }
        { //::SireCAS::ArcSech::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ArcSech::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcSech::evaluate );
            
            ArcSech_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Complex evaluation" );
        
        }
        ArcSech_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::ArcSech::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::ArcSech::typeName );
            
            ArcSech_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireCAS::ArcSech::what
        
            typedef char const * ( ::SireCAS::ArcSech::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::ArcSech::what );
            
            ArcSech_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ArcSech_exposer.staticmethod( "typeName" );
        ArcSech_exposer.def( "__copy__", &__copy__);
        ArcSech_exposer.def( "__deepcopy__", &__copy__);
        ArcSech_exposer.def( "clone", &__copy__);
        ArcSech_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::ArcSech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ArcSech_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::ArcSech >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ArcSech_exposer.def( "__setstate__", &__setstate__base64< ::SireCAS::ArcSech > );
        ArcSech_exposer.def( "__getstate__", &__getstate__base64< ::SireCAS::ArcSech > );
        ArcSech_exposer.def( "__str__", &__str__< ::SireCAS::ArcSech > );
        ArcSech_exposer.def( "__repr__", &__str__< ::SireCAS::ArcSech > );
        ArcSech_exposer.def( "__hash__", &::SireCAS::ArcSech::hash );
    }

}
