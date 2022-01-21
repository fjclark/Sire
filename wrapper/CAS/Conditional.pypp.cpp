// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Conditional.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "complexvalues.h"

#include "conditional.h"

#include "expressions.h"

#include "functions.h"

#include "identities.h"

#include "symbols.h"

#include "values.h"

#include "conditional.h"

SireCAS::Conditional __copy__(const SireCAS::Conditional &other){ return SireCAS::Conditional(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Conditional_class(){

    { //::SireCAS::Conditional
        typedef bp::class_< SireCAS::Conditional, bp::bases< SireCAS::ExBase > > Conditional_exposer_t;
        Conditional_exposer_t Conditional_exposer = Conditional_exposer_t( "Conditional", "This is a conditional expression. If the condition is true,\nthen true_expression is evaluated, else if the condition\nis false then false_expression is evaluate\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Conditional_scope( Conditional_exposer );
        Conditional_exposer.def( bp::init< SireCAS::Condition const &, SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("condition"), bp::arg("true_expression"), bp::arg("false_expression") ), "Construct a conditional where if condition is true, then\ntrue_expression is evaluated, while if condition is false,\nthen false_expression is evaluated") );
        Conditional_exposer.def( bp::init< SireCAS::Conditional const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Conditional::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::Conditional::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::Conditional::children );
            
            Conditional_exposer.def( 
                "children"
                , children_function_value
                , "Return the children of this expression" );
        
        }
        { //::SireCAS::Conditional::condition
        
            typedef ::SireCAS::Condition const & ( ::SireCAS::Conditional::*condition_function_type)(  ) const;
            condition_function_type condition_function_value( &::SireCAS::Conditional::condition );
            
            Conditional_exposer.def( 
                "condition"
                , condition_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the condition" );
        
        }
        { //::SireCAS::Conditional::conjugate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*conjugate_function_type)(  ) const;
            conjugate_function_type conjugate_function_value( &::SireCAS::Conditional::conjugate );
            
            Conditional_exposer.def( 
                "conjugate"
                , conjugate_function_value
                , "Return the complex conjugate of this expression" );
        
        }
        { //::SireCAS::Conditional::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::Conditional::differentiate );
            
            Conditional_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , "Return the differential of this expression" );
        
        }
        { //::SireCAS::Conditional::evaluate
        
            typedef double ( ::SireCAS::Conditional::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Conditional::evaluate );
            
            Conditional_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this expression for the passed values" );
        
        }
        { //::SireCAS::Conditional::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Conditional::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Conditional::evaluate );
            
            Conditional_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this expresion for the passed values" );
        
        }
        { //::SireCAS::Conditional::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::Conditional::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::Conditional::expand );
            
            Conditional_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , "Expand this expression in terms of symbol" );
        
        }
        { //::SireCAS::Conditional::falseExpression
        
            typedef ::SireCAS::Expression const & ( ::SireCAS::Conditional::*falseExpression_function_type)(  ) const;
            falseExpression_function_type falseExpression_function_value( &::SireCAS::Conditional::falseExpression );
            
            Conditional_exposer.def( 
                "falseExpression"
                , falseExpression_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the expression to be evaluated if the condition is false" );
        
        }
        { //::SireCAS::Conditional::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::Conditional::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::Conditional::functions );
            
            Conditional_exposer.def( 
                "functions"
                , functions_function_value
                , "Return the functions used in this expression" );
        
        }
        { //::SireCAS::Conditional::hash
        
            typedef ::uint ( ::SireCAS::Conditional::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireCAS::Conditional::hash );
            
            Conditional_exposer.def( 
                "hash"
                , hash_function_value
                , "Hash this conditional" );
        
        }
        { //::SireCAS::Conditional::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*integrate_function_type)( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::Conditional::integrate );
            
            Conditional_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") )
                , "Return the integral of this expression" );
        
        }
        { //::SireCAS::Conditional::isComplex
        
            typedef bool ( ::SireCAS::Conditional::*isComplex_function_type)(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::Conditional::isComplex );
            
            Conditional_exposer.def( 
                "isComplex"
                , isComplex_function_value
                , "Is this a complex expression?" );
        
        }
        { //::SireCAS::Conditional::isCompound
        
            typedef bool ( ::SireCAS::Conditional::*isCompound_function_type)(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::Conditional::isCompound );
            
            Conditional_exposer.def( 
                "isCompound"
                , isCompound_function_value
                , "Is this a compound expression?" );
        
        }
        { //::SireCAS::Conditional::isConstant
        
            typedef bool ( ::SireCAS::Conditional::*isConstant_function_type)(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::Conditional::isConstant );
            
            Conditional_exposer.def( 
                "isConstant"
                , isConstant_function_value
                , "Return whether or not this is constant" );
        
        }
        { //::SireCAS::Conditional::isFunction
        
            typedef bool ( ::SireCAS::Conditional::*isFunction_function_type)( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::Conditional::isFunction );
            
            Conditional_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("arg0") )
                , "Return whether or not this is a function of symbol" );
        
        }
        { //::SireCAS::Conditional::isNull
        
            typedef bool ( ::SireCAS::Conditional::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireCAS::Conditional::isNull );
            
            Conditional_exposer.def( 
                "isNull"
                , isNull_function_value
                , "Return whether or not this is null" );
        
        }
        { //::SireCAS::Conditional::operator=
        
            typedef ::SireCAS::Conditional & ( ::SireCAS::Conditional::*assign_function_type)( ::SireCAS::Conditional const & ) ;
            assign_function_type assign_function_value( &::SireCAS::Conditional::operator= );
            
            Conditional_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Conditional_exposer.def( bp::self == bp::self );
        Conditional_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Conditional::series
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*series_function_type)( ::SireCAS::Symbol const &,int ) const;
            series_function_type series_function_value( &::SireCAS::Conditional::series );
            
            Conditional_exposer.def( 
                "series"
                , series_function_value
                , ( bp::arg("symbol"), bp::arg("n") )
                , "Return the series expansion of this product with respect to symbol, to order n" );
        
        }
        { //::SireCAS::Conditional::simplify
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*simplify_function_type)( int ) const;
            simplify_function_type simplify_function_value( &::SireCAS::Conditional::simplify );
            
            Conditional_exposer.def( 
                "simplify"
                , simplify_function_value
                , ( bp::arg("options")=(int)(0) )
                , "Try to simplify this condition" );
        
        }
        { //::SireCAS::Conditional::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::Conditional::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::Conditional::substitute );
            
            Conditional_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , "Substitute identities into this expression" );
        
        }
        { //::SireCAS::Conditional::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::Conditional::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Conditional::symbols );
            
            Conditional_exposer.def( 
                "symbols"
                , symbols_function_value
                , "Return the symbols used in this expression" );
        
        }
        { //::SireCAS::Conditional::toString
        
            typedef ::QString ( ::SireCAS::Conditional::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Conditional::toString );
            
            Conditional_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this conditional" );
        
        }
        { //::SireCAS::Conditional::trueExpression
        
            typedef ::SireCAS::Expression const & ( ::SireCAS::Conditional::*trueExpression_function_type)(  ) const;
            trueExpression_function_type trueExpression_function_value( &::SireCAS::Conditional::trueExpression );
            
            Conditional_exposer.def( 
                "trueExpression"
                , trueExpression_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the expression to be evaluated if the condition is true" );
        
        }
        { //::SireCAS::Conditional::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Conditional::typeName );
            
            Conditional_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireCAS::Conditional::what
        
            typedef char const * ( ::SireCAS::Conditional::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Conditional::what );
            
            Conditional_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        Conditional_exposer.staticmethod( "typeName" );
        Conditional_exposer.def( "__copy__", &__copy__);
        Conditional_exposer.def( "__deepcopy__", &__copy__);
        Conditional_exposer.def( "clone", &__copy__);
        Conditional_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Conditional >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Conditional_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Conditional >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Conditional_exposer.def( "__setstate__", &__setstate__base64< ::SireCAS::Conditional > );
        Conditional_exposer.def( "__getstate__", &__getstate__base64< ::SireCAS::Conditional > );
        Conditional_exposer.def( "__str__", &__str__< ::SireCAS::Conditional > );
        Conditional_exposer.def( "__repr__", &__str__< ::SireCAS::Conditional > );
        Conditional_exposer.def( "__hash__", &::SireCAS::Conditional::hash );
    }

}
