// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Product.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "expression.h"

#include "expressions.h"

#include "functions.h"

#include "i.h"

#include "identities.h"

#include "integrationconstant.h"

#include "power.h"

#include "product.h"

#include "sum.h"

#include "symbol.h"

#include "symbols.h"

#include <QDebug>

#include <boost/assert.hpp>

#include "product.h"

SireCAS::Product __copy__(const SireCAS::Product &other){ return SireCAS::Product(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Product_class(){

    { //::SireCAS::Product
        typedef bp::class_< SireCAS::Product, bp::bases< SireCAS::ExBase > > Product_exposer_t;
        Product_exposer_t Product_exposer = Product_exposer_t( "Product", "\nThis class holds a collection of expressions that are to be multiplied (or divided)\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty, zero Product") );
        bp::scope Product_scope( Product_exposer );
        Product_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("ex0"), bp::arg("ex1") ), "Construct the product of two expressions") );
        Product_exposer.def( bp::init< SireCAS::Expressions const & >(( bp::arg("expressions") ), "Construct the product of all of the expressions in expressions") );
        Product_exposer.def( bp::init< SireCAS::Product const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Product::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::Product::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::Product::children );
            
            Product_exposer.def( 
                "children"
                , children_function_value
                , bp::release_gil_policy()
                , "Return the child expressions of this product" );
        
        }
        { //::SireCAS::Product::conjugate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*conjugate_function_type)(  ) const;
            conjugate_function_type conjugate_function_value( &::SireCAS::Product::conjugate );
            
            Product_exposer.def( 
                "conjugate"
                , conjugate_function_value
                , bp::release_gil_policy()
                , "Return the complex conjugate of this product" );
        
        }
        { //::SireCAS::Product::denominator
        
            typedef ::SireCAS::Product ( ::SireCAS::Product::*denominator_function_type)(  ) const;
            denominator_function_type denominator_function_value( &::SireCAS::Product::denominator );
            
            Product_exposer.def( 
                "denominator"
                , denominator_function_value
                , bp::release_gil_policy()
                , "Return the Product of expressions on the denominator of this Product" );
        
        }
        { //::SireCAS::Product::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::Product::differentiate );
            
            Product_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the differential of this product..." );
        
        }
        { //::SireCAS::Product::evaluate
        
            typedef double ( ::SireCAS::Product::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Product::evaluate );
            
            Product_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Evaluate this product" );
        
        }
        { //::SireCAS::Product::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Product::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Product::evaluate );
            
            Product_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Evaluate this product" );
        
        }
        { //::SireCAS::Product::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::Product::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::Product::expand );
            
            Product_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Product::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::Product::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::Product::functions );
            
            Product_exposer.def( 
                "functions"
                , functions_function_value
                , bp::release_gil_policy()
                , "Return all of the functions used in this product" );
        
        }
        { //::SireCAS::Product::hash
        
            typedef ::uint ( ::SireCAS::Product::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireCAS::Product::hash );
            
            Product_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "Return a hash for this product" );
        
        }
        { //::SireCAS::Product::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*integrate_function_type)( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::Product::integrate );
            
            Product_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Need to write integral of a product..." );
        
        }
        { //::SireCAS::Product::isComplex
        
            typedef bool ( ::SireCAS::Product::*isComplex_function_type)(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::Product::isComplex );
            
            Product_exposer.def( 
                "isComplex"
                , isComplex_function_value
                , bp::release_gil_policy()
                , "Return whether or not this Product contains any complex terms" );
        
        }
        { //::SireCAS::Product::isCompound
        
            typedef bool ( ::SireCAS::Product::*isCompound_function_type)(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::Product::isCompound );
            
            Product_exposer.def( 
                "isCompound"
                , isCompound_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is compound (requires brackets when printing)" );
        
        }
        { //::SireCAS::Product::isConstant
        
            typedef bool ( ::SireCAS::Product::*isConstant_function_type)(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::Product::isConstant );
            
            Product_exposer.def( 
                "isConstant"
                , isConstant_function_value
                , bp::release_gil_policy()
                , "Return whether this is a constant" );
        
        }
        { //::SireCAS::Product::isFunction
        
            typedef bool ( ::SireCAS::Product::*isFunction_function_type)( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::Product::isFunction );
            
            Product_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("arg0") )
                , bp::release_gil_policy()
                , "Return whether or not this is a function of symbol" );
        
        }
        { //::SireCAS::Product::numerator
        
            typedef ::SireCAS::Product ( ::SireCAS::Product::*numerator_function_type)(  ) const;
            numerator_function_type numerator_function_value( &::SireCAS::Product::numerator );
            
            Product_exposer.def( 
                "numerator"
                , numerator_function_value
                , bp::release_gil_policy()
                , "Return the Product of expressions on the numerator of this Product" );
        
        }
        Product_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Product::reduce
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*reduce_function_type)(  ) const;
            reduce_function_type reduce_function_value( &::SireCAS::Product::reduce );
            
            Product_exposer.def( 
                "reduce"
                , reduce_function_value
                , bp::release_gil_policy()
                , "Reduce a Product down to a simple form. This will not attempt to collapse\ncommon factors - if you want to do this then call the collapse function." );
        
        }
        { //::SireCAS::Product::series
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*series_function_type)( ::SireCAS::Symbol const &,int ) const;
            series_function_type series_function_value( &::SireCAS::Product::series );
            
            Product_exposer.def( 
                "series"
                , series_function_value
                , ( bp::arg("symbol"), bp::arg("n") )
                , bp::release_gil_policy()
                , "Return the series expansion of this product with respect to symbol, to order n" );
        
        }
        { //::SireCAS::Product::simplify
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*simplify_function_type)( int ) const;
            simplify_function_type simplify_function_value( &::SireCAS::Product::simplify );
            
            Product_exposer.def( 
                "simplify"
                , simplify_function_value
                , ( bp::arg("options")=(int)(0) )
                , "Try to simplify this product" );
        
        }
        { //::SireCAS::Product::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::Product::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::Product::substitute );
            
            Product_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , bp::release_gil_policy()
                , "Return the product with the identities in identities substituted in" );
        
        }
        { //::SireCAS::Product::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::Product::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Product::symbols );
            
            Product_exposer.def( 
                "symbols"
                , symbols_function_value
                , bp::release_gil_policy()
                , "Return all of the symbols used in this product" );
        
        }
        { //::SireCAS::Product::toOpenMMString
        
            typedef ::QString ( ::SireCAS::Product::*toOpenMMString_function_type)(  ) const;
            toOpenMMString_function_type toOpenMMString_function_value( &::SireCAS::Product::toOpenMMString );
            
            Product_exposer.def( 
                "toOpenMMString"
                , toOpenMMString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this Product in the OpenMM syntax" );
        
        }
        { //::SireCAS::Product::toString
        
            typedef ::QString ( ::SireCAS::Product::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Product::toString );
            
            Product_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this Product" );
        
        }
        { //::SireCAS::Product::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Product::typeName );
            
            Product_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Product::what
        
            typedef char const * ( ::SireCAS::Product::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Product::what );
            
            Product_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Product_exposer.staticmethod( "typeName" );
        Product_exposer.def( "__copy__", &__copy__);
        Product_exposer.def( "__deepcopy__", &__copy__);
        Product_exposer.def( "clone", &__copy__);
        Product_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Product >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Product_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Product >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Product_exposer.def_pickle(sire_pickle_suite< ::SireCAS::Product >());
        Product_exposer.def( "__str__", &__str__< ::SireCAS::Product > );
        Product_exposer.def( "__repr__", &__str__< ::SireCAS::Product > );
        Product_exposer.def( "__hash__", &::SireCAS::Product::hash );
    }

}
