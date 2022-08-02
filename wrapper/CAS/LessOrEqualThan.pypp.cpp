// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "LessOrEqualThan.pypp.hpp"

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

SireCAS::LessOrEqualThan __copy__(const SireCAS::LessOrEqualThan &other){ return SireCAS::LessOrEqualThan(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_LessOrEqualThan_class(){

    { //::SireCAS::LessOrEqualThan
        typedef bp::class_< SireCAS::LessOrEqualThan, bp::bases< SireCAS::Condition, SireCAS::ExBase > > LessOrEqualThan_exposer_t;
        LessOrEqualThan_exposer_t LessOrEqualThan_exposer = LessOrEqualThan_exposer_t( "LessOrEqualThan", "This is a conditional function that returns whether or\nnot the first expression is less than or equal to\nthe second expression\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope LessOrEqualThan_scope( LessOrEqualThan_exposer );
        LessOrEqualThan_exposer.def( bp::init< SireCAS::Expression const &, SireCAS::Expression const & >(( bp::arg("left_hand_side"), bp::arg("right_hand_side") ), "Construct to compare left_hand_side with right_hand_side") );
        LessOrEqualThan_exposer.def( bp::init< SireCAS::LessOrEqualThan const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::LessOrEqualThan::alwaysFalse
        
            typedef bool ( ::SireCAS::LessOrEqualThan::*alwaysFalse_function_type)(  ) const;
            alwaysFalse_function_type alwaysFalse_function_value( &::SireCAS::LessOrEqualThan::alwaysFalse );
            
            LessOrEqualThan_exposer.def( 
                "alwaysFalse"
                , alwaysFalse_function_value
                , bp::release_gil_policy()
                , "Return whether or not we can be absolutely sure that this\ncondition will always be false. Note that this doesnt try\ntoo hard, so some things that are always false may not\nbe reported here as being always false, e.g. x > x + 1" );
        
        }
        { //::SireCAS::LessOrEqualThan::alwaysTrue
        
            typedef bool ( ::SireCAS::LessOrEqualThan::*alwaysTrue_function_type)(  ) const;
            alwaysTrue_function_type alwaysTrue_function_value( &::SireCAS::LessOrEqualThan::alwaysTrue );
            
            LessOrEqualThan_exposer.def( 
                "alwaysTrue"
                , alwaysTrue_function_value
                , bp::release_gil_policy()
                , "Return whether or not we can be absolutely sure that this\ncondition will always be true. Note that this doesnt try\ntoo hard, so some things that are always true may not\nbe reported here as being always true, e.g. x + 1 > x" );
        
        }
        { //::SireCAS::LessOrEqualThan::operator=
        
            typedef ::SireCAS::LessOrEqualThan & ( ::SireCAS::LessOrEqualThan::*assign_function_type)( ::SireCAS::LessOrEqualThan const & ) ;
            assign_function_type assign_function_value( &::SireCAS::LessOrEqualThan::operator= );
            
            LessOrEqualThan_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        LessOrEqualThan_exposer.def( bp::self == bp::self );
        LessOrEqualThan_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::LessOrEqualThan::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::LessOrEqualThan::typeName );
            
            LessOrEqualThan_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::LessOrEqualThan::what
        
            typedef char const * ( ::SireCAS::LessOrEqualThan::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::LessOrEqualThan::what );
            
            LessOrEqualThan_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        LessOrEqualThan_exposer.staticmethod( "typeName" );
        LessOrEqualThan_exposer.def( "__copy__", &__copy__);
        LessOrEqualThan_exposer.def( "__deepcopy__", &__copy__);
        LessOrEqualThan_exposer.def( "clone", &__copy__);
        LessOrEqualThan_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::LessOrEqualThan >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LessOrEqualThan_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::LessOrEqualThan >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        LessOrEqualThan_exposer.def_pickle(sire_pickle_suite< ::SireCAS::LessOrEqualThan >());
        LessOrEqualThan_exposer.def( "__str__", &__str__< ::SireCAS::LessOrEqualThan > );
        LessOrEqualThan_exposer.def( "__repr__", &__str__< ::SireCAS::LessOrEqualThan > );
        LessOrEqualThan_exposer.def( "__hash__", &::SireCAS::LessOrEqualThan::hash );
    }

}
