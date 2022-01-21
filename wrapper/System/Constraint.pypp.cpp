// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Constraint.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireStream/streamdata.hpp"

#include "SireSystem/errors.h"

#include "constraint.h"

#include "delta.h"

#include "system.h"

#include <QDebug>

#include "constraint.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Constraint_class(){

    { //::SireSystem::Constraint
        typedef bp::class_< SireSystem::Constraint, bp::bases< SireBase::Property >, boost::noncopyable > Constraint_exposer_t;
        Constraint_exposer_t Constraint_exposer = Constraint_exposer_t( "Constraint", "This is the base class of all constraints. A constraint is an object\nthat is added to a System that tries to ensure that a condition is\nenforced. For example, a constraint could be used to change\nthe geometry of molecules with respect to lambda, or to change\nforcefield parameters with respect to alpha\n\nA Constraint class does its best to enforce a constraint - if it\ncant, then an exception is raised when a violation of the\nconstraint is detected (SireSystem::constraint_error)\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Constraint_scope( Constraint_exposer );
        { //::SireSystem::Constraint::apply
        
            typedef ::SireSystem::System ( ::SireSystem::Constraint::*apply_function_type)( ::SireSystem::System const & ) ;
            apply_function_type apply_function_value( &::SireSystem::Constraint::apply );
            
            Constraint_exposer.def( 
                "apply"
                , apply_function_value
                , ( bp::arg("system") )
                , "Apply this constraint to the passed system, returning\na new system in which the constraint is satisfied" );
        
        }
        { //::SireSystem::Constraint::assertSatisfied
        
            typedef void ( ::SireSystem::Constraint::*assertSatisfied_function_type)( ::SireSystem::System const & ) const;
            assertSatisfied_function_type assertSatisfied_function_value( &::SireSystem::Constraint::assertSatisfied );
            
            Constraint_exposer.def( 
                "assertSatisfied"
                , assertSatisfied_function_value
                , ( bp::arg("system") )
                , "Assert that the constraint is satisfied in the passed system\nThrow: SireSystem::constraint_error\n" );
        
        }
        { //::SireSystem::Constraint::isSatisfied
        
            typedef bool ( ::SireSystem::Constraint::*isSatisfied_function_type)( ::SireSystem::System const & ) const;
            isSatisfied_function_type isSatisfied_function_value( &::SireSystem::Constraint::isSatisfied );
            
            Constraint_exposer.def( 
                "isSatisfied"
                , isSatisfied_function_value
                , ( bp::arg("system") )
                , "Return whether or not this constraint is satisfied for\nthe passed system" );
        
        }
        { //::SireSystem::Constraint::mayAffect
        
            typedef bool ( ::SireSystem::Constraint::*mayAffect_function_type)( ::SireSystem::Delta const & ) const;
            mayAffect_function_type mayAffect_function_value( &::SireSystem::Constraint::mayAffect );
            
            Constraint_exposer.def( 
                "mayAffect"
                , mayAffect_function_value
                , ( bp::arg("delta") )
                , "Return whether or not this constraint may affect the passed delta" );
        
        }
        { //::SireSystem::Constraint::null
        
            typedef ::SireSystem::NullConstraint const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireSystem::Constraint::null );
            
            Constraint_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireSystem::Constraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::Constraint::typeName );
            
            Constraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        Constraint_exposer.staticmethod( "null" );
        Constraint_exposer.staticmethod( "typeName" );
        Constraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::Constraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::Constraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constraint_exposer.def_pickle(sire_pickle_suite< ::SireSystem::Constraint >());
        Constraint_exposer.def( "__str__", &__str__< ::SireSystem::Constraint > );
        Constraint_exposer.def( "__repr__", &__str__< ::SireSystem::Constraint > );
    }

}
