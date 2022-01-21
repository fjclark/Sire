// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ResIdxAtomMCSMatcher.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireUnits/units.h"

#include "atom.h"

#include "atomidentifier.h"

#include "atomidx.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "atomname.h"

#include "atomselection.h"

#include "evaluator.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.h"

#include "selector.hpp"

#include "tostring.h"

#include "atommatchers.h"

SireMol::ResIdxAtomMCSMatcher __copy__(const SireMol::ResIdxAtomMCSMatcher &other){ return SireMol::ResIdxAtomMCSMatcher(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResIdxAtomMCSMatcher_class(){

    { //::SireMol::ResIdxAtomMCSMatcher
        typedef bp::class_< SireMol::ResIdxAtomMCSMatcher, bp::bases< SireMol::AtomMatcher, SireBase::Property > > ResIdxAtomMCSMatcher_exposer_t;
        ResIdxAtomMCSMatcher_exposer_t ResIdxAtomMCSMatcher_exposer = ResIdxAtomMCSMatcher_exposer_t( "ResIdxAtomMCSMatcher", "Match atoms by name MCS within each residue.\n\nAuthor: Lester Hedges\n", bp::init< >("Constructor") );
        bp::scope ResIdxAtomMCSMatcher_scope( ResIdxAtomMCSMatcher_exposer );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< bool >(( bp::arg("verbose") ), "Constructor") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireUnits::Dimension::Time const &, bool >(( bp::arg("timeout"), bp::arg("verbose") ), "Construct specifying the timeout for the MCS match") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireMol::AtomMatcher const &, bool >(( bp::arg("prematcher"), bp::arg("verbose") ), "Construct specifying the prematcher for the MCS match") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireMol::AtomMatcher const &, SireUnits::Dimension::Time const &, bool >(( bp::arg("prematcher"), bp::arg("timeout"), bp::arg("verbose") ), "Construct specifying the timeout and prematcher for the MCS match") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< bool, bool >(( bp::arg("match_light_atoms"), bp::arg("verbose") ), "Constructor, specifying whether or not to match light atoms") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireUnits::Dimension::Time const &, bool, bool >(( bp::arg("timeout"), bp::arg("match_light_atoms"), bp::arg("verbose") ), "Construct specifying the timeout for the MCS match, and specifying whether or not\nto match light atoms") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireMol::AtomMatcher const &, bool, bool >(( bp::arg("prematcher"), bp::arg("match_light_atoms"), bp::arg("verbose") ), "Construct specifying the prematcher for the MCS match,\nand specifying whether or not to match light atoms\n") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireMol::AtomMatcher const &, SireUnits::Dimension::Time const &, bool, bool >(( bp::arg("prematcher"), bp::arg("timeout"), bp::arg("match_light_atoms"), bp::arg("verbose") ), "Construct specifying the timeout and prematcher for the MCS match,\nand specifying whether or not to match light atoms\n") );
        ResIdxAtomMCSMatcher_exposer.def( bp::init< SireMol::ResIdxAtomMCSMatcher const & >(( bp::arg("arg0") ), "Copy constructor") );
        { //::SireMol::ResIdxAtomMCSMatcher::isVerbose
        
            typedef bool ( ::SireMol::ResIdxAtomMCSMatcher::*isVerbose_function_type)(  ) const;
            isVerbose_function_type isVerbose_function_value( &::SireMol::ResIdxAtomMCSMatcher::isVerbose );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "isVerbose"
                , isVerbose_function_value
                , "Return whether or not this report progress to stdout." );
        
        }
        { //::SireMol::ResIdxAtomMCSMatcher::matchingLightAtoms
        
            typedef bool ( ::SireMol::ResIdxAtomMCSMatcher::*matchingLightAtoms_function_type)(  ) const;
            matchingLightAtoms_function_type matchingLightAtoms_function_value( &::SireMol::ResIdxAtomMCSMatcher::matchingLightAtoms );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "matchingLightAtoms"
                , matchingLightAtoms_function_value
                , "Return whether or not this will include light atoms (e.g. hydrogen)\nwhen searching for the maximum common substructure" );
        
        }
        ResIdxAtomMCSMatcher_exposer.def( bp::self != bp::self );
        { //::SireMol::ResIdxAtomMCSMatcher::operator=
        
            typedef ::SireMol::ResIdxAtomMCSMatcher & ( ::SireMol::ResIdxAtomMCSMatcher::*assign_function_type)( ::SireMol::ResIdxAtomMCSMatcher const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResIdxAtomMCSMatcher::operator= );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResIdxAtomMCSMatcher_exposer.def( bp::self == bp::self );
        { //::SireMol::ResIdxAtomMCSMatcher::preMatcher
        
            typedef ::SireMol::AtomMatcher const & ( ::SireMol::ResIdxAtomMCSMatcher::*preMatcher_function_type)(  ) const;
            preMatcher_function_type preMatcher_function_value( &::SireMol::ResIdxAtomMCSMatcher::preMatcher );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "preMatcher"
                , preMatcher_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the prematcher (if any) that is used to pre-match atoms\nbefore the MCS match" );
        
        }
        { //::SireMol::ResIdxAtomMCSMatcher::timeout
        
            typedef ::SireUnits::Dimension::Time ( ::SireMol::ResIdxAtomMCSMatcher::*timeout_function_type)(  ) const;
            timeout_function_type timeout_function_value( &::SireMol::ResIdxAtomMCSMatcher::timeout );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "timeout"
                , timeout_function_value
                , "Return the timeout before the MCS match is abandoned" );
        
        }
        { //::SireMol::ResIdxAtomMCSMatcher::toString
        
            typedef ::QString ( ::SireMol::ResIdxAtomMCSMatcher::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResIdxAtomMCSMatcher::toString );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::ResIdxAtomMCSMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResIdxAtomMCSMatcher::typeName );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ResIdxAtomMCSMatcher::what
        
            typedef char const * ( ::SireMol::ResIdxAtomMCSMatcher::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ResIdxAtomMCSMatcher::what );
            
            ResIdxAtomMCSMatcher_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ResIdxAtomMCSMatcher_exposer.staticmethod( "typeName" );
        ResIdxAtomMCSMatcher_exposer.def( "__copy__", &__copy__);
        ResIdxAtomMCSMatcher_exposer.def( "__deepcopy__", &__copy__);
        ResIdxAtomMCSMatcher_exposer.def( "clone", &__copy__);
        ResIdxAtomMCSMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResIdxAtomMCSMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdxAtomMCSMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResIdxAtomMCSMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdxAtomMCSMatcher_exposer.def( "__getstate_manages_dict__", true);
        ResIdxAtomMCSMatcher_exposer.def( "__safe_for_unpickling__", true);
        ResIdxAtomMCSMatcher_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::ResIdxAtomMCSMatcher > );
        ResIdxAtomMCSMatcher_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::ResIdxAtomMCSMatcher > );
        ResIdxAtomMCSMatcher_exposer.def( "__str__", &__str__< ::SireMol::ResIdxAtomMCSMatcher > );
        ResIdxAtomMCSMatcher_exposer.def( "__repr__", &__str__< ::SireMol::ResIdxAtomMCSMatcher > );
    }

}
