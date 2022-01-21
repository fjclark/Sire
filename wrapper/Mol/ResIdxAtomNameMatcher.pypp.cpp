// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResIdxAtomNameMatcher.pypp.hpp"

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

SireMol::ResIdxAtomNameMatcher __copy__(const SireMol::ResIdxAtomNameMatcher &other){ return SireMol::ResIdxAtomNameMatcher(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResIdxAtomNameMatcher_class(){

    { //::SireMol::ResIdxAtomNameMatcher
        typedef bp::class_< SireMol::ResIdxAtomNameMatcher, bp::bases< SireMol::AtomMatcher, SireBase::Property > > ResIdxAtomNameMatcher_exposer_t;
        ResIdxAtomNameMatcher_exposer_t ResIdxAtomNameMatcher_exposer = ResIdxAtomNameMatcher_exposer_t( "ResIdxAtomNameMatcher", "Match atoms by name within each residue.\n\nAuthor: Lester Hedges\n", bp::init< >("Constructor") );
        bp::scope ResIdxAtomNameMatcher_scope( ResIdxAtomNameMatcher_exposer );
        ResIdxAtomNameMatcher_exposer.def( bp::init< SireMol::ResIdxAtomNameMatcher const & >(( bp::arg("arg0") ), "Copy constructor") );
        ResIdxAtomNameMatcher_exposer.def( bp::self != bp::self );
        { //::SireMol::ResIdxAtomNameMatcher::operator=
        
            typedef ::SireMol::ResIdxAtomNameMatcher & ( ::SireMol::ResIdxAtomNameMatcher::*assign_function_type)( ::SireMol::ResIdxAtomNameMatcher const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResIdxAtomNameMatcher::operator= );
            
            ResIdxAtomNameMatcher_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResIdxAtomNameMatcher_exposer.def( bp::self == bp::self );
        { //::SireMol::ResIdxAtomNameMatcher::toString
        
            typedef ::QString ( ::SireMol::ResIdxAtomNameMatcher::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResIdxAtomNameMatcher::toString );
            
            ResIdxAtomNameMatcher_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::ResIdxAtomNameMatcher::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResIdxAtomNameMatcher::typeName );
            
            ResIdxAtomNameMatcher_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ResIdxAtomNameMatcher::what
        
            typedef char const * ( ::SireMol::ResIdxAtomNameMatcher::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ResIdxAtomNameMatcher::what );
            
            ResIdxAtomNameMatcher_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ResIdxAtomNameMatcher_exposer.staticmethod( "typeName" );
        ResIdxAtomNameMatcher_exposer.def( "__copy__", &__copy__);
        ResIdxAtomNameMatcher_exposer.def( "__deepcopy__", &__copy__);
        ResIdxAtomNameMatcher_exposer.def( "clone", &__copy__);
        ResIdxAtomNameMatcher_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResIdxAtomNameMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdxAtomNameMatcher_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResIdxAtomNameMatcher >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResIdxAtomNameMatcher_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::ResIdxAtomNameMatcher > );
        ResIdxAtomNameMatcher_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::ResIdxAtomNameMatcher > );
        ResIdxAtomNameMatcher_exposer.def( "__str__", &__str__< ::SireMol::ResIdxAtomNameMatcher > );
        ResIdxAtomNameMatcher_exposer.def( "__repr__", &__str__< ::SireMol::ResIdxAtomNameMatcher > );
    }

}
