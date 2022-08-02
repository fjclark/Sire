// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CLJCutoffFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireBase/lengthproperty.h"

#include "SireBase/numberproperty.h"

#include "SireBase/properties.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireMaths/multidouble.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/cartesian.h"

#include "SireVol/gridinfo.h"

#include "SireVol/periodicbox.h"

#include "SireVol/triclinicbox.h"

#include "cljboxes.h"

#include "cljfunction.h"

#include "switchingfunction.h"

#include "tbb/blocked_range.h"

#include "tbb/parallel_for.h"

#include "tostring.h"

#include <QElapsedTimer>

#include "cljfunction.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CLJCutoffFunction_class(){

    { //::SireMM::CLJCutoffFunction
        typedef bp::class_< SireMM::CLJCutoffFunction, bp::bases< SireMM::CLJFunction, SireBase::Property >, boost::noncopyable > CLJCutoffFunction_exposer_t;
        CLJCutoffFunction_exposer_t CLJCutoffFunction_exposer = CLJCutoffFunction_exposer_t( "CLJCutoffFunction", "This is the base class of all CLJ functions that have a cutoff\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope CLJCutoffFunction_scope( CLJCutoffFunction_exposer );
        { //::SireMM::CLJCutoffFunction::containsProperty
        
            typedef bool ( ::SireMM::CLJCutoffFunction::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CLJCutoffFunction::containsProperty );
            
            CLJCutoffFunction_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return whether or not this function contains a property called name" );
        
        }
        { //::SireMM::CLJCutoffFunction::coulombCutoff
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::CLJCutoffFunction::*coulombCutoff_function_type)(  ) const;
            coulombCutoff_function_type coulombCutoff_function_value( &::SireMM::CLJCutoffFunction::coulombCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "coulombCutoff"
                , coulombCutoff_function_value
                , bp::release_gil_policy()
                , "Return the coulomb cutoff distance" );
        
        }
        { //::SireMM::CLJCutoffFunction::hasCutoff
        
            typedef bool ( ::SireMM::CLJCutoffFunction::*hasCutoff_function_type)(  ) const;
            hasCutoff_function_type hasCutoff_function_value( &::SireMM::CLJCutoffFunction::hasCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "hasCutoff"
                , hasCutoff_function_value
                , bp::release_gil_policy()
                , "Return whether or not this function has a cutoff" );
        
        }
        { //::SireMM::CLJCutoffFunction::ljCutoff
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::CLJCutoffFunction::*ljCutoff_function_type)(  ) const;
            ljCutoff_function_type ljCutoff_function_value( &::SireMM::CLJCutoffFunction::ljCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "ljCutoff"
                , ljCutoff_function_value
                , bp::release_gil_policy()
                , "Return the LJ cutoff distance" );
        
        }
        { //::SireMM::CLJCutoffFunction::properties
        
            typedef ::SireBase::Properties ( ::SireMM::CLJCutoffFunction::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CLJCutoffFunction::properties );
            
            CLJCutoffFunction_exposer.def( 
                "properties"
                , properties_function_value
                , bp::release_gil_policy()
                , "Return the properties that can be set in this function" );
        
        }
        { //::SireMM::CLJCutoffFunction::property
        
            typedef ::SireBase::PropertyPtr ( ::SireMM::CLJCutoffFunction::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CLJCutoffFunction::property );
            
            CLJCutoffFunction_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Return the value of the property with name name" );
        
        }
        { //::SireMM::CLJCutoffFunction::setCoulombCutoff
        
            typedef void ( ::SireMM::CLJCutoffFunction::*setCoulombCutoff_function_type)( ::SireUnits::Dimension::Length ) ;
            setCoulombCutoff_function_type setCoulombCutoff_function_value( &::SireMM::CLJCutoffFunction::setCoulombCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "setCoulombCutoff"
                , setCoulombCutoff_function_value
                , ( bp::arg("distance") )
                , bp::release_gil_policy()
                , "Set the coulomb cutoff to the specified distance" );
        
        }
        { //::SireMM::CLJCutoffFunction::setCutoff
        
            typedef void ( ::SireMM::CLJCutoffFunction::*setCutoff_function_type)( ::SireUnits::Dimension::Length ) ;
            setCutoff_function_type setCutoff_function_value( &::SireMM::CLJCutoffFunction::setCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "setCutoff"
                , setCutoff_function_value
                , ( bp::arg("distance") )
                , bp::release_gil_policy()
                , "Set the coulomb and LJ cutoff distances to distance" );
        
        }
        { //::SireMM::CLJCutoffFunction::setCutoff
        
            typedef void ( ::SireMM::CLJCutoffFunction::*setCutoff_function_type)( ::SireUnits::Dimension::Length,::SireUnits::Dimension::Length ) ;
            setCutoff_function_type setCutoff_function_value( &::SireMM::CLJCutoffFunction::setCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "setCutoff"
                , setCutoff_function_value
                , ( bp::arg("coulomb_cutoff"), bp::arg("lj_cutoff") )
                , bp::release_gil_policy()
                , "Set the coulomb and LJ cutoff distances to the specified values" );
        
        }
        { //::SireMM::CLJCutoffFunction::setLJCutoff
        
            typedef void ( ::SireMM::CLJCutoffFunction::*setLJCutoff_function_type)( ::SireUnits::Dimension::Length ) ;
            setLJCutoff_function_type setLJCutoff_function_value( &::SireMM::CLJCutoffFunction::setLJCutoff );
            
            CLJCutoffFunction_exposer.def( 
                "setLJCutoff"
                , setLJCutoff_function_value
                , ( bp::arg("distance") )
                , bp::release_gil_policy()
                , "Set the LJ cutoff to the specified distance" );
        
        }
        { //::SireMM::CLJCutoffFunction::setProperty
        
            typedef ::SireMM::CLJFunctionPtr ( ::SireMM::CLJCutoffFunction::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) const;
            setProperty_function_type setProperty_function_value( &::SireMM::CLJCutoffFunction::setProperty );
            
            CLJCutoffFunction_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the property with name name to value value" );
        
        }
        { //::SireMM::CLJCutoffFunction::toString
        
            typedef ::QString ( ::SireMM::CLJCutoffFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJCutoffFunction::toString );
            
            CLJCutoffFunction_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::CLJCutoffFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJCutoffFunction::typeName );
            
            CLJCutoffFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CLJCutoffFunction_exposer.staticmethod( "typeName" );
        CLJCutoffFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJCutoffFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJCutoffFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJCutoffFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJCutoffFunction_exposer.def_pickle(sire_pickle_suite< ::SireMM::CLJCutoffFunction >());
        CLJCutoffFunction_exposer.def( "__str__", &__str__< ::SireMM::CLJCutoffFunction > );
        CLJCutoffFunction_exposer.def( "__repr__", &__str__< ::SireMM::CLJCutoffFunction > );
    }

}
