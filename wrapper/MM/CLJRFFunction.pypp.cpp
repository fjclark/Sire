// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CLJRFFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/numberproperty.h"

#include "SireError/errors.h"

#include "SireMaths/multidouble.h"

#include "SireMaths/multifloat.h"

#include "SireMaths/multiint.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "SireVol/gridinfo.h"

#include "cljrffunction.h"

#include <QDebug>

#include <QElapsedTimer>

#include "cljrffunction.h"

SireMM::CLJRFFunction __copy__(const SireMM::CLJRFFunction &other){ return SireMM::CLJRFFunction(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CLJRFFunction_class(){

    { //::SireMM::CLJRFFunction
        typedef bp::class_< SireMM::CLJRFFunction, bp::bases< SireMM::CLJCutoffFunction, SireMM::CLJFunction, SireBase::Property > > CLJRFFunction_exposer_t;
        CLJRFFunction_exposer_t CLJRFFunction_exposer = CLJRFFunction_exposer_t( "CLJRFFunction", "This CLJFunction calculates the intermolecular coulomb and LJ energy of the passed\nCLJAtoms using a reaction field cutoff\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope CLJRFFunction_scope( CLJRFFunction_exposer );
        CLJRFFunction_exposer.def( bp::init< SireUnits::Dimension::Length >(( bp::arg("cutoff") ), "Copy constructor") );
        CLJRFFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireUnits::Dimension::Length >(( bp::arg("coul_cutoff"), bp::arg("lj_cutoff") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length >(( bp::arg("space"), bp::arg("cutoff") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireUnits::Dimension::Length >(( bp::arg("space"), bp::arg("coul_cutoff"), bp::arg("lj_cutoff") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("cutoff"), bp::arg("combining_rules") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireUnits::Dimension::Length, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("coul_cutoff"), bp::arg("lj_cutoff"), bp::arg("combining_rules") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireVol::Space const &, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("combining_rules") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("cutoff"), bp::arg("combining_rules") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireVol::Space const &, SireUnits::Dimension::Length, SireUnits::Dimension::Length, SireMM::CLJFunction::COMBINING_RULES >(( bp::arg("space"), bp::arg("coul_cutoff"), bp::arg("lj_cutoff"), bp::arg("combining_rules") ), "") );
        CLJRFFunction_exposer.def( bp::init< SireMM::CLJRFFunction const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::CLJRFFunction::containsProperty
        
            typedef bool ( ::SireMM::CLJRFFunction::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireMM::CLJRFFunction::containsProperty );
            
            CLJRFFunction_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "Return whether or not this function contains a property called name" );
        
        }
        { //::SireMM::CLJRFFunction::defaultRFFunction
        
            typedef ::SireMM::CLJFunctionPtr ( *defaultRFFunction_function_type )(  );
            defaultRFFunction_function_type defaultRFFunction_function_value( &::SireMM::CLJRFFunction::defaultRFFunction );
            
            CLJRFFunction_exposer.def( 
                "defaultRFFunction"
                , defaultRFFunction_function_value
                , "" );
        
        }
        { //::SireMM::CLJRFFunction::dielectric
        
            typedef float ( ::SireMM::CLJRFFunction::*dielectric_function_type)(  ) const;
            dielectric_function_type dielectric_function_value( &::SireMM::CLJRFFunction::dielectric );
            
            CLJRFFunction_exposer.def( 
                "dielectric"
                , dielectric_function_value
                , "Return the value of the dielectric constant" );
        
        }
        CLJRFFunction_exposer.def( bp::self != bp::self );
        { //::SireMM::CLJRFFunction::operator=
        
            typedef ::SireMM::CLJRFFunction & ( ::SireMM::CLJRFFunction::*assign_function_type)( ::SireMM::CLJRFFunction const & ) ;
            assign_function_type assign_function_value( &::SireMM::CLJRFFunction::operator= );
            
            CLJRFFunction_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CLJRFFunction_exposer.def( bp::self == bp::self );
        { //::SireMM::CLJRFFunction::properties
        
            typedef ::SireBase::Properties ( ::SireMM::CLJRFFunction::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMM::CLJRFFunction::properties );
            
            CLJRFFunction_exposer.def( 
                "properties"
                , properties_function_value
                , "Return the properties of this function" );
        
        }
        { //::SireMM::CLJRFFunction::property
        
            typedef ::SireBase::PropertyPtr ( ::SireMM::CLJRFFunction::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireMM::CLJRFFunction::property );
            
            CLJRFFunction_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , "Return the value of the property with name name" );
        
        }
        { //::SireMM::CLJRFFunction::setDielectric
        
            typedef void ( ::SireMM::CLJRFFunction::*setDielectric_function_type)( float ) ;
            setDielectric_function_type setDielectric_function_value( &::SireMM::CLJRFFunction::setDielectric );
            
            CLJRFFunction_exposer.def( 
                "setDielectric"
                , setDielectric_function_value
                , ( bp::arg("dielectric") )
                , "Set the dielectric constant to dielectric" );
        
        }
        { //::SireMM::CLJRFFunction::setProperty
        
            typedef ::SireMM::CLJFunctionPtr ( ::SireMM::CLJRFFunction::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) const;
            setProperty_function_type setProperty_function_value( &::SireMM::CLJRFFunction::setProperty );
            
            CLJRFFunction_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("value") )
                , "Return a copy of this function where the property name has been set to the\nvalue value" );
        
        }
        { //::SireMM::CLJRFFunction::supportsGridCalculation
        
            typedef bool ( ::SireMM::CLJRFFunction::*supportsGridCalculation_function_type)(  ) const;
            supportsGridCalculation_function_type supportsGridCalculation_function_value( &::SireMM::CLJRFFunction::supportsGridCalculation );
            
            CLJRFFunction_exposer.def( 
                "supportsGridCalculation"
                , supportsGridCalculation_function_value
                , "This function does support calculations using a grid" );
        
        }
        { //::SireMM::CLJRFFunction::toString
        
            typedef ::QString ( ::SireMM::CLJRFFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::CLJRFFunction::toString );
            
            CLJRFFunction_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMM::CLJRFFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::CLJRFFunction::typeName );
            
            CLJRFFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::CLJRFFunction::what
        
            typedef char const * ( ::SireMM::CLJRFFunction::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::CLJRFFunction::what );
            
            CLJRFFunction_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CLJRFFunction_exposer.staticmethod( "defaultRFFunction" );
        CLJRFFunction_exposer.staticmethod( "typeName" );
        CLJRFFunction_exposer.def( "__copy__", &__copy__);
        CLJRFFunction_exposer.def( "__deepcopy__", &__copy__);
        CLJRFFunction_exposer.def( "clone", &__copy__);
        CLJRFFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::CLJRFFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJRFFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::CLJRFFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CLJRFFunction_exposer.def( "__getstate_manages_dict__", true);
        CLJRFFunction_exposer.def( "__safe_for_unpickling__", true);
        CLJRFFunction_exposer.def( "__setstate__", &__setstate__base64< ::SireMM::CLJRFFunction > );
        CLJRFFunction_exposer.def( "__getstate__", &__getstate__base64< ::SireMM::CLJRFFunction > );
        CLJRFFunction_exposer.def( "__str__", &__str__< ::SireMM::CLJRFFunction > );
        CLJRFFunction_exposer.def( "__repr__", &__str__< ::SireMM::CLJRFFunction > );
    }

}
