// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AngleComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/ff.h"

#include "SireStream/datastream.h"

#include "internalcomponent.h"

#include "internalcomponent.h"

SireMM::AngleComponent __copy__(const SireMM::AngleComponent &other){ return SireMM::AngleComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AngleComponent_class(){

    { //::SireMM::AngleComponent
        typedef bp::class_< SireMM::AngleComponent, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > AngleComponent_exposer_t;
        AngleComponent_exposer_t AngleComponent_exposer = AngleComponent_exposer_t( "AngleComponent", "This class represents a Angle component of a forcefield", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("ffname")=SireFF::FFName() ), "Constructor") );
        bp::scope AngleComponent_scope( AngleComponent_exposer );
        AngleComponent_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("symbol") ), "Construct from a symbol\nThrow: SireError::incompatible_error\n") );
        AngleComponent_exposer.def( bp::init< SireMM::AngleComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::AngleComponent::changeEnergy
        
            typedef void ( ::SireMM::AngleComponent::*changeEnergy_function_type)( ::SireFF::FF &,::SireMM::AngleEnergy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireMM::AngleComponent::changeEnergy );
            
            AngleComponent_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("angnrg") )
                , "Change the component of the energy in the forcefield ff\nby delta" );
        
        }
        { //::SireMM::AngleComponent::setEnergy
        
            typedef void ( ::SireMM::AngleComponent::*setEnergy_function_type)( ::SireFF::FF &,::SireMM::AngleEnergy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireMM::AngleComponent::setEnergy );
            
            AngleComponent_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("angnrg") )
                , "Set the component of the energy in the forcefield ff\nto be equal to the passed energy" );
        
        }
        { //::SireMM::AngleComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::AngleComponent::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireMM::AngleComponent::symbols );
            
            AngleComponent_exposer.def( 
                "symbols"
                , symbols_function_value
                , "" );
        
        }
        { //::SireMM::AngleComponent::total
        
            typedef ::SireMM::AngleComponent const & ( ::SireMM::AngleComponent::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireMM::AngleComponent::total );
            
            AngleComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireMM::AngleComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::AngleComponent::typeName );
            
            AngleComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::AngleComponent::what
        
            typedef char const * ( ::SireMM::AngleComponent::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMM::AngleComponent::what );
            
            AngleComponent_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        AngleComponent_exposer.staticmethod( "typeName" );
        AngleComponent_exposer.def( "__copy__", &__copy__);
        AngleComponent_exposer.def( "__deepcopy__", &__copy__);
        AngleComponent_exposer.def( "clone", &__copy__);
        AngleComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AngleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AngleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleComponent_exposer.def( "__getstate_manages_dict__", true);
        AngleComponent_exposer.def( "__safe_for_unpickling__", true);
        AngleComponent_exposer.def( "__setstate__", &__setstate__base64< ::SireMM::AngleComponent > );
        AngleComponent_exposer.def( "__getstate__", &__getstate__base64< ::SireMM::AngleComponent > );
        AngleComponent_exposer.def( "__str__", &__str__< ::SireMM::AngleComponent > );
        AngleComponent_exposer.def( "__repr__", &__str__< ::SireMM::AngleComponent > );
        AngleComponent_exposer.def( "__hash__", &::SireMM::AngleComponent::hash );
    }

}
