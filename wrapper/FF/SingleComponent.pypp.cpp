// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "SingleComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "ff.h"

#include "ffcomponent.h"

#include <QRegExp>

#include "ffcomponent.h"

SireFF::SingleComponent __copy__(const SireFF::SingleComponent &other){ return SireFF::SingleComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SingleComponent_class(){

    { //::SireFF::SingleComponent
        typedef bp::class_< SireFF::SingleComponent, bp::bases< SireFF::FFComponent, SireCAS::Symbol, SireCAS::ExBase > > SingleComponent_exposer_t;
        SingleComponent_exposer_t SingleComponent_exposer = SingleComponent_exposer_t( "SingleComponent", "This class represents the single component of a single component forcefield.\nThis is provides a simple default class for simple, single-component\nforcefields\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< SireFF::FFName const & > >(( bp::arg("ffname")=SireFF::FFName() ), "Constructor") );
        bp::scope SingleComponent_scope( SingleComponent_exposer );
        SingleComponent_exposer.def( bp::init< SireFF::FFName const &, QString const & >(( bp::arg("ffname"), bp::arg("suffix") ), "Construct using the passed forcefield name and suffix") );
        SingleComponent_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("symbol") ), "Construct from a symbol\nThrow: SireError::incompatible_error\n") );
        SingleComponent_exposer.def( bp::init< SireFF::SingleComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::SingleComponent::changeEnergy
        
            typedef void ( ::SireFF::SingleComponent::*changeEnergy_function_type)( ::SireFF::FF &,::SireFF::SingleEnergy const & ) const;
            changeEnergy_function_type changeEnergy_function_value( &::SireFF::SingleComponent::changeEnergy );
            
            SingleComponent_exposer.def( 
                "changeEnergy"
                , changeEnergy_function_value
                , ( bp::arg("ff"), bp::arg("nrg") )
                , "Change the energy in the forcefield ff by delta" );
        
        }
        { //::SireFF::SingleComponent::setEnergy
        
            typedef void ( ::SireFF::SingleComponent::*setEnergy_function_type)( ::SireFF::FF &,::SireFF::SingleEnergy const & ) const;
            setEnergy_function_type setEnergy_function_value( &::SireFF::SingleComponent::setEnergy );
            
            SingleComponent_exposer.def( 
                "setEnergy"
                , setEnergy_function_value
                , ( bp::arg("ff"), bp::arg("nrg") )
                , "Set the energy in the forcefield ff to equal to the passed SingleEnergy" );
        
        }
        { //::SireFF::SingleComponent::symbols
        
            typedef ::SireCAS::Symbols ( ::SireFF::SingleComponent::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireFF::SingleComponent::symbols );
            
            SingleComponent_exposer.def( 
                "symbols"
                , symbols_function_value
                , "" );
        
        }
        { //::SireFF::SingleComponent::total
        
            typedef ::SireFF::SingleComponent const & ( ::SireFF::SingleComponent::*total_function_type)(  ) const;
            total_function_type total_function_value( &::SireFF::SingleComponent::total );
            
            SingleComponent_exposer.def( 
                "total"
                , total_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::SingleComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::SingleComponent::typeName );
            
            SingleComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::SingleComponent::what
        
            typedef char const * ( ::SireFF::SingleComponent::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::SingleComponent::what );
            
            SingleComponent_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        SingleComponent_exposer.staticmethod( "typeName" );
        SingleComponent_exposer.def( "__copy__", &__copy__);
        SingleComponent_exposer.def( "__deepcopy__", &__copy__);
        SingleComponent_exposer.def( "clone", &__copy__);
        SingleComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::SingleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SingleComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::SingleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SingleComponent_exposer.def_pickle(sire_pickle_suite< ::SireFF::SingleComponent >());
        SingleComponent_exposer.def( "__str__", &__str__< ::SireFF::SingleComponent > );
        SingleComponent_exposer.def( "__repr__", &__str__< ::SireFF::SingleComponent > );
        SingleComponent_exposer.def( "__hash__", &::SireFF::SingleComponent::hash );
    }

}
