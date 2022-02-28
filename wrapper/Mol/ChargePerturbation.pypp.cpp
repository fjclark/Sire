// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChargePerturbation.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/values.h"

#include "SireStream/datastream.h"

#include "atomcharges.h"

#include "chargeperturbation.h"

#include "molecule.h"

#include "moleditor.h"

#include "mover.hpp"

#include "chargeperturbation.h"

SireMol::ChargePerturbation __copy__(const SireMol::ChargePerturbation &other){ return SireMol::ChargePerturbation(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChargePerturbation_class(){

    { //::SireMol::ChargePerturbation
        typedef bp::class_< SireMol::ChargePerturbation, bp::bases< SireMol::Perturbation, SireBase::Property > > ChargePerturbation_exposer_t;
        ChargePerturbation_exposer_t ChargePerturbation_exposer = ChargePerturbation_exposer_t( "ChargePerturbation", "This perturbation is used to scale charges from one value\nto another as a function of lambda\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor - this creates a charge perturbation that\nperturbs from charges in initial_charge to charges in\nfinal_charge, placing the current charges in charge,\nand using Perturbation::defaultEquation() to map the\ncharges") );
        bp::scope ChargePerturbation_scope( ChargePerturbation_exposer );
        ChargePerturbation_exposer.def( bp::init< SireBase::PropertyMap const & >(( bp::arg("map") ), "Construct, using the passed map to find the properties used\nby this perturbation") );
        ChargePerturbation_exposer.def( bp::init< SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("mapping_function"), bp::arg("map")=SireBase::PropertyMap() ), "Construct, using the passed map to find the properties used\nby this perturbation and the passed mapping function to map\nthe charges between the states") );
        ChargePerturbation_exposer.def( bp::init< SireMol::ChargePerturbation const & >(( bp::arg("other") ), "Copy constructor") );
        ChargePerturbation_exposer.def( bp::self != bp::self );
        { //::SireMol::ChargePerturbation::operator=
        
            typedef ::SireMol::ChargePerturbation & ( ::SireMol::ChargePerturbation::*assign_function_type)( ::SireMol::ChargePerturbation const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChargePerturbation::operator= );
            
            ChargePerturbation_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChargePerturbation_exposer.def( bp::self == bp::self );
        { //::SireMol::ChargePerturbation::requiredProperties
        
            typedef ::QSet< QString > ( ::SireMol::ChargePerturbation::*requiredProperties_function_type)(  ) const;
            requiredProperties_function_type requiredProperties_function_value( &::SireMol::ChargePerturbation::requiredProperties );
            
            ChargePerturbation_exposer.def( 
                "requiredProperties"
                , requiredProperties_function_value
                , "Return the properties required or changed by this perturbation" );
        
        }
        { //::SireMol::ChargePerturbation::toString
        
            typedef ::QString ( ::SireMol::ChargePerturbation::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChargePerturbation::toString );
            
            ChargePerturbation_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::ChargePerturbation::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChargePerturbation::typeName );
            
            ChargePerturbation_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ChargePerturbation::wouldChange
        
            typedef bool ( ::SireMol::ChargePerturbation::*wouldChange_function_type)( ::SireMol::Molecule const &,::SireCAS::Values const & ) const;
            wouldChange_function_type wouldChange_function_value( &::SireMol::ChargePerturbation::wouldChange );
            
            ChargePerturbation_exposer.def( 
                "wouldChange"
                , wouldChange_function_value
                , ( bp::arg("molecule"), bp::arg("values") )
                , "Return whether or not this perturbation with the passed values would\nchange the molecule molecule" );
        
        }
        ChargePerturbation_exposer.staticmethod( "typeName" );
        ChargePerturbation_exposer.def( "__copy__", &__copy__);
        ChargePerturbation_exposer.def( "__deepcopy__", &__copy__);
        ChargePerturbation_exposer.def( "clone", &__copy__);
        ChargePerturbation_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChargePerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChargePerturbation_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChargePerturbation >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChargePerturbation_exposer.def_pickle(sire_pickle_suite< ::SireMol::ChargePerturbation >());
        ChargePerturbation_exposer.def( "__str__", &__str__< ::SireMol::ChargePerturbation > );
        ChargePerturbation_exposer.def( "__repr__", &__str__< ::SireMol::ChargePerturbation > );
    }

}
