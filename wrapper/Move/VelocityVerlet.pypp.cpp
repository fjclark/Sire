// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "VelocityVerlet.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/forcetable.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atommasses.h"

#include "SireMol/core.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/convert.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "velocityverlet.h"

#include "velocityverlet.h"

SireMove::VelocityVerlet __copy__(const SireMove::VelocityVerlet &other){ return SireMove::VelocityVerlet(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_VelocityVerlet_class(){

    { //::SireMove::VelocityVerlet
        typedef bp::class_< SireMove::VelocityVerlet, bp::bases< SireMove::Integrator, SireBase::Property > > VelocityVerlet_exposer_t;
        VelocityVerlet_exposer_t VelocityVerlet_exposer = VelocityVerlet_exposer_t( "VelocityVerlet", "This class implements an atomistic velocity verlet dynamics integrator\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< bool > >(( bp::arg("frequent_save_velocities")=(bool)(false) ), "Constructor") );
        bp::scope VelocityVerlet_scope( VelocityVerlet_exposer );
        VelocityVerlet_exposer.def( bp::init< SireMove::VelocityVerlet const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::VelocityVerlet::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::VelocityVerlet::*createWorkspace_function_type)( ::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::VelocityVerlet::createWorkspace );
            
            VelocityVerlet_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() )
                , "Create an empty workspace" );
        
        }
        { //::SireMove::VelocityVerlet::createWorkspace
        
            typedef ::SireMove::IntegratorWorkspacePtr ( ::SireMove::VelocityVerlet::*createWorkspace_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            createWorkspace_function_type createWorkspace_function_value( &::SireMove::VelocityVerlet::createWorkspace );
            
            VelocityVerlet_exposer.def( 
                "createWorkspace"
                , createWorkspace_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , "Create a workspace for this integrator for the molecule group molgroup" );
        
        }
        { //::SireMove::VelocityVerlet::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::VelocityVerlet::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::VelocityVerlet::ensemble );
            
            VelocityVerlet_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "Return the ensemble of this integrator" );
        
        }
        { //::SireMove::VelocityVerlet::integrate
        
            typedef void ( ::SireMove::VelocityVerlet::*integrate_function_type)( ::SireMove::IntegratorWorkspace &,::SireCAS::Symbol const &,::SireUnits::Dimension::Time,int,bool ) ;
            integrate_function_type integrate_function_value( &::SireMove::VelocityVerlet::integrate );
            
            VelocityVerlet_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("workspace"), bp::arg("nrg_component"), bp::arg("timestep"), bp::arg("nmoves"), bp::arg("record_stats") )
                , bp::release_gil_policy()
                , "Integrate the coordinates of the atoms in the molecules in molgroup\nusing the forces in forcetable, using the optionally supplied\nproperty map to find the necessary molecular properties\nThrow: SireMol::missing_molecule\nThrow: SireBase::missing_property\nThrow: SireError:invalid_cast\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::VelocityVerlet::isTimeReversible
        
            typedef bool ( ::SireMove::VelocityVerlet::*isTimeReversible_function_type)(  ) const;
            isTimeReversible_function_type isTimeReversible_function_value( &::SireMove::VelocityVerlet::isTimeReversible );
            
            VelocityVerlet_exposer.def( 
                "isTimeReversible"
                , isTimeReversible_function_value
                , bp::release_gil_policy()
                , "Return whether or not this integrator is time-reversible" );
        
        }
        VelocityVerlet_exposer.def( bp::self != bp::self );
        { //::SireMove::VelocityVerlet::operator=
        
            typedef ::SireMove::VelocityVerlet & ( ::SireMove::VelocityVerlet::*assign_function_type)( ::SireMove::VelocityVerlet const & ) ;
            assign_function_type assign_function_value( &::SireMove::VelocityVerlet::operator= );
            
            VelocityVerlet_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        VelocityVerlet_exposer.def( bp::self == bp::self );
        { //::SireMove::VelocityVerlet::toString
        
            typedef ::QString ( ::SireMove::VelocityVerlet::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::VelocityVerlet::toString );
            
            VelocityVerlet_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this integrator" );
        
        }
        { //::SireMove::VelocityVerlet::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::VelocityVerlet::typeName );
            
            VelocityVerlet_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        VelocityVerlet_exposer.staticmethod( "typeName" );
        VelocityVerlet_exposer.def( "__copy__", &__copy__);
        VelocityVerlet_exposer.def( "__deepcopy__", &__copy__);
        VelocityVerlet_exposer.def( "clone", &__copy__);
        VelocityVerlet_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::VelocityVerlet >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VelocityVerlet_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::VelocityVerlet >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VelocityVerlet_exposer.def_pickle(sire_pickle_suite< ::SireMove::VelocityVerlet >());
        VelocityVerlet_exposer.def( "__str__", &__str__< ::SireMove::VelocityVerlet > );
        VelocityVerlet_exposer.def( "__repr__", &__str__< ::SireMove::VelocityVerlet > );
    }

}
