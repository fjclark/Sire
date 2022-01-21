// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MTSMC.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/symbol.h"

#include "SireError/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "mtsmc.h"

#include "mtsmc.h"

SireMove::MTSMC __copy__(const SireMove::MTSMC &other){ return SireMove::MTSMC(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_MTSMC_class(){

    { //::SireMove::MTSMC
        typedef bp::class_< SireMove::MTSMC, bp::bases< SireMove::MonteCarlo, SireMove::Move, SireBase::Property > > MTSMC_exposer_t;
        MTSMC_exposer_t MTSMC_exposer = MTSMC_exposer_t( "MTSMC", "This is a multiple-time-step Monte Carlo moves. This uses\nthe Metropolis-Hamilton acceptance test to perform\nnfast Monte Carlo moves on a fast Hamiltonian,\nand then testing whether the resulting configuration\nis suitable for inclusion in the ensemble generated\nusing a slow Hamiltonian.\n\nFor an application of this method to enhance sampling\nin QMMM simulations, see Woods, Manby and Mulholland,\nJ. Chem. Phys. 2008.\n\nAuthor: Christopher Woods\n", bp::init< bp::optional< SireBase::PropertyMap const & > >(( bp::arg("map")=SireBase::PropertyMap() ), "Null constructor") );
        bp::scope MTSMC_scope( MTSMC_exposer );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\n1 move using fastmoves for every slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, int, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("nfastmoves"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\nnfastmoves moves using fastmoves for every slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireCAS::Symbol const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("fastcomponent"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\n1 move using fastmoves operating on\nthe energy component fastcomponent for every slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireCAS::Symbol const &, int, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("fastcomponent"), bp::arg("nfastmoves"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\nnfastmoves moves using fastmoves operating on\nthe energy component fastcomponent for every slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireSystem::Constraints const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("slow_constraints"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\n1 move using fastmoves, then applies the constraints\nin slow_constraints for each slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireSystem::Constraints const &, int, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("slow_constraints"), bp::arg("nfastmoves"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\nnfastmoves moves using fastmoves, then applies the constraints\nin slow_constraints for each slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireCAS::Symbol const &, SireSystem::Constraints const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("fastcomponent"), bp::arg("slow_constraints"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\n1 move using fastmoves, using the energy component\nfastcomponent, then applies the constraints\nin slow_constraints for each slow move") );
        MTSMC_exposer.def( bp::init< SireMove::Moves const &, SireCAS::Symbol const &, SireSystem::Constraints const &, int, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("fastmoves"), bp::arg("fastcomponent"), bp::arg("slow_constraints"), bp::arg("nfastmoves"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a multiple time step Monte Carlo move that performs\nnfastmoves moves using fastmoves, using the energy component\nfastcomponent, then applies the constraints\nin slow_constraints for each slow move") );
        MTSMC_exposer.def( bp::init< SireMove::MTSMC const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::MTSMC::addSlowConstraint
        
            typedef void ( ::SireMove::MTSMC::*addSlowConstraint_function_type)( ::SireSystem::Constraint const & ) ;
            addSlowConstraint_function_type addSlowConstraint_function_value( &::SireMove::MTSMC::addSlowConstraint );
            
            MTSMC_exposer.def( 
                "addSlowConstraint"
                , addSlowConstraint_function_value
                , ( bp::arg("constraint") )
                , "Add a constraint that is applied at the end of each block\nof fast moves (i.e. before the slow move is tested)" );
        
        }
        { //::SireMove::MTSMC::clearStatistics
        
            typedef void ( ::SireMove::MTSMC::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::MTSMC::clearStatistics );
            
            MTSMC_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , "Completely clear all of the move statistics" );
        
        }
        { //::SireMove::MTSMC::fastEnergyComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::MTSMC::*fastEnergyComponent_function_type)(  ) const;
            fastEnergyComponent_function_type fastEnergyComponent_function_value( &::SireMove::MTSMC::fastEnergyComponent );
            
            MTSMC_exposer.def( 
                "fastEnergyComponent"
                , fastEnergyComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the energy component on which the fast moves will operate" );
        
        }
        { //::SireMove::MTSMC::fastMoves
        
            typedef ::SireMove::Moves const & ( ::SireMove::MTSMC::*fastMoves_function_type)(  ) const;
            fastMoves_function_type fastMoves_function_value( &::SireMove::MTSMC::fastMoves );
            
            MTSMC_exposer.def( 
                "fastMoves"
                , fastMoves_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the fast moves that will be performed for every\nslow move" );
        
        }
        { //::SireMove::MTSMC::move
        
            typedef void ( ::SireMove::MTSMC::*move_function_type)( ::SireSystem::System &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::MTSMC::move );
            
            MTSMC_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Perform the move - this will perform nfastmoves using fastmoves,\nand will then accept or reject the result based on the difference\nin the difference in energy between the fast and slow energies\nbefore and after the moves" );
        
        }
        { //::SireMove::MTSMC::nFastMoves
        
            typedef int ( ::SireMove::MTSMC::*nFastMoves_function_type)(  ) const;
            nFastMoves_function_type nFastMoves_function_value( &::SireMove::MTSMC::nFastMoves );
            
            MTSMC_exposer.def( 
                "nFastMoves"
                , nFastMoves_function_value
                , "Return the number of fast moves to perform per slow move" );
        
        }
        MTSMC_exposer.def( bp::self != bp::self );
        { //::SireMove::MTSMC::operator=
        
            typedef ::SireMove::MTSMC & ( ::SireMove::MTSMC::*assign_function_type)( ::SireMove::MTSMC const & ) ;
            assign_function_type assign_function_value( &::SireMove::MTSMC::operator= );
            
            MTSMC_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MTSMC_exposer.def( bp::self == bp::self );
        { //::SireMove::MTSMC::removeSlowConstraints
        
            typedef void ( ::SireMove::MTSMC::*removeSlowConstraints_function_type)(  ) ;
            removeSlowConstraints_function_type removeSlowConstraints_function_value( &::SireMove::MTSMC::removeSlowConstraints );
            
            MTSMC_exposer.def( 
                "removeSlowConstraints"
                , removeSlowConstraints_function_value
                , "Remove all of the slow constraints" );
        
        }
        { //::SireMove::MTSMC::setFastEnergyComponent
        
            typedef void ( ::SireMove::MTSMC::*setFastEnergyComponent_function_type)( ::SireCAS::Symbol const & ) ;
            setFastEnergyComponent_function_type setFastEnergyComponent_function_value( &::SireMove::MTSMC::setFastEnergyComponent );
            
            MTSMC_exposer.def( 
                "setFastEnergyComponent"
                , setFastEnergyComponent_function_value
                , ( bp::arg("component") )
                , "Set the energy component to be used for the fast moves" );
        
        }
        { //::SireMove::MTSMC::setFastMoves
        
            typedef void ( ::SireMove::MTSMC::*setFastMoves_function_type)( ::SireMove::Moves const & ) ;
            setFastMoves_function_type setFastMoves_function_value( &::SireMove::MTSMC::setFastMoves );
            
            MTSMC_exposer.def( 
                "setFastMoves"
                , setFastMoves_function_value
                , ( bp::arg("fastmoves") )
                , "Set the moves to be performed using the fast energy component.\nNote that these moves will be performd using the current\nfast energy component, which will override any energy\ncomponent currently set for these moves" );
        
        }
        { //::SireMove::MTSMC::setGenerator
        
            typedef void ( ::SireMove::MTSMC::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::MTSMC::setGenerator );
            
            MTSMC_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , "Set the random number generator used by this and all of the\ncontained moves" );
        
        }
        { //::SireMove::MTSMC::setNFastMoves
        
            typedef void ( ::SireMove::MTSMC::*setNFastMoves_function_type)( int ) ;
            setNFastMoves_function_type setNFastMoves_function_value( &::SireMove::MTSMC::setNFastMoves );
            
            MTSMC_exposer.def( 
                "setNFastMoves"
                , setNFastMoves_function_value
                , ( bp::arg("nfast") )
                , "Return the number of fast moves to perform per slow move" );
        
        }
        { //::SireMove::MTSMC::setSlowConstraints
        
            typedef void ( ::SireMove::MTSMC::*setSlowConstraints_function_type)( ::SireSystem::Constraints const & ) ;
            setSlowConstraints_function_type setSlowConstraints_function_value( &::SireMove::MTSMC::setSlowConstraints );
            
            MTSMC_exposer.def( 
                "setSlowConstraints"
                , setSlowConstraints_function_value
                , ( bp::arg("constraints") )
                , "Set the constraints that are applied at the end of each block\nof fast moves (i.e. before the slow move is tested).\nThis replaces any existing slow constraints" );
        
        }
        { //::SireMove::MTSMC::setSlowEnergyComponent
        
            typedef void ( ::SireMove::MTSMC::*setSlowEnergyComponent_function_type)( ::SireCAS::Symbol const & ) ;
            setSlowEnergyComponent_function_type setSlowEnergyComponent_function_value( &::SireMove::MTSMC::setSlowEnergyComponent );
            
            MTSMC_exposer.def( 
                "setSlowEnergyComponent"
                , setSlowEnergyComponent_function_value
                , ( bp::arg("component") )
                , "Set the energy component that will be used for the slow moves" );
        
        }
        { //::SireMove::MTSMC::slowConstraints
        
            typedef ::SireSystem::Constraints const & ( ::SireMove::MTSMC::*slowConstraints_function_type)(  ) const;
            slowConstraints_function_type slowConstraints_function_value( &::SireMove::MTSMC::slowConstraints );
            
            MTSMC_exposer.def( 
                "slowConstraints"
                , slowConstraints_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the constraints that are applied at the end of each\nblock of fast moves" );
        
        }
        { //::SireMove::MTSMC::slowEnergyComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::MTSMC::*slowEnergyComponent_function_type)(  ) const;
            slowEnergyComponent_function_type slowEnergyComponent_function_value( &::SireMove::MTSMC::slowEnergyComponent );
            
            MTSMC_exposer.def( 
                "slowEnergyComponent"
                , slowEnergyComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the energy component that will ultimately be used\nto generate the ensemble" );
        
        }
        { //::SireMove::MTSMC::toString
        
            typedef ::QString ( ::SireMove::MTSMC::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::MTSMC::toString );
            
            MTSMC_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this move" );
        
        }
        { //::SireMove::MTSMC::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::MTSMC::typeName );
            
            MTSMC_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        MTSMC_exposer.staticmethod( "typeName" );
        MTSMC_exposer.def( "__copy__", &__copy__);
        MTSMC_exposer.def( "__deepcopy__", &__copy__);
        MTSMC_exposer.def( "clone", &__copy__);
        MTSMC_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::MTSMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MTSMC_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::MTSMC >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MTSMC_exposer.def( "__getstate_manages_dict__", true);
        MTSMC_exposer.def( "__safe_for_unpickling__", true);
        MTSMC_exposer.def( "__setstate__", &__setstate__base64< ::SireMove::MTSMC > );
        MTSMC_exposer.def( "__getstate__", &__getstate__base64< ::SireMove::MTSMC > );
        MTSMC_exposer.def( "__str__", &__str__< ::SireMove::MTSMC > );
        MTSMC_exposer.def( "__repr__", &__str__< ::SireMove::MTSMC > );
    }

}
