// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Moves.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/constraints.h"

#include "SireSystem/errors.h"

#include "SireSystem/system.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "ensemble.h"

#include "moves.h"

#include "simulation.h"

#include <QDebug>

#include <QElapsedTimer>

#include <QMutex>

#include "moves.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Moves_class(){

    { //::SireMove::Moves
        typedef bp::class_< SireMove::Moves, bp::bases< SireBase::Property >, boost::noncopyable > Moves_exposer_t;
        Moves_exposer_t Moves_exposer = Moves_exposer_t( "Moves", "This is the base class of all Moves objects. These are objects\nthat contain a collection of moves that are applied to a system\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Moves_scope( Moves_exposer );
        { //::SireMove::Moves::acceptableDelta
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::Moves::*acceptableDelta_function_type)(  ) const;
            acceptableDelta_function_type acceptableDelta_function_value( &::SireMove::Moves::acceptableDelta );
            
            Moves_exposer.def( 
                "acceptableDelta"
                , acceptableDelta_function_value
                , bp::release_gil_policy()
                , "Return the acceptable level of drift in the running total from the recalculated total energy" );
        
        }
        { //::SireMove::Moves::checkingRunningTotal
        
            typedef bool ( ::SireMove::Moves::*checkingRunningTotal_function_type)(  ) const;
            checkingRunningTotal_function_type checkingRunningTotal_function_value( &::SireMove::Moves::checkingRunningTotal );
            
            Moves_exposer.def( 
                "checkingRunningTotal"
                , checkingRunningTotal_function_value
                , bp::release_gil_policy()
                , "Return whether or not we are checking of the running total energy" );
        
        }
        { //::SireMove::Moves::chemicalPotential
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::Moves::*chemicalPotential_function_type)(  ) const;
            chemicalPotential_function_type chemicalPotential_function_value( &::SireMove::Moves::chemicalPotential );
            
            Moves_exposer.def( 
                "chemicalPotential"
                , chemicalPotential_function_value
                , bp::release_gil_policy()
                , "Return the constant chemical potential that these moves sample\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::clearStatistics
        
            typedef void ( ::SireMove::Moves::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::Moves::clearStatistics );
            
            Moves_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::clearTiming
        
            typedef void ( ::SireMove::Moves::*clearTiming_function_type)(  ) ;
            clearTiming_function_type clearTiming_function_value( &::SireMove::Moves::clearTiming );
            
            Moves_exposer.def( 
                "clearTiming"
                , clearTiming_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::count
        
            typedef int ( ::SireMove::Moves::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMove::Moves::count );
            
            Moves_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of different move types in this set of Moves\n(this is the maximum number of Move objects that can be iterated over)" );
        
        }
        { //::SireMove::Moves::energy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::Moves::*energy_function_type)( ::SireSystem::System & ) const;
            energy_function_type energy_function_value( &::SireMove::Moves::energy );
            
            Moves_exposer.def( 
                "energy"
                , energy_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Return the energy of the system that is seen by these moves" );
        
        }
        { //::SireMove::Moves::energyComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::Moves::*energyComponent_function_type)(  ) const;
            energyComponent_function_type energyComponent_function_value( &::SireMove::Moves::energyComponent );
            
            Moves_exposer.def( 
                "energyComponent"
                , energyComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMove::Moves::ensemble
        
            typedef ::SireMove::Ensemble ( ::SireMove::Moves::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::Moves::ensemble );
            
            Moves_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::release_gil_policy()
                , "Return the ensemble that will be sampled by these moves" );
        
        }
        { //::SireMove::Moves::fugacity
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::Moves::*fugacity_function_type)(  ) const;
            fugacity_function_type fugacity_function_value( &::SireMove::Moves::fugacity );
            
            Moves_exposer.def( 
                "fugacity"
                , fugacity_function_value
                , bp::release_gil_policy()
                , "Return the constant fugacity that these moves sample\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::isConstantChemicalPotential
        
            typedef bool ( ::SireMove::Moves::*isConstantChemicalPotential_function_type)(  ) const;
            isConstantChemicalPotential_function_type isConstantChemicalPotential_function_value( &::SireMove::Moves::isConstantChemicalPotential );
            
            Moves_exposer.def( 
                "isConstantChemicalPotential"
                , isConstantChemicalPotential_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the chemical potential constant" );
        
        }
        { //::SireMove::Moves::isConstantEnergy
        
            typedef bool ( ::SireMove::Moves::*isConstantEnergy_function_type)(  ) const;
            isConstantEnergy_function_type isConstantEnergy_function_value( &::SireMove::Moves::isConstantEnergy );
            
            Moves_exposer.def( 
                "isConstantEnergy"
                , isConstantEnergy_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keeps the total energy constant" );
        
        }
        { //::SireMove::Moves::isConstantFugacity
        
            typedef bool ( ::SireMove::Moves::*isConstantFugacity_function_type)(  ) const;
            isConstantFugacity_function_type isConstantFugacity_function_value( &::SireMove::Moves::isConstantFugacity );
            
            Moves_exposer.def( 
                "isConstantFugacity"
                , isConstantFugacity_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the fugacity (related to the chemical\npotential) constant" );
        
        }
        { //::SireMove::Moves::isConstantLambda
        
            typedef bool ( ::SireMove::Moves::*isConstantLambda_function_type)( ::SireCAS::Symbol const & ) const;
            isConstantLambda_function_type isConstantLambda_function_value( &::SireMove::Moves::isConstantLambda );
            
            Moves_exposer.def( 
                "isConstantLambda"
                , isConstantLambda_function_value
                , ( bp::arg("lam") )
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the symbol lam constant" );
        
        }
        { //::SireMove::Moves::isConstantPressure
        
            typedef bool ( ::SireMove::Moves::*isConstantPressure_function_type)(  ) const;
            isConstantPressure_function_type isConstantPressure_function_value( &::SireMove::Moves::isConstantPressure );
            
            Moves_exposer.def( 
                "isConstantPressure"
                , isConstantPressure_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the pressure constant" );
        
        }
        { //::SireMove::Moves::isConstantTemperature
        
            typedef bool ( ::SireMove::Moves::*isConstantTemperature_function_type)(  ) const;
            isConstantTemperature_function_type isConstantTemperature_function_value( &::SireMove::Moves::isConstantTemperature );
            
            Moves_exposer.def( 
                "isConstantTemperature"
                , isConstantTemperature_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the temperature constant" );
        
        }
        { //::SireMove::Moves::isConstantVolume
        
            typedef bool ( ::SireMove::Moves::*isConstantVolume_function_type)(  ) const;
            isConstantVolume_function_type isConstantVolume_function_value( &::SireMove::Moves::isConstantVolume );
            
            Moves_exposer.def( 
                "isConstantVolume"
                , isConstantVolume_function_value
                , bp::release_gil_policy()
                , "Return whether or not these moves keep the volume constant" );
        
        }
        { //::SireMove::Moves::move
        
            typedef ::SireSystem::System ( ::SireMove::Moves::*move_function_type)( ::SireSystem::System const &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::Moves::move );
            
            Moves_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves")=(int)(1), bp::arg("record_stats")=(bool)(false) )
                , "" );
        
        }
        { //::SireMove::Moves::moves
        
            typedef ::QList< SireBase::PropPtr< SireMove::Move > > ( ::SireMove::Moves::*moves_function_type)(  ) const;
            moves_function_type moves_function_value( &::SireMove::Moves::moves );
            
            Moves_exposer.def( 
                "moves"
                , moves_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::nMoveTypes
        
            typedef int ( ::SireMove::Moves::*nMoveTypes_function_type)(  ) const;
            nMoveTypes_function_type nMoveTypes_function_value( &::SireMove::Moves::nMoveTypes );
            
            Moves_exposer.def( 
                "nMoveTypes"
                , nMoveTypes_function_value
                , bp::release_gil_policy()
                , "Return the number of different move types in this set of Moves\n(this is the maximum number of Move objects that can be iterated over)" );
        
        }
        { //::SireMove::Moves::nMoves
        
            typedef int ( ::SireMove::Moves::*nMoves_function_type)(  ) const;
            nMoves_function_type nMoves_function_value( &::SireMove::Moves::nMoves );
            
            Moves_exposer.def( 
                "nMoves"
                , nMoves_function_value
                , bp::release_gil_policy()
                , "Return the total number of moves performed using the moves\nin this set" );
        
        }
        { //::SireMove::Moves::null
        
            typedef ::SireMove::SameMoves const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMove::Moves::null );
            
            Moves_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMove::Moves::operator[]
        
            typedef ::SireMove::Move const & ( ::SireMove::Moves::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMove::Moves::operator[] );
            
            Moves_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "" );
        
        }
        { //::SireMove::Moves::pressure
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::Moves::*pressure_function_type)(  ) const;
            pressure_function_type pressure_function_value( &::SireMove::Moves::pressure );
            
            Moves_exposer.def( 
                "pressure"
                , pressure_function_value
                , bp::release_gil_policy()
                , "Return the constant pressure that these moves sample\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::setAcceptableDelta
        
            typedef void ( ::SireMove::Moves::*setAcceptableDelta_function_type)( ::SireUnits::Dimension::MolarEnergy ) ;
            setAcceptableDelta_function_type setAcceptableDelta_function_value( &::SireMove::Moves::setAcceptableDelta );
            
            Moves_exposer.def( 
                "setAcceptableDelta"
                , setAcceptableDelta_function_value
                , ( bp::arg("delta") )
                , bp::release_gil_policy()
                , "Set the acceptable level of drift in the running total from the recalculated total energy" );
        
        }
        { //::SireMove::Moves::setCheckRunningTotal
        
            typedef void ( ::SireMove::Moves::*setCheckRunningTotal_function_type)( bool ) ;
            setCheckRunningTotal_function_type setCheckRunningTotal_function_value( &::SireMove::Moves::setCheckRunningTotal );
            
            Moves_exposer.def( 
                "setCheckRunningTotal"
                , setCheckRunningTotal_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Switch on or off the checking of the running total energy" );
        
        }
        { //::SireMove::Moves::setChemicalPotential
        
            typedef void ( ::SireMove::Moves::*setChemicalPotential_function_type)( ::SireUnits::Dimension::MolarEnergy const & ) ;
            setChemicalPotential_function_type setChemicalPotential_function_value( &::SireMove::Moves::setChemicalPotential );
            
            Moves_exposer.def( 
                "setChemicalPotential"
                , setChemicalPotential_function_value
                , ( bp::arg("chemical_potential") )
                , bp::release_gil_policy()
                , "Set the chemical potential that these constant chemical potential moves sample\nto chemical_potential\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::setEnergyComponent
        
            typedef void ( ::SireMove::Moves::*setEnergyComponent_function_type)( ::SireCAS::Symbol const & ) ;
            setEnergyComponent_function_type setEnergyComponent_function_value( &::SireMove::Moves::setEnergyComponent );
            
            Moves_exposer.def( 
                "setEnergyComponent"
                , setEnergyComponent_function_value
                , ( bp::arg("nrg_component") )
                , bp::release_gil_policy()
                , "Set the energy component used by these moves - by default\nthis will raise an error as setting the energy component is\nnot supported\nThrow: SireError::unsupported\n" );
        
        }
        { //::SireMove::Moves::setFugacity
        
            typedef void ( ::SireMove::Moves::*setFugacity_function_type)( ::SireUnits::Dimension::Pressure const & ) ;
            setFugacity_function_type setFugacity_function_value( &::SireMove::Moves::setFugacity );
            
            Moves_exposer.def( 
                "setFugacity"
                , setFugacity_function_value
                , ( bp::arg("fugacity") )
                , bp::release_gil_policy()
                , "Set the fugacity that these constant fugacity moves sample\nto fugacity\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::setGenerator
        
            typedef void ( ::SireMove::Moves::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::Moves::setGenerator );
            
            Moves_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::setPressure
        
            typedef void ( ::SireMove::Moves::*setPressure_function_type)( ::SireUnits::Dimension::Pressure const & ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::Moves::setPressure );
            
            Moves_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("pressure") )
                , bp::release_gil_policy()
                , "Set the pressure that these constant pressure moves sample\nto pressure\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::setSpaceProperty
        
            typedef void ( ::SireMove::Moves::*setSpaceProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setSpaceProperty_function_type setSpaceProperty_function_value( &::SireMove::Moves::setSpaceProperty );
            
            Moves_exposer.def( 
                "setSpaceProperty"
                , setSpaceProperty_function_value
                , ( bp::arg("space_property") )
                , bp::release_gil_policy()
                , "Set the space property used by these moves - by default\nthis will raise an error as setting the space property is\nnot supported\nThrow: SireError::unsupported\n" );
        
        }
        { //::SireMove::Moves::setTemperature
        
            typedef void ( ::SireMove::Moves::*setTemperature_function_type)( ::SireUnits::Dimension::Temperature const & ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::Moves::setTemperature );
            
            Moves_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("temperature") )
                , bp::release_gil_policy()
                , "Set the temperature that these constant temperature moves sample\nto temperature\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::size
        
            typedef int ( ::SireMove::Moves::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMove::Moves::size );
            
            Moves_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "Return the number of different move types in this set of Moves\n(this is the maximum number of Move objects that can be iterated over)" );
        
        }
        { //::SireMove::Moves::spaceProperty
        
            typedef ::SireBase::PropertyName const & ( ::SireMove::Moves::*spaceProperty_function_type)(  ) const;
            spaceProperty_function_type spaceProperty_function_value( &::SireMove::Moves::spaceProperty );
            
            Moves_exposer.def( 
                "spaceProperty"
                , spaceProperty_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMove::Moves::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::Moves::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireMove::Moves::temperature );
            
            Moves_exposer.def( 
                "temperature"
                , temperature_function_value
                , bp::release_gil_policy()
                , "Return the constant temperature that these moves sample\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Moves::timing
        
            typedef ::SireUnits::Dimension::Time ( ::SireMove::Moves::*timing_function_type)( int ) const;
            timing_function_type timing_function_value( &::SireMove::Moves::timing );
            
            Moves_exposer.def( 
                "timing"
                , timing_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return the average time per move for the ith Move object in this set\nof Moves. This returns 0 if no timing data is available" );
        
        }
        { //::SireMove::Moves::timing
        
            typedef ::QList< SireUnits::Dimension::PhysUnit< 0, 0, 1, 0, 0, 0, 0 > > ( ::SireMove::Moves::*timing_function_type)(  ) const;
            timing_function_type timing_function_value( &::SireMove::Moves::timing );
            
            Moves_exposer.def( 
                "timing"
                , timing_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::toString
        
            typedef ::QString ( ::SireMove::Moves::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::Moves::toString );
            
            Moves_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Moves::typeName );
            
            Moves_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMove::Moves::volume
        
            typedef ::SireUnits::Dimension::Volume ( ::SireMove::Moves::*volume_function_type)( ::SireSystem::System const & ) const;
            volume_function_type volume_function_value( &::SireMove::Moves::volume );
            
            Moves_exposer.def( 
                "volume"
                , volume_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Return the volume of the system as it is seen by these moves" );
        
        }
        Moves_exposer.staticmethod( "null" );
        Moves_exposer.staticmethod( "typeName" );
        Moves_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Moves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Moves_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Moves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Moves_exposer.def_pickle(sire_pickle_suite< ::SireMove::Moves >());
        Moves_exposer.def( "__str__", &__str__< ::SireMove::Moves > );
        Moves_exposer.def( "__repr__", &__str__< ::SireMove::Moves > );
        Moves_exposer.def( "__len__", &__len_size< ::SireMove::Moves > );
    }

}
