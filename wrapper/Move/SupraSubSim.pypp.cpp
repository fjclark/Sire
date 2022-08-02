// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SupraSubSim.pypp.hpp"

namespace bp = boost::python;

#include "SireCluster/nodes.h"

#include "SireError/errors.h"

#include "suprasubmoves.h"

#include "suprasubsim.h"

#include "suprasubsim.h"

SireMove::SupraSubSim __copy__(const SireMove::SupraSubSim &other){ return SireMove::SupraSubSim(other); }

const char* pvt_get_name(const SireMove::SupraSubSim&){ return "SireMove::SupraSubSim";}

#include "Helpers/release_gil_policy.hpp"

void register_SupraSubSim_class(){

    { //::SireMove::SupraSubSim
        typedef bp::class_< SireMove::SupraSubSim > SupraSubSim_exposer_t;
        SupraSubSim_exposer_t SupraSubSim_exposer = SupraSubSim_exposer_t( "SupraSubSim", "This class is used to start and manage an active\nsub-simulation of a supra-simulation.\n\nA supra-simulation consists of a collection of\nsupra-moves that are applied to a supra-system,\nwhile a supra-sub-simulation consists of a collection\nof sub-moves that are applied to a sub-system\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope SupraSubSim_scope( SupraSubSim_exposer );
        SupraSubSim_exposer.def( bp::init< SireMove::SupraSubSim const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::SupraSubSim::abort
        
            typedef void ( ::SireMove::SupraSubSim::*abort_function_type)(  ) ;
            abort_function_type abort_function_value( &::SireMove::SupraSubSim::abort );
            
            SupraSubSim_exposer.def( 
                "abort"
                , abort_function_value
                , bp::release_gil_policy()
                , "Abort the running simulation" );
        
        }
        { //::SireMove::SupraSubSim::hasFinished
        
            typedef bool ( ::SireMove::SupraSubSim::*hasFinished_function_type)(  ) ;
            hasFinished_function_type hasFinished_function_value( &::SireMove::SupraSubSim::hasFinished );
            
            SupraSubSim_exposer.def( 
                "hasFinished"
                , hasFinished_function_value
                , bp::release_gil_policy()
                , "Return whether or not the simulation has finished\n(completed all of the moves)" );
        
        }
        { //::SireMove::SupraSubSim::initialMoves
        
            typedef ::SireMove::SupraSubMovesPtr ( ::SireMove::SupraSubSim::*initialMoves_function_type)(  ) ;
            initialMoves_function_type initialMoves_function_value( &::SireMove::SupraSubSim::initialMoves );
            
            SupraSubSim_exposer.def( 
                "initialMoves"
                , initialMoves_function_value
                , bp::release_gil_policy()
                , "Return the Moves in the state they were in before the simulation started" );
        
        }
        { //::SireMove::SupraSubSim::initialSystem
        
            typedef ::SireMove::SupraSubSystemPtr ( ::SireMove::SupraSubSim::*initialSystem_function_type)(  ) ;
            initialSystem_function_type initialSystem_function_value( &::SireMove::SupraSubSim::initialSystem );
            
            SupraSubSim_exposer.def( 
                "initialSystem"
                , initialSystem_function_value
                , bp::release_gil_policy()
                , "Return the sub-system in the state it was in before the simulation started" );
        
        }
        { //::SireMove::SupraSubSim::input
        
            typedef ::SireMove::SupraSubSimPacket ( ::SireMove::SupraSubSim::*input_function_type)(  ) ;
            input_function_type input_function_value( &::SireMove::SupraSubSim::input );
            
            SupraSubSim_exposer.def( 
                "input"
                , input_function_value
                , bp::release_gil_policy()
                , "Return the initial input simulation WorkPacket" );
        
        }
        { //::SireMove::SupraSubSim::interimMoves
        
            typedef ::SireMove::SupraSubMovesPtr ( ::SireMove::SupraSubSim::*interimMoves_function_type)(  ) ;
            interimMoves_function_type interimMoves_function_value( &::SireMove::SupraSubSim::interimMoves );
            
            SupraSubSim_exposer.def( 
                "interimMoves"
                , interimMoves_function_value
                , bp::release_gil_policy()
                , "Return the current state of the moves (updated while the simulation\nis running). This will throw an exception if the system hits an\nerror state" );
        
        }
        { //::SireMove::SupraSubSim::interimResult
        
            typedef ::SireMove::SupraSubSimPacket ( ::SireMove::SupraSubSim::*interimResult_function_type)(  ) ;
            interimResult_function_type interimResult_function_value( &::SireMove::SupraSubSim::interimResult );
            
            SupraSubSim_exposer.def( 
                "interimResult"
                , interimResult_function_value
                , bp::release_gil_policy()
                , "Return the simulation WorkPacket from an intermediate point along\nthe simulation. This will throw an error if the simulation is in an\nerror state, and the initial packet if the simulation\nwas aborted" );
        
        }
        { //::SireMove::SupraSubSim::interimSystem
        
            typedef ::SireMove::SupraSubSystemPtr ( ::SireMove::SupraSubSim::*interimSystem_function_type)(  ) ;
            interimSystem_function_type interimSystem_function_value( &::SireMove::SupraSubSim::interimSystem );
            
            SupraSubSim_exposer.def( 
                "interimSystem"
                , interimSystem_function_value
                , bp::release_gil_policy()
                , "Return the current state of the System (updated while the simulation\nis running). This will throw an exception if the system hits an\nerror state" );
        
        }
        { //::SireMove::SupraSubSim::isError
        
            typedef bool ( ::SireMove::SupraSubSim::*isError_function_type)(  ) ;
            isError_function_type isError_function_value( &::SireMove::SupraSubSim::isError );
            
            SupraSubSim_exposer.def( 
                "isError"
                , isError_function_value
                , bp::release_gil_policy()
                , "Return whether or not this simulation is in an error state" );
        
        }
        { //::SireMove::SupraSubSim::isRunning
        
            typedef bool ( ::SireMove::SupraSubSim::*isRunning_function_type)(  ) ;
            isRunning_function_type isRunning_function_value( &::SireMove::SupraSubSim::isRunning );
            
            SupraSubSim_exposer.def( 
                "isRunning"
                , isRunning_function_value
                , bp::release_gil_policy()
                , "Return whether or not this simulation is running" );
        
        }
        { //::SireMove::SupraSubSim::moves
        
            typedef ::SireMove::SupraSubMovesPtr ( ::SireMove::SupraSubSim::*moves_function_type)(  ) ;
            moves_function_type moves_function_value( &::SireMove::SupraSubSim::moves );
            
            SupraSubSim_exposer.def( 
                "moves"
                , moves_function_value
                , bp::release_gil_policy()
                , "Return the final state of the moves after the simulation. This\nblocks until the simulation has finished and will throw an\nexception if the system hits an error state" );
        
        }
        SupraSubSim_exposer.def( bp::self != bp::self );
        { //::SireMove::SupraSubSim::operator=
        
            typedef ::SireMove::SupraSubSim & ( ::SireMove::SupraSubSim::*assign_function_type)( ::SireMove::SupraSubSim const & ) ;
            assign_function_type assign_function_value( &::SireMove::SupraSubSim::operator= );
            
            SupraSubSim_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SupraSubSim_exposer.def( bp::self == bp::self );
        { //::SireMove::SupraSubSim::progress
        
            typedef float ( ::SireMove::SupraSubSim::*progress_function_type)(  ) ;
            progress_function_type progress_function_value( &::SireMove::SupraSubSim::progress );
            
            SupraSubSim_exposer.def( 
                "progress"
                , progress_function_value
                , bp::release_gil_policy()
                , "Return the progress of the simulation (as a percentage)" );
        
        }
        { //::SireMove::SupraSubSim::result
        
            typedef ::SireMove::SupraSubSimPacket ( ::SireMove::SupraSubSim::*result_function_type)(  ) ;
            result_function_type result_function_value( &::SireMove::SupraSubSim::result );
            
            SupraSubSim_exposer.def( 
                "result"
                , result_function_value
                , bp::release_gil_policy()
                , "Return the final result of the simulation. This blocks until\nthe simulation has stopped, and will throw an exception if the\nsimulation is in an error state. This returns the initial\nsimulation WorkPacket if the simulation was aborted" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireMove::SupraSubSystem const &,::SireMove::SupraSubMoves const &,int,bool );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("system"), bp::arg("moves"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Run the sub-system simulation consisting of nmoves moves from moves\non the sub-system system, recording statistics if record_stats\nis true. Run the simulation in the current thread" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireMove::SupraSubSystem const &,::SireMove::SupraSubMove const &,int,bool );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("system"), bp::arg("move"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Run the sub-system simulation consisting of nmoves moves from moves\non the sub-system system, recording statistics if record_stats\nis true. Run the simulation in the current thread" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireMove::SupraSubSimPacket const & );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("simpacket") )
                , bp::release_gil_policy()
                , "Run the sub-system simulation described\nin simpacket in the current thread" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireCluster::Node &,::SireMove::SupraSubSystem const &,::SireMove::SupraSubMoves const &,int,bool );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("node"), bp::arg("system"), bp::arg("moves"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Run the sub-system simulation consisting of nmoves moves from moves\non the sub-system system, recording statistics if record_stats\nis true. Run the simulation on the node node, returning a handle to\nthe running simulation" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireCluster::Node &,::SireMove::SupraSubSystem const &,::SireMove::SupraSubMove const &,int,bool );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("node"), bp::arg("system"), bp::arg("move"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Run the sub-system simulation consisting of nmoves moves from moves\non the sub-system system, recording statistics if record_stats\nis true. Run the simulation on the node node, returning a handle to\nthe running simulation" );
        
        }
        { //::SireMove::SupraSubSim::run
        
            typedef ::SireMove::SupraSubSim ( *run_function_type )( ::SireCluster::Node &,::SireMove::SupraSubSimPacket const & );
            run_function_type run_function_value( &::SireMove::SupraSubSim::run );
            
            SupraSubSim_exposer.def( 
                "run"
                , run_function_value
                , ( bp::arg("node"), bp::arg("simpacket") )
                , bp::release_gil_policy()
                , "Run the sub-system simulation described in simpacket on the node node\nand return a handle to the running simulation" );
        
        }
        { //::SireMove::SupraSubSim::stop
        
            typedef void ( ::SireMove::SupraSubSim::*stop_function_type)(  ) ;
            stop_function_type stop_function_value( &::SireMove::SupraSubSim::stop );
            
            SupraSubSim_exposer.def( 
                "stop"
                , stop_function_value
                , bp::release_gil_policy()
                , "Stop the running simulation" );
        
        }
        { //::SireMove::SupraSubSim::system
        
            typedef ::SireMove::SupraSubSystemPtr ( ::SireMove::SupraSubSim::*system_function_type)(  ) ;
            system_function_type system_function_value( &::SireMove::SupraSubSim::system );
            
            SupraSubSim_exposer.def( 
                "system"
                , system_function_value
                , bp::release_gil_policy()
                , "Return the final state of the system after the simulation. This\nblocks until the simulation has finished and will throw an\nexception if the system hits an error state" );
        
        }
        { //::SireMove::SupraSubSim::throwError
        
            typedef void ( ::SireMove::SupraSubSim::*throwError_function_type)(  ) ;
            throwError_function_type throwError_function_value( &::SireMove::SupraSubSim::throwError );
            
            SupraSubSim_exposer.def( 
                "throwError"
                , throwError_function_value
                , bp::release_gil_policy()
                , "Throw any error associated with this simulation - this does\nnothing if we are not in an error state" );
        
        }
        { //::SireMove::SupraSubSim::wait
        
            typedef void ( ::SireMove::SupraSubSim::*wait_function_type)(  ) ;
            wait_function_type wait_function_value( &::SireMove::SupraSubSim::wait );
            
            SupraSubSim_exposer.def( 
                "wait"
                , wait_function_value
                , bp::release_gil_policy()
                , "Wait for the simulation to complete" );
        
        }
        { //::SireMove::SupraSubSim::wait
        
            typedef bool ( ::SireMove::SupraSubSim::*wait_function_type)( int ) ;
            wait_function_type wait_function_value( &::SireMove::SupraSubSim::wait );
            
            SupraSubSim_exposer.def( 
                "wait"
                , wait_function_value
                , ( bp::arg("timeout") )
                , bp::release_gil_policy()
                , "Wait for the simulation to stop running, or for timeout\nmilliseconds to pass, whichever comes soonest. This returns\nwhether or not the simulation has stopped" );
        
        }
        { //::SireMove::SupraSubSim::wasAborted
        
            typedef bool ( ::SireMove::SupraSubSim::*wasAborted_function_type)(  ) ;
            wasAborted_function_type wasAborted_function_value( &::SireMove::SupraSubSim::wasAborted );
            
            SupraSubSim_exposer.def( 
                "wasAborted"
                , wasAborted_function_value
                , bp::release_gil_policy()
                , "Return whether or not the simulation was aborted" );
        
        }
        { //::SireMove::SupraSubSim::wasStopped
        
            typedef bool ( ::SireMove::SupraSubSim::*wasStopped_function_type)(  ) ;
            wasStopped_function_type wasStopped_function_value( &::SireMove::SupraSubSim::wasStopped );
            
            SupraSubSim_exposer.def( 
                "wasStopped"
                , wasStopped_function_value
                , bp::release_gil_policy()
                , "Return whether or not the simulation was stopped" );
        
        }
        SupraSubSim_exposer.staticmethod( "run" );
        SupraSubSim_exposer.def( "__copy__", &__copy__);
        SupraSubSim_exposer.def( "__deepcopy__", &__copy__);
        SupraSubSim_exposer.def( "clone", &__copy__);
        SupraSubSim_exposer.def( "__str__", &pvt_get_name);
        SupraSubSim_exposer.def( "__repr__", &pvt_get_name);
    }

}
