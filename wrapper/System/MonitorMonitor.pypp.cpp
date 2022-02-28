// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MonitorMonitor.pypp.hpp"

namespace bp = boost::python;

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "monitormonitor.h"

#include "monitormonitor.h"

SireSystem::MonitorMonitor __copy__(const SireSystem::MonitorMonitor &other){ return SireSystem::MonitorMonitor(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_MonitorMonitor_class(){

    { //::SireSystem::MonitorMonitor
        typedef bp::class_< SireSystem::MonitorMonitor, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > MonitorMonitor_exposer_t;
        MonitorMonitor_exposer_t MonitorMonitor_exposer = MonitorMonitor_exposer_t( "MonitorMonitor", "This is a monitor that can be used to monitor other monitors.\nIt is useful to use to save the current state of a monitor,\ne.g. snap-shotting the average during the simulation, or\nsnap-shotting the PDB. Equally, it can be used to transfer\nmonitor data between systems, or from a system to a supra-system.\n\nThe MonitorMonitor can be non-destructive (the original monitor\nis just copied), or destructive (the original monitor is copied,\nthen cleared), or annihilative (the original monitor is copied,\nthen completely removed from the system)\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope MonitorMonitor_scope( MonitorMonitor_exposer );
        MonitorMonitor_exposer.def( bp::init< SireSystem::MonitorID const &, bp::optional< bool, bool > >(( bp::arg("id"), bp::arg("clear_original")=(bool)(false), bp::arg("remove_original")=(bool)(false) ), "Construct to monitor the SystemMonitor identified by id,\noptionally clearing the original monitor if clear_original is\ntrue, and optionally removing the original monitor\nif remove_original is true") );
        MonitorMonitor_exposer.def( bp::init< SireSystem::MonitorMonitor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::MonitorMonitor::at
        
            typedef ::SireSystem::SystemMonitor const & ( ::SireSystem::MonitorMonitor::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireSystem::MonitorMonitor::at );
            
            MonitorMonitor_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ith state monitored\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireSystem::MonitorMonitor::clearOriginal
        
            typedef bool ( ::SireSystem::MonitorMonitor::*clearOriginal_function_type)(  ) const;
            clearOriginal_function_type clearOriginal_function_value( &::SireSystem::MonitorMonitor::clearOriginal );
            
            MonitorMonitor_exposer.def( 
                "clearOriginal"
                , clearOriginal_function_value
                , "Return whether or not this MonitorMonitor will clear the\noriginal monitor whenever it takes a copy" );
        
        }
        { //::SireSystem::MonitorMonitor::clearStatistics
        
            typedef void ( ::SireSystem::MonitorMonitor::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::MonitorMonitor::clearStatistics );
            
            MonitorMonitor_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , "Completely clear statistics" );
        
        }
        { //::SireSystem::MonitorMonitor::count
        
            typedef int ( ::SireSystem::MonitorMonitor::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireSystem::MonitorMonitor::count );
            
            MonitorMonitor_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of states monitored so far" );
        
        }
        { //::SireSystem::MonitorMonitor::monitor
        
            typedef void ( ::SireSystem::MonitorMonitor::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::MonitorMonitor::monitor );
            
            MonitorMonitor_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , "Monitor the passed system" );
        
        }
        { //::SireSystem::MonitorMonitor::nStates
        
            typedef int ( ::SireSystem::MonitorMonitor::*nStates_function_type)(  ) const;
            nStates_function_type nStates_function_value( &::SireSystem::MonitorMonitor::nStates );
            
            MonitorMonitor_exposer.def( 
                "nStates"
                , nStates_function_value
                , "Return the number of states monitored so far" );
        
        }
        MonitorMonitor_exposer.def( bp::self != bp::self );
        { //::SireSystem::MonitorMonitor::operator=
        
            typedef ::SireSystem::MonitorMonitor & ( ::SireSystem::MonitorMonitor::*assign_function_type)( ::SireSystem::MonitorMonitor const & ) ;
            assign_function_type assign_function_value( &::SireSystem::MonitorMonitor::operator= );
            
            MonitorMonitor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MonitorMonitor_exposer.def( bp::self == bp::self );
        { //::SireSystem::MonitorMonitor::operator[]
        
            typedef ::SireSystem::SystemMonitor const & ( ::SireSystem::MonitorMonitor::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireSystem::MonitorMonitor::operator[] );
            
            MonitorMonitor_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireSystem::MonitorMonitor::removeOriginal
        
            typedef bool ( ::SireSystem::MonitorMonitor::*removeOriginal_function_type)(  ) const;
            removeOriginal_function_type removeOriginal_function_value( &::SireSystem::MonitorMonitor::removeOriginal );
            
            MonitorMonitor_exposer.def( 
                "removeOriginal"
                , removeOriginal_function_value
                , "Return whether or not this MonitorMonitor will remove\nthe original monitor from the system whenever it\ntakes a copy" );
        
        }
        { //::SireSystem::MonitorMonitor::setClearOriginal
        
            typedef void ( ::SireSystem::MonitorMonitor::*setClearOriginal_function_type)( bool ) ;
            setClearOriginal_function_type setClearOriginal_function_value( &::SireSystem::MonitorMonitor::setClearOriginal );
            
            MonitorMonitor_exposer.def( 
                "setClearOriginal"
                , setClearOriginal_function_value
                , ( bp::arg("clear") )
                , "Set whether or not to clear the statistics of the original\nmonitor when this monitor takes the copy" );
        
        }
        { //::SireSystem::MonitorMonitor::setRemoveOriginal
        
            typedef void ( ::SireSystem::MonitorMonitor::*setRemoveOriginal_function_type)( bool ) ;
            setRemoveOriginal_function_type setRemoveOriginal_function_value( &::SireSystem::MonitorMonitor::setRemoveOriginal );
            
            MonitorMonitor_exposer.def( 
                "setRemoveOriginal"
                , setRemoveOriginal_function_value
                , ( bp::arg("remove") )
                , "Set whether or not to remove the original monitor when\nthis monitor takes a copy (effectively thus moving\nthe monitor from the system to this MonitorMonitor)" );
        
        }
        { //::SireSystem::MonitorMonitor::size
        
            typedef int ( ::SireSystem::MonitorMonitor::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireSystem::MonitorMonitor::size );
            
            MonitorMonitor_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of states monitored so far" );
        
        }
        { //::SireSystem::MonitorMonitor::states
        
            typedef ::QList< SireBase::PropPtr< SireSystem::SystemMonitor > > const & ( ::SireSystem::MonitorMonitor::*states_function_type)(  ) const;
            states_function_type states_function_value( &::SireSystem::MonitorMonitor::states );
            
            MonitorMonitor_exposer.def( 
                "states"
                , states_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the states monitored" );
        
        }
        { //::SireSystem::MonitorMonitor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::MonitorMonitor::typeName );
            
            MonitorMonitor_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        MonitorMonitor_exposer.staticmethod( "typeName" );
        MonitorMonitor_exposer.def( "__copy__", &__copy__);
        MonitorMonitor_exposer.def( "__deepcopy__", &__copy__);
        MonitorMonitor_exposer.def( "clone", &__copy__);
        MonitorMonitor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::MonitorMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorMonitor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::MonitorMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorMonitor_exposer.def_pickle(sire_pickle_suite< ::SireSystem::MonitorMonitor >());
        MonitorMonitor_exposer.def( "__str__", &__str__< ::SireSystem::MonitorMonitor > );
        MonitorMonitor_exposer.def( "__repr__", &__str__< ::SireSystem::MonitorMonitor > );
        MonitorMonitor_exposer.def( "__len__", &__len_size< ::SireSystem::MonitorMonitor > );
    }

}
