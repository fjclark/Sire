// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "WorkPacket.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "workpacket.h"

#include <QDebug>

#include <QTextStream>

#include "workpacket.h"

SireCluster::WorkPacket __copy__(const SireCluster::WorkPacket &other){ return SireCluster::WorkPacket(other); }

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireCluster::WorkPacket&){ return "SireCluster::WorkPacket";}

void register_WorkPacket_class(){

    { //::SireCluster::WorkPacket
        typedef bp::class_< SireCluster::WorkPacket > WorkPacket_exposer_t;
        WorkPacket_exposer_t WorkPacket_exposer = WorkPacket_exposer_t( "WorkPacket", bp::init< >() );
        bp::scope WorkPacket_scope( WorkPacket_exposer );
        WorkPacket_exposer.def( bp::init< SireCluster::WorkPacketBase const & >(( bp::arg("work") )) );
        WorkPacket_exposer.def( bp::init< SireCluster::WorkPacket const & >(( bp::arg("other") )) );
        { //::SireCluster::WorkPacket::abort
        
            typedef void ( ::SireCluster::WorkPacket::*abort_function_type)(  ) ;
            abort_function_type abort_function_value( &::SireCluster::WorkPacket::abort );
            
            WorkPacket_exposer.def( 
                "abort"
                , abort_function_value );
        
        }
        { //::SireCluster::WorkPacket::base
        
            typedef ::SireCluster::WorkPacketBase const & ( ::SireCluster::WorkPacket::*base_function_type)(  ) const;
            base_function_type base_function_value( &::SireCluster::WorkPacket::base );
            
            WorkPacket_exposer.def( 
                "base"
                , base_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireCluster::WorkPacket::hasFinished
        
            typedef bool ( ::SireCluster::WorkPacket::*hasFinished_function_type)(  ) const;
            hasFinished_function_type hasFinished_function_value( &::SireCluster::WorkPacket::hasFinished );
            
            WorkPacket_exposer.def( 
                "hasFinished"
                , hasFinished_function_value );
        
        }
        { //::SireCluster::WorkPacket::isError
        
            typedef bool ( ::SireCluster::WorkPacket::*isError_function_type)(  ) const;
            isError_function_type isError_function_value( &::SireCluster::WorkPacket::isError );
            
            WorkPacket_exposer.def( 
                "isError"
                , isError_function_value );
        
        }
        { //::SireCluster::WorkPacket::isNull
        
            typedef bool ( ::SireCluster::WorkPacket::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireCluster::WorkPacket::isNull );
            
            WorkPacket_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireCluster::WorkPacket::operator=
        
            typedef ::SireCluster::WorkPacket & ( ::SireCluster::WorkPacket::*assign_function_type)( ::SireCluster::WorkPacket const & ) ;
            assign_function_type assign_function_value( &::SireCluster::WorkPacket::operator= );
            
            WorkPacket_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireCluster::WorkPacket::pack
        
            typedef ::QByteArray ( ::SireCluster::WorkPacket::*pack_function_type)(  ) const;
            pack_function_type pack_function_value( &::SireCluster::WorkPacket::pack );
            
            WorkPacket_exposer.def( 
                "pack"
                , pack_function_value );
        
        }
        { //::SireCluster::WorkPacket::progress
        
            typedef float ( ::SireCluster::WorkPacket::*progress_function_type)(  ) const;
            progress_function_type progress_function_value( &::SireCluster::WorkPacket::progress );
            
            WorkPacket_exposer.def( 
                "progress"
                , progress_function_value );
        
        }
        { //::SireCluster::WorkPacket::runChunk
        
            typedef void ( ::SireCluster::WorkPacket::*runChunk_function_type)(  ) ;
            runChunk_function_type runChunk_function_value( &::SireCluster::WorkPacket::runChunk );
            
            WorkPacket_exposer.def( 
                "runChunk"
                , runChunk_function_value );
        
        }
        { //::SireCluster::WorkPacket::shouldPack
        
            typedef bool ( ::SireCluster::WorkPacket::*shouldPack_function_type)(  ) const;
            shouldPack_function_type shouldPack_function_value( &::SireCluster::WorkPacket::shouldPack );
            
            WorkPacket_exposer.def( 
                "shouldPack"
                , shouldPack_function_value );
        
        }
        { //::SireCluster::WorkPacket::throwError
        
            typedef void ( ::SireCluster::WorkPacket::*throwError_function_type)(  ) const;
            throwError_function_type throwError_function_value( &::SireCluster::WorkPacket::throwError );
            
            WorkPacket_exposer.def( 
                "throwError"
                , throwError_function_value );
        
        }
        { //::SireCluster::WorkPacket::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCluster::WorkPacket::typeName );
            
            WorkPacket_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCluster::WorkPacket::unpack
        
            typedef ::SireCluster::WorkPacket ( *unpack_function_type )( ::QByteArray const & );
            unpack_function_type unpack_function_value( &::SireCluster::WorkPacket::unpack );
            
            WorkPacket_exposer.def( 
                "unpack"
                , unpack_function_value
                , ( bp::arg("data") ) );
        
        }
        { //::SireCluster::WorkPacket::wasAborted
        
            typedef bool ( ::SireCluster::WorkPacket::*wasAborted_function_type)(  ) const;
            wasAborted_function_type wasAborted_function_value( &::SireCluster::WorkPacket::wasAborted );
            
            WorkPacket_exposer.def( 
                "wasAborted"
                , wasAborted_function_value );
        
        }
        { //::SireCluster::WorkPacket::what
        
            typedef char const * ( ::SireCluster::WorkPacket::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCluster::WorkPacket::what );
            
            WorkPacket_exposer.def( 
                "what"
                , what_function_value );
        
        }
        WorkPacket_exposer.staticmethod( "typeName" );
        WorkPacket_exposer.staticmethod( "unpack" );
        WorkPacket_exposer.def( "__copy__", &__copy__);
        WorkPacket_exposer.def( "__deepcopy__", &__copy__);
        WorkPacket_exposer.def( "clone", &__copy__);
        WorkPacket_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCluster::WorkPacket >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        WorkPacket_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCluster::WorkPacket >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        WorkPacket_exposer.def( "__str__", &pvt_get_name);
        WorkPacket_exposer.def( "__repr__", &pvt_get_name);
    }

}
