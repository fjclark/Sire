// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FFIdx.pypp.hpp"

namespace bp = boost::python;

#include "ffidx.h"

#include "ffidx.h"

#include "forcefields.h"

#include "ffidx.h"

SireFF::FFIdx __copy__(const SireFF::FFIdx &other){ return SireFF::FFIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_FFIdx_class(){

    { //::SireFF::FFIdx
        typedef bp::class_< SireFF::FFIdx, bp::bases< SireFF::FFID, SireID::ID, SireID::IndexBase > > FFIdx_exposer_t;
        FFIdx_exposer_t FFIdx_exposer = FFIdx_exposer_t( "FFIdx", "This is an ID object that is used to index forcefields (e.g. index\nin a list or array).\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope FFIdx_scope( FFIdx_exposer );
        FFIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") ), "") );
        FFIdx_exposer.def( bp::init< SireFF::FFIdx const & >(( bp::arg("other") ), "") );
        { //::SireFF::FFIdx::hash
        
            typedef ::uint ( ::SireFF::FFIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireFF::FFIdx::hash );
            
            FFIdx_exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireFF::FFIdx::isNull
        
            typedef bool ( ::SireFF::FFIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireFF::FFIdx::isNull );
            
            FFIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireFF::FFIdx::map
        
            typedef ::QList< SireFF::FFIdx > ( ::SireFF::FFIdx::*map_function_type)( ::SireFF::ForceFields const & ) const;
            map_function_type map_function_value( &::SireFF::FFIdx::map );
            
            FFIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("ffields") )
                , "Short cut function to map this index to the index of the\nmatching forcefield in the passed ForceFields object\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireFF::FFIdx::null
        
            typedef ::SireFF::FFIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireFF::FFIdx::null );
            
            FFIdx_exposer.def( 
                "null"
                , null_function_value
                , "" );
        
        }
        { //::SireFF::FFIdx::operator=
        
            typedef ::SireFF::FFIdx & ( ::SireFF::FFIdx::*assign_function_type)( ::SireFF::FFIdx const & ) ;
            assign_function_type assign_function_value( &::SireFF::FFIdx::operator= );
            
            FFIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::FFIdx::toString
        
            typedef ::QString ( ::SireFF::FFIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireFF::FFIdx::toString );
            
            FFIdx_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireFF::FFIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::FFIdx::typeName );
            
            FFIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::FFIdx::what
        
            typedef char const * ( ::SireFF::FFIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::FFIdx::what );
            
            FFIdx_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        FFIdx_exposer.staticmethod( "null" );
        FFIdx_exposer.staticmethod( "typeName" );
        FFIdx_exposer.def( "__copy__", &__copy__);
        FFIdx_exposer.def( "__deepcopy__", &__copy__);
        FFIdx_exposer.def( "clone", &__copy__);
        FFIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::FFIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::FFIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFIdx_exposer.def( "__getstate_manages_dict__", true);
        FFIdx_exposer.def( "__safe_for_unpickling__", true);
        FFIdx_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::FFIdx > );
        FFIdx_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::FFIdx > );
        FFIdx_exposer.def( "__str__", &__str__< ::SireFF::FFIdx > );
        FFIdx_exposer.def( "__repr__", &__str__< ::SireFF::FFIdx > );
        FFIdx_exposer.def( "__hash__", &::SireFF::FFIdx::hash );
    }

}
