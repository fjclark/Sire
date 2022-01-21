// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SameSupraMoves.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "supramoves.h"

#include "suprasystem.h"

#include "supramoves.h"

SireMove::SameSupraMoves __copy__(const SireMove::SameSupraMoves &other){ return SireMove::SameSupraMoves(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_SameSupraMoves_class(){

    { //::SireMove::SameSupraMoves
        typedef bp::class_< SireMove::SameSupraMoves, bp::bases< SireMove::SupraMoves, SireBase::Property > > SameSupraMoves_exposer_t;
        SameSupraMoves_exposer_t SameSupraMoves_exposer = SameSupraMoves_exposer_t( "SameSupraMoves", "This SupraMoves object is used to apply the same SupraMove\nrepeatedly to a SupraSystem\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope SameSupraMoves_scope( SameSupraMoves_exposer );
        SameSupraMoves_exposer.def( bp::init< SireMove::SupraMove const & >(( bp::arg("move") ), "Construct to run the move move repeatedly") );
        SameSupraMoves_exposer.def( bp::init< SireMove::SameSupraMoves const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::SameSupraMoves::clearStatistics
        
            typedef void ( ::SireMove::SameSupraMoves::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireMove::SameSupraMoves::clearStatistics );
            
            SameSupraMoves_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , "Clear all move statistics" );
        
        }
        { //::SireMove::SameSupraMoves::move
        
            typedef void ( ::SireMove::SameSupraMoves::*move_function_type)( ::SireMove::SupraSystem &,int,bool ) ;
            move_function_type move_function_value( &::SireMove::SameSupraMoves::move );
            
            SameSupraMoves_exposer.def( 
                "move"
                , move_function_value
                , ( bp::arg("system"), bp::arg("nmoves"), bp::arg("record_stats")=(bool)(true) )
                , "Perform the moves nmoves times" );
        
        }
        { //::SireMove::SameSupraMoves::moves
        
            typedef ::QList< SireBase::PropPtr< SireMove::SupraMove > > ( ::SireMove::SameSupraMoves::*moves_function_type)(  ) const;
            moves_function_type moves_function_value( &::SireMove::SameSupraMoves::moves );
            
            SameSupraMoves_exposer.def( 
                "moves"
                , moves_function_value
                , "Return a list of all of the moves" );
        
        }
        { //::SireMove::SameSupraMoves::nMoves
        
            typedef int ( ::SireMove::SameSupraMoves::*nMoves_function_type)(  ) const;
            nMoves_function_type nMoves_function_value( &::SireMove::SameSupraMoves::nMoves );
            
            SameSupraMoves_exposer.def( 
                "nMoves"
                , nMoves_function_value
                , "Return the total number of moves that have been performed" );
        
        }
        SameSupraMoves_exposer.def( bp::self != bp::self );
        { //::SireMove::SameSupraMoves::operator=
        
            typedef ::SireMove::SameSupraMoves & ( ::SireMove::SameSupraMoves::*assign_function_type)( ::SireMove::SameSupraMoves const & ) ;
            assign_function_type assign_function_value( &::SireMove::SameSupraMoves::operator= );
            
            SameSupraMoves_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SameSupraMoves_exposer.def( bp::self == bp::self );
        { //::SireMove::SameSupraMoves::toString
        
            typedef ::QString ( ::SireMove::SameSupraMoves::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMove::SameSupraMoves::toString );
            
            SameSupraMoves_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this moves set" );
        
        }
        { //::SireMove::SameSupraMoves::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::SameSupraMoves::typeName );
            
            SameSupraMoves_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        SameSupraMoves_exposer.staticmethod( "typeName" );
        SameSupraMoves_exposer.def( "__copy__", &__copy__);
        SameSupraMoves_exposer.def( "__deepcopy__", &__copy__);
        SameSupraMoves_exposer.def( "clone", &__copy__);
        SameSupraMoves_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::SameSupraMoves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SameSupraMoves_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::SameSupraMoves >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SameSupraMoves_exposer.def_pickle(sire_pickle_suite< ::SireMove::SameSupraMoves >());
        SameSupraMoves_exposer.def( "__str__", &__str__< ::SireMove::SameSupraMoves > );
        SameSupraMoves_exposer.def( "__repr__", &__str__< ::SireMove::SameSupraMoves > );
        SameSupraMoves_exposer.def( "__len__", &__len_size< ::SireMove::SameSupraMoves > );
    }

}
