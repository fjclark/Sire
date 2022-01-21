// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FEP.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "fep.h"

#include "tostring.h"

#include "fep.h"

SireAnalysis::FEP __copy__(const SireAnalysis::FEP &other){ return SireAnalysis::FEP(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_FEP_class(){

    { //::SireAnalysis::FEP
        typedef bp::class_< SireAnalysis::FEP, bp::bases< SireBase::Property > > FEP_exposer_t;
        FEP_exposer_t FEP_exposer = FEP_exposer_t( "FEP", "This class is used to analyse the free energies that are\ncalculated during a free energy perturbation (FEP) simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope FEP_scope( FEP_exposer );
        FEP_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("deltas") ), "Construct to use the passed set of windows, with the free energy deltas from\neach window to the window above") );
        FEP_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::FreeEnergyAverage > const &, QMap< double, SireMaths::FreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_deltas"), bp::arg("backwards_deltas") ), "Construct to use the passed windows, with the free energy deltas from\neach window to the window above in forwards_deltas and from the window\nbelow to each window in backwards_deltas") );
        FEP_exposer.def( bp::init< SireAnalysis::FEPDeltas const & >(( bp::arg("deltas") ), "Construct to use the passed FEP deltas") );
        FEP_exposer.def( bp::init< SireAnalysis::FEP const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::FEP::add
        
            typedef void ( ::SireAnalysis::FEP::*add_function_type)( ::QList< double > const &,::QMap< double, SireMaths::FreeEnergyAverage > const & ) ;
            add_function_type add_function_value( &::SireAnalysis::FEP::add );
            
            FEP_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("windows"), bp::arg("deltas") )
                , "Add the data for the next iteration, which contains the deltas for the passed windows,\nwith the free energy being for each window to the next window" );
        
        }
        { //::SireAnalysis::FEP::add
        
            typedef void ( ::SireAnalysis::FEP::*add_function_type)( ::QList< double > const &,::QMap< double, SireMaths::FreeEnergyAverage > const &,::QMap< double, SireMaths::FreeEnergyAverage > const & ) ;
            add_function_type add_function_value( &::SireAnalysis::FEP::add );
            
            FEP_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("windows"), bp::arg("forwards_deltas"), bp::arg("backwards_deltas") )
                , "Add the data for the next iteration, which contains the deltas for the passed windows,\nwith forwards_deltas containing the free energy from each window to the next window,\nand backwards_deltas containing the free energy from the previous window to each window" );
        
        }
        { //::SireAnalysis::FEP::add
        
            typedef void ( ::SireAnalysis::FEP::*add_function_type)( ::SireAnalysis::FEPDeltas const & ) ;
            add_function_type add_function_value( &::SireAnalysis::FEP::add );
            
            FEP_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("deltas") )
                , "Add the data for the next iteration" );
        
        }
        { //::SireAnalysis::FEP::at
        
            typedef ::SireAnalysis::FEPDeltas ( ::SireAnalysis::FEP::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireAnalysis::FEP::at );
            
            FEP_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , "Return the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::FEP::clear
        
            typedef void ( ::SireAnalysis::FEP::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireAnalysis::FEP::clear );
            
            FEP_exposer.def( 
                "clear"
                , clear_function_value
                , "Remove all values from the histogram" );
        
        }
        { //::SireAnalysis::FEP::count
        
            typedef int ( ::SireAnalysis::FEP::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireAnalysis::FEP::count );
            
            FEP_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::FEP::deltas
        
            typedef ::QList< SireAnalysis::FEPDeltas > ( ::SireAnalysis::FEP::*deltas_function_type)(  ) const;
            deltas_function_type deltas_function_value( &::SireAnalysis::FEP::deltas );
            
            FEP_exposer.def( 
                "deltas"
                , deltas_function_value
                , "Return the deltas for all iterations" );
        
        }
        { //::SireAnalysis::FEP::lambdaValues
        
            typedef ::QList< double > ( ::SireAnalysis::FEP::*lambdaValues_function_type)(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::SireAnalysis::FEP::lambdaValues );
            
            FEP_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value
                , "Return the values of all windows" );
        
        }
        { //::SireAnalysis::FEP::merge
        
            typedef ::SireAnalysis::FEPDeltas ( ::SireAnalysis::FEP::*merge_function_type)( int,int ) const;
            merge_function_type merge_function_value( &::SireAnalysis::FEP::merge );
            
            FEP_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "Merge the deltas for iterations start->end" );
        
        }
        { //::SireAnalysis::FEP::merge
        
            typedef ::SireAnalysis::FEPDeltas ( ::SireAnalysis::FEP::*merge_function_type)( ::QList< int > ) const;
            merge_function_type merge_function_value( &::SireAnalysis::FEP::merge );
            
            FEP_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("indicies") )
                , "Merge the deltas at the passed indicies" );
        
        }
        { //::SireAnalysis::FEP::nIterations
        
            typedef int ( ::SireAnalysis::FEP::*nIterations_function_type)(  ) const;
            nIterations_function_type nIterations_function_value( &::SireAnalysis::FEP::nIterations );
            
            FEP_exposer.def( 
                "nIterations"
                , nIterations_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::FEP::nLambdaValues
        
            typedef int ( ::SireAnalysis::FEP::*nLambdaValues_function_type)(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::SireAnalysis::FEP::nLambdaValues );
            
            FEP_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value
                , "Return the number of lambda values (windows)" );
        
        }
        { //::SireAnalysis::FEP::nSamples
        
            typedef ::qint64 ( ::SireAnalysis::FEP::*nSamples_function_type)(  ) const;
            nSamples_function_type nSamples_function_value( &::SireAnalysis::FEP::nSamples );
            
            FEP_exposer.def( 
                "nSamples"
                , nSamples_function_value
                , "Return the total number of samples in the simulation" );
        
        }
        { //::SireAnalysis::FEP::nWindows
        
            typedef int ( ::SireAnalysis::FEP::*nWindows_function_type)(  ) const;
            nWindows_function_type nWindows_function_value( &::SireAnalysis::FEP::nWindows );
            
            FEP_exposer.def( 
                "nWindows"
                , nWindows_function_value
                , "Return the number of windows" );
        
        }
        FEP_exposer.def( bp::self != bp::self );
        { //::SireAnalysis::FEP::operator=
        
            typedef ::SireAnalysis::FEP & ( ::SireAnalysis::FEP::*assign_function_type)( ::SireAnalysis::FEP const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::FEP::operator= );
            
            FEP_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FEP_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::FEP::operator[]
        
            typedef ::SireAnalysis::FEPDeltas ( ::SireAnalysis::FEP::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireAnalysis::FEP::operator[] );
            
            FEP_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireAnalysis::FEP::removeAt
        
            typedef void ( ::SireAnalysis::FEP::*removeAt_function_type)( int ) ;
            removeAt_function_type removeAt_function_value( &::SireAnalysis::FEP::removeAt );
            
            FEP_exposer.def( 
                "removeAt"
                , removeAt_function_value
                , ( bp::arg("i") )
                , "Remove the data for iteration i" );
        
        }
        { //::SireAnalysis::FEP::removeRange
        
            typedef void ( ::SireAnalysis::FEP::*removeRange_function_type)( int,int ) ;
            removeRange_function_type removeRange_function_value( &::SireAnalysis::FEP::removeRange );
            
            FEP_exposer.def( 
                "removeRange"
                , removeRange_function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "Remove every iteration from start to end (inclusively)" );
        
        }
        { //::SireAnalysis::FEP::rollingAverage
        
            typedef ::QList< SireAnalysis::FEPDeltas > ( ::SireAnalysis::FEP::*rollingAverage_function_type)( int ) const;
            rollingAverage_function_type rollingAverage_function_value( &::SireAnalysis::FEP::rollingAverage );
            
            FEP_exposer.def( 
                "rollingAverage"
                , rollingAverage_function_value
                , ( bp::arg("niterations") )
                , "Return a list of Gradients that represents the rolling average over niterations\niterations over this TI data set. If this data set contains 100 iterations, and\nwe calculate the rolling average over 50 iterations, then the returned Gradients\nwill be the average from 1-50, then 2-51, 3-52.....51-100" );
        
        }
        { //::SireAnalysis::FEP::set
        
            typedef void ( ::SireAnalysis::FEP::*set_function_type)( int,::QList< double > const &,::QMap< double, SireMaths::FreeEnergyAverage > const & ) ;
            set_function_type set_function_value( &::SireAnalysis::FEP::set );
            
            FEP_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("windows"), bp::arg("deltas") )
                , "Set the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::FEP::set
        
            typedef void ( ::SireAnalysis::FEP::*set_function_type)( int,::QList< double > const &,::QMap< double, SireMaths::FreeEnergyAverage > const &,::QMap< double, SireMaths::FreeEnergyAverage > const & ) ;
            set_function_type set_function_value( &::SireAnalysis::FEP::set );
            
            FEP_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("windows"), bp::arg("forwards_deltas"), bp::arg("backwards_deltas") )
                , "Set the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::FEP::set
        
            typedef void ( ::SireAnalysis::FEP::*set_function_type)( int,::SireAnalysis::FEPDeltas const & ) ;
            set_function_type set_function_value( &::SireAnalysis::FEP::set );
            
            FEP_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("deltas") )
                , "Set the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::FEP::size
        
            typedef int ( ::SireAnalysis::FEP::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireAnalysis::FEP::size );
            
            FEP_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::FEP::toString
        
            typedef ::QString ( ::SireAnalysis::FEP::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::FEP::toString );
            
            FEP_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireAnalysis::FEP::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::FEP::typeName );
            
            FEP_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireAnalysis::FEP::what
        
            typedef char const * ( ::SireAnalysis::FEP::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::FEP::what );
            
            FEP_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        { //::SireAnalysis::FEP::windows
        
            typedef ::QList< double > ( ::SireAnalysis::FEP::*windows_function_type)(  ) const;
            windows_function_type windows_function_value( &::SireAnalysis::FEP::windows );
            
            FEP_exposer.def( 
                "windows"
                , windows_function_value
                , "Return the value of all windows" );
        
        }
        FEP_exposer.staticmethod( "typeName" );
        FEP_exposer.def( "__copy__", &__copy__);
        FEP_exposer.def( "__deepcopy__", &__copy__);
        FEP_exposer.def( "clone", &__copy__);
        FEP_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::FEP >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEP_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::FEP >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FEP_exposer.def( "__setstate__", &__setstate__base64< ::SireAnalysis::FEP > );
        FEP_exposer.def( "__getstate__", &__getstate__base64< ::SireAnalysis::FEP > );
        FEP_exposer.def( "__str__", &__str__< ::SireAnalysis::FEP > );
        FEP_exposer.def( "__repr__", &__str__< ::SireAnalysis::FEP > );
        FEP_exposer.def( "__len__", &__len_size< ::SireAnalysis::FEP > );
    }

}
