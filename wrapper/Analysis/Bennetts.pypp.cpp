// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Bennetts.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/maths.h"

#include "SireStream/registeralternativename.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/temperature.h"

#include "SireUnits/units.h"

#include "bennetts.h"

#include "tostring.h"

#include "bennetts.h"

SireAnalysis::Bennetts __copy__(const SireAnalysis::Bennetts &other){ return SireAnalysis::Bennetts(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_Bennetts_class(){

    { //::SireAnalysis::Bennetts
        typedef bp::class_< SireAnalysis::Bennetts, bp::bases< SireBase::Property > > Bennetts_exposer_t;
        Bennetts_exposer_t Bennetts_exposer = Bennetts_exposer_t( "Bennetts", "This class is used to analyse the free energies that are\ncalculated during a Bennetts Acceptance Ratio simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope Bennetts_scope( Bennetts_exposer );
        Bennetts_exposer.def( bp::init< QList< double > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const &, QMap< double, SireMaths::BennettsFreeEnergyAverage > const & >(( bp::arg("windows"), bp::arg("forwards_ratios"), bp::arg("backwards_ratios") ), "Construct to use the passed windows, with the free energy ratios from\neach window to the window above in forwards_ratios and from the window\nbelow to each window in backwards_ratios") );
        Bennetts_exposer.def( bp::init< SireAnalysis::BennettsRatios const & >(( bp::arg("ratios") ), "Construct to use the passed Bennetts ratios") );
        Bennetts_exposer.def( bp::init< SireAnalysis::Bennetts const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireAnalysis::Bennetts::add
        
            typedef void ( ::SireAnalysis::Bennetts::*add_function_type)( ::QList< double > const &,::QMap< double, SireMaths::BennettsFreeEnergyAverage > const &,::QMap< double, SireMaths::BennettsFreeEnergyAverage > const & ) ;
            add_function_type add_function_value( &::SireAnalysis::Bennetts::add );
            
            Bennetts_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("windows"), bp::arg("forwards_ratios"), bp::arg("backwards_ratios") )
                , "Add the data for the next iteration, which contains the ratios for the passed windows,\nwith forwards_ratios containing the free energy from each window to the next window,\nand backwards_ratios containing the free energy from the previous window to each window" );
        
        }
        { //::SireAnalysis::Bennetts::add
        
            typedef void ( ::SireAnalysis::Bennetts::*add_function_type)( ::SireAnalysis::BennettsRatios const & ) ;
            add_function_type add_function_value( &::SireAnalysis::Bennetts::add );
            
            Bennetts_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("ratios") )
                , "Add the data for the next iteration" );
        
        }
        { //::SireAnalysis::Bennetts::at
        
            typedef ::SireAnalysis::BennettsRatios ( ::SireAnalysis::Bennetts::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireAnalysis::Bennetts::at );
            
            Bennetts_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , "Return the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::Bennetts::clear
        
            typedef void ( ::SireAnalysis::Bennetts::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireAnalysis::Bennetts::clear );
            
            Bennetts_exposer.def( 
                "clear"
                , clear_function_value
                , "Remove all values from the histogram" );
        
        }
        { //::SireAnalysis::Bennetts::count
        
            typedef int ( ::SireAnalysis::Bennetts::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireAnalysis::Bennetts::count );
            
            Bennetts_exposer.def( 
                "count"
                , count_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::Bennetts::lambdaValues
        
            typedef ::QList< double > ( ::SireAnalysis::Bennetts::*lambdaValues_function_type)(  ) const;
            lambdaValues_function_type lambdaValues_function_value( &::SireAnalysis::Bennetts::lambdaValues );
            
            Bennetts_exposer.def( 
                "lambdaValues"
                , lambdaValues_function_value
                , "Return the values of all windows" );
        
        }
        { //::SireAnalysis::Bennetts::merge
        
            typedef ::SireAnalysis::BennettsRatios ( ::SireAnalysis::Bennetts::*merge_function_type)( int,int ) const;
            merge_function_type merge_function_value( &::SireAnalysis::Bennetts::merge );
            
            Bennetts_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "Merge the deltas for iterations start->end" );
        
        }
        { //::SireAnalysis::Bennetts::merge
        
            typedef ::SireAnalysis::BennettsRatios ( ::SireAnalysis::Bennetts::*merge_function_type)( ::QList< int > ) const;
            merge_function_type merge_function_value( &::SireAnalysis::Bennetts::merge );
            
            Bennetts_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("indicies") )
                , "Merge the deltas at the passed indicies" );
        
        }
        { //::SireAnalysis::Bennetts::nIterations
        
            typedef int ( ::SireAnalysis::Bennetts::*nIterations_function_type)(  ) const;
            nIterations_function_type nIterations_function_value( &::SireAnalysis::Bennetts::nIterations );
            
            Bennetts_exposer.def( 
                "nIterations"
                , nIterations_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::Bennetts::nLambdaValues
        
            typedef int ( ::SireAnalysis::Bennetts::*nLambdaValues_function_type)(  ) const;
            nLambdaValues_function_type nLambdaValues_function_value( &::SireAnalysis::Bennetts::nLambdaValues );
            
            Bennetts_exposer.def( 
                "nLambdaValues"
                , nLambdaValues_function_value
                , "Return the number of lambda values (windows)" );
        
        }
        { //::SireAnalysis::Bennetts::nSamples
        
            typedef ::qint64 ( ::SireAnalysis::Bennetts::*nSamples_function_type)(  ) const;
            nSamples_function_type nSamples_function_value( &::SireAnalysis::Bennetts::nSamples );
            
            Bennetts_exposer.def( 
                "nSamples"
                , nSamples_function_value
                , "Return the total number of samples in the simulation" );
        
        }
        { //::SireAnalysis::Bennetts::nWindows
        
            typedef int ( ::SireAnalysis::Bennetts::*nWindows_function_type)(  ) const;
            nWindows_function_type nWindows_function_value( &::SireAnalysis::Bennetts::nWindows );
            
            Bennetts_exposer.def( 
                "nWindows"
                , nWindows_function_value
                , "Return the number of windows" );
        
        }
        Bennetts_exposer.def( bp::self != bp::self );
        { //::SireAnalysis::Bennetts::operator=
        
            typedef ::SireAnalysis::Bennetts & ( ::SireAnalysis::Bennetts::*assign_function_type)( ::SireAnalysis::Bennetts const & ) ;
            assign_function_type assign_function_value( &::SireAnalysis::Bennetts::operator= );
            
            Bennetts_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Bennetts_exposer.def( bp::self == bp::self );
        { //::SireAnalysis::Bennetts::operator[]
        
            typedef ::SireAnalysis::BennettsRatios ( ::SireAnalysis::Bennetts::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireAnalysis::Bennetts::operator[] );
            
            Bennetts_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireAnalysis::Bennetts::ratios
        
            typedef ::QList< SireAnalysis::BennettsRatios > ( ::SireAnalysis::Bennetts::*ratios_function_type)(  ) const;
            ratios_function_type ratios_function_value( &::SireAnalysis::Bennetts::ratios );
            
            Bennetts_exposer.def( 
                "ratios"
                , ratios_function_value
                , "Return the deltas for all iterations" );
        
        }
        { //::SireAnalysis::Bennetts::removeAt
        
            typedef void ( ::SireAnalysis::Bennetts::*removeAt_function_type)( int ) ;
            removeAt_function_type removeAt_function_value( &::SireAnalysis::Bennetts::removeAt );
            
            Bennetts_exposer.def( 
                "removeAt"
                , removeAt_function_value
                , ( bp::arg("i") )
                , "Remove the data for iteration i" );
        
        }
        { //::SireAnalysis::Bennetts::removeRange
        
            typedef void ( ::SireAnalysis::Bennetts::*removeRange_function_type)( int,int ) ;
            removeRange_function_type removeRange_function_value( &::SireAnalysis::Bennetts::removeRange );
            
            Bennetts_exposer.def( 
                "removeRange"
                , removeRange_function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "Remove every iteration from start to end (inclusively)" );
        
        }
        { //::SireAnalysis::Bennetts::rollingAverage
        
            typedef ::QList< SireAnalysis::BennettsRatios > ( ::SireAnalysis::Bennetts::*rollingAverage_function_type)( int ) const;
            rollingAverage_function_type rollingAverage_function_value( &::SireAnalysis::Bennetts::rollingAverage );
            
            Bennetts_exposer.def( 
                "rollingAverage"
                , rollingAverage_function_value
                , ( bp::arg("niterations") )
                , "Return a list of Gradients that represents the rolling average over niterations\niterations over this TI data set. If this data set contains 100 iterations, and\nwe calculate the rolling average over 50 iterations, then the returned Gradients\nwill be the average from 1-50, then 2-51, 3-52.....51-100" );
        
        }
        { //::SireAnalysis::Bennetts::set
        
            typedef void ( ::SireAnalysis::Bennetts::*set_function_type)( int,::QList< double > const &,::QMap< double, SireMaths::BennettsFreeEnergyAverage > const &,::QMap< double, SireMaths::BennettsFreeEnergyAverage > const & ) ;
            set_function_type set_function_value( &::SireAnalysis::Bennetts::set );
            
            Bennetts_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("windows"), bp::arg("forwards_ratios"), bp::arg("backwards_ratios") )
                , "Set the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::Bennetts::set
        
            typedef void ( ::SireAnalysis::Bennetts::*set_function_type)( int,::SireAnalysis::BennettsRatios const & ) ;
            set_function_type set_function_value( &::SireAnalysis::Bennetts::set );
            
            Bennetts_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("ratios") )
                , "Set the deltas for the ith iteration" );
        
        }
        { //::SireAnalysis::Bennetts::size
        
            typedef int ( ::SireAnalysis::Bennetts::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireAnalysis::Bennetts::size );
            
            Bennetts_exposer.def( 
                "size"
                , size_function_value
                , "Return the number of iterations" );
        
        }
        { //::SireAnalysis::Bennetts::toString
        
            typedef ::QString ( ::SireAnalysis::Bennetts::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireAnalysis::Bennetts::toString );
            
            Bennetts_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireAnalysis::Bennetts::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireAnalysis::Bennetts::typeName );
            
            Bennetts_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireAnalysis::Bennetts::what
        
            typedef char const * ( ::SireAnalysis::Bennetts::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireAnalysis::Bennetts::what );
            
            Bennetts_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        { //::SireAnalysis::Bennetts::windows
        
            typedef ::QList< double > ( ::SireAnalysis::Bennetts::*windows_function_type)(  ) const;
            windows_function_type windows_function_value( &::SireAnalysis::Bennetts::windows );
            
            Bennetts_exposer.def( 
                "windows"
                , windows_function_value
                , "Return the value of all windows" );
        
        }
        Bennetts_exposer.staticmethod( "typeName" );
        Bennetts_exposer.def( "__copy__", &__copy__);
        Bennetts_exposer.def( "__deepcopy__", &__copy__);
        Bennetts_exposer.def( "clone", &__copy__);
        Bennetts_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireAnalysis::Bennetts >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bennetts_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireAnalysis::Bennetts >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bennetts_exposer.def( "__setstate__", &__setstate__base64< ::SireAnalysis::Bennetts > );
        Bennetts_exposer.def( "__getstate__", &__getstate__base64< ::SireAnalysis::Bennetts > );
        Bennetts_exposer.def( "__str__", &__str__< ::SireAnalysis::Bennetts > );
        Bennetts_exposer.def( "__repr__", &__str__< ::SireAnalysis::Bennetts > );
        Bennetts_exposer.def( "__len__", &__len_size< ::SireAnalysis::Bennetts > );
    }

}
