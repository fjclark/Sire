// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "QMProgram.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/molecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "latticecharges.h"

#include "qmprogram.h"

#include <QMutex>

#include "qmprogram.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_QMProgram_class(){

    { //::Squire::QMProgram
        typedef bp::class_< Squire::QMProgram, bp::bases< SireBase::Property >, boost::noncopyable > QMProgram_exposer_t;
        QMProgram_exposer_t QMProgram_exposer = QMProgram_exposer_t( "QMProgram", "This is the base class of all QM programs. These are wrappers that\nprovide the functionality to calculate QM energies and forces\nby calling separate QM programs\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope QMProgram_scope( QMProgram_exposer );
        { //::Squire::QMProgram::calculateCharges
        
            typedef ::SireMol::AtomCharges ( ::Squire::QMProgram::*calculateCharges_function_type)( ::SireMol::Molecule const &,::SireBase::PropertyMap const & ) const;
            calculateCharges_function_type calculateCharges_function_value( &::Squire::QMProgram::calculateCharges );
            
            QMProgram_exposer.def( 
                "calculateCharges"
                , calculateCharges_function_value
                , ( bp::arg("molecule"), bp::arg("map") )
                , "Calculate the charges on the molecule molecule using the properties\nspecified in the passed property map" );
        
        }
        { //::Squire::QMProgram::calculateCharges
        
            typedef ::SireMol::AtomCharges ( ::Squire::QMProgram::*calculateCharges_function_type)( ::SireMol::Molecule const & ) const;
            calculateCharges_function_type calculateCharges_function_value( &::Squire::QMProgram::calculateCharges );
            
            QMProgram_exposer.def( 
                "calculateCharges"
                , calculateCharges_function_value
                , ( bp::arg("molecule") )
                , "Calculate the charges on the molecule molecule using the default\nproperty locations" );
        
        }
        { //::Squire::QMProgram::chargeCommandFile
        
            typedef ::QString ( ::Squire::QMProgram::*chargeCommandFile_function_type)( ::SireMol::Molecule const & ) const;
            chargeCommandFile_function_type chargeCommandFile_function_value( &::Squire::QMProgram::chargeCommandFile );
            
            QMProgram_exposer.def( 
                "chargeCommandFile"
                , chargeCommandFile_function_value
                , ( bp::arg("molecule") )
                , "Return the command file that would be used to calculate the atomic\npartial charges of the passed molecule" );
        
        }
        { //::Squire::QMProgram::chargeCommandFile
        
            typedef ::QString ( ::Squire::QMProgram::*chargeCommandFile_function_type)( ::SireMol::Molecule const &,::SireBase::PropertyMap const & ) const;
            chargeCommandFile_function_type chargeCommandFile_function_value( &::Squire::QMProgram::chargeCommandFile );
            
            QMProgram_exposer.def( 
                "chargeCommandFile"
                , chargeCommandFile_function_value
                , ( bp::arg("molecule"), bp::arg("map") )
                , "Return the command file that would be used to calculate the atomic\npartial charges of the passed molecule" );
        
        }
        { //::Squire::QMProgram::null
        
            typedef ::Squire::NullQM const & ( *null_function_type )(  );
            null_function_type null_function_value( &::Squire::QMProgram::null );
            
            QMProgram_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::Squire::QMProgram::numberOfMMAtomsLimit
        
            typedef int ( ::Squire::QMProgram::*numberOfMMAtomsLimit_function_type)(  ) const;
            numberOfMMAtomsLimit_function_type numberOfMMAtomsLimit_function_value( &::Squire::QMProgram::numberOfMMAtomsLimit );
            
            QMProgram_exposer.def( 
                "numberOfMMAtomsLimit"
                , numberOfMMAtomsLimit_function_value
                , "Return the maximum number of MM atoms supported by this QM program. This\nreturns -1 if there is no limit" );
        
        }
        { //::Squire::QMProgram::numberOfMMAtomsLimit
        
            typedef int ( ::Squire::QMProgram::*numberOfMMAtomsLimit_function_type)( int ) const;
            numberOfMMAtomsLimit_function_type numberOfMMAtomsLimit_function_value( &::Squire::QMProgram::numberOfMMAtomsLimit );
            
            QMProgram_exposer.def( 
                "numberOfMMAtomsLimit"
                , numberOfMMAtomsLimit_function_value
                , ( bp::arg("num_qm_atoms") )
                , "Return the maximum number of QM atoms supported by this QM program,\ngiven the supplied number of MM atoms. This returns -1 if there is\nno limit" );
        
        }
        { //::Squire::QMProgram::supportsGaussianCharges
        
            typedef bool ( ::Squire::QMProgram::*supportsGaussianCharges_function_type)(  ) const;
            supportsGaussianCharges_function_type supportsGaussianCharges_function_value( &::Squire::QMProgram::supportsGaussianCharges );
            
            QMProgram_exposer.def( 
                "supportsGaussianCharges"
                , supportsGaussianCharges_function_value
                , "Return whether or not this QM program supports the use\nof gaussian lattice charges (which can polarise the QM wavefunction)" );
        
        }
        { //::Squire::QMProgram::supportsLatticeCharges
        
            typedef bool ( ::Squire::QMProgram::*supportsLatticeCharges_function_type)(  ) const;
            supportsLatticeCharges_function_type supportsLatticeCharges_function_value( &::Squire::QMProgram::supportsLatticeCharges );
            
            QMProgram_exposer.def( 
                "supportsLatticeCharges"
                , supportsLatticeCharges_function_value
                , "Return whether or not this QM program supports the use\nof point lattice charges (which can polarise the QM wavefunction)" );
        
        }
        { //::Squire::QMProgram::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::QMProgram::typeName );
            
            QMProgram_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        QMProgram_exposer.staticmethod( "null" );
        QMProgram_exposer.staticmethod( "typeName" );
        QMProgram_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::QMProgram >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMProgram_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::QMProgram >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        QMProgram_exposer.def( "__getstate_manages_dict__", true);
        QMProgram_exposer.def( "__safe_for_unpickling__", true);
        QMProgram_exposer.def( "__setstate__", &__setstate__base64< ::Squire::QMProgram > );
        QMProgram_exposer.def( "__getstate__", &__getstate__base64< ::Squire::QMProgram > );
        QMProgram_exposer.def( "__str__", &__str__< ::Squire::QMProgram > );
        QMProgram_exposer.def( "__repr__", &__str__< ::Squire::QMProgram > );
    }

}
