// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "PDB2.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/pdb2.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleditor.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "pdb2.h"

#include <QFile>

#include <QtMath>

#include "pdb2.h"

SireIO::PDB2 __copy__(const SireIO::PDB2 &other){ return SireIO::PDB2(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_PDB2_class(){

    { //::SireIO::PDB2
        typedef bp::class_< SireIO::PDB2, bp::bases< SireIO::MoleculeParser, SireBase::Property > > PDB2_exposer_t;
        PDB2_exposer_t PDB2_exposer = PDB2_exposer_t( "PDB2", "This class holds a parser for reading and writing\nProtein Data Bank (PDB) files\n\nAuthor: Lester Hedges\n", bp::init< >("Constructor") );
        bp::scope PDB2_scope( PDB2_exposer );
        PDB2_exposer.def( bp::init< QString const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("filename"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the file called filename. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        PDB2_exposer.def( bp::init< QStringList const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("lines"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to read in the data from the passed text lines. The\npassed property map can be used to pass extra parameters to control\nthe parsing") );
        PDB2_exposer.def( bp::init< SireSystem::System const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("system"), bp::arg("map")=SireBase::PropertyMap() ), "Construct this parser by extracting all necessary information from the\npassed SireSystem::System, looking for the properties that are specified\nin the passed property map") );
        PDB2_exposer.def( bp::init< SireIO::PDB2 const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::PDB2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDB2::*construct_function_type)( ::QString const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDB2::construct );
            
            PDB2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("filename"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\nfile using the passed properties" );
        
        }
        { //::SireIO::PDB2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDB2::*construct_function_type)( ::QStringList const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDB2::construct );
            
            PDB2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("lines"), bp::arg("map") )
                , "Return the parser that has been constructed by reading in the passed\ntext lines using the passed properties" );
        
        }
        { //::SireIO::PDB2::construct
        
            typedef ::SireIO::MoleculeParserPtr ( ::SireIO::PDB2::*construct_function_type)( ::SireSystem::System const &,::SireBase::PropertyMap const & ) const;
            construct_function_type construct_function_value( &::SireIO::PDB2::construct );
            
            PDB2_exposer.def( 
                "construct"
                , construct_function_value
                , ( bp::arg("system"), bp::arg("map") )
                , "Return the parser that has been constructed by extract all necessary\ndata from the passed SireSystem::System using the specified properties" );
        
        }
        { //::SireIO::PDB2::formatDescription
        
            typedef ::QString ( ::SireIO::PDB2::*formatDescription_function_type)(  ) const;
            formatDescription_function_type formatDescription_function_value( &::SireIO::PDB2::formatDescription );
            
            PDB2_exposer.def( 
                "formatDescription"
                , formatDescription_function_value
                , "Return a description of the file format" );
        
        }
        { //::SireIO::PDB2::formatName
        
            typedef ::QString ( ::SireIO::PDB2::*formatName_function_type)(  ) const;
            formatName_function_type formatName_function_value( &::SireIO::PDB2::formatName );
            
            PDB2_exposer.def( 
                "formatName"
                , formatName_function_value
                , "Return the format name that is used to identify this file format within Sire" );
        
        }
        { //::SireIO::PDB2::formatSuffix
        
            typedef ::QStringList ( ::SireIO::PDB2::*formatSuffix_function_type)(  ) const;
            formatSuffix_function_type formatSuffix_function_value( &::SireIO::PDB2::formatSuffix );
            
            PDB2_exposer.def( 
                "formatSuffix"
                , formatSuffix_function_value
                , "Return the suffixes that these files are normally associated with" );
        
        }
        { //::SireIO::PDB2::isLead
        
            typedef bool ( ::SireIO::PDB2::*isLead_function_type)(  ) const;
            isLead_function_type isLead_function_value( &::SireIO::PDB2::isLead );
            
            PDB2_exposer.def( 
                "isLead"
                , isLead_function_value
                , "Return whether or not this is a lead parser. The lead parser is responsible\nfor starting the process of turning the parsed file into the System. There\nmust be one and one-only lead parser in a set of parsers creating a System" );
        
        }
        { //::SireIO::PDB2::nAtoms
        
            typedef int ( ::SireIO::PDB2::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::PDB2::nAtoms );
            
            PDB2_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "Return the total number of atoms." );
        
        }
        { //::SireIO::PDB2::nAtoms
        
            typedef int ( ::SireIO::PDB2::*nAtoms_function_type)( int ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::PDB2::nAtoms );
            
            PDB2_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("i") )
                , "Return the number of atoms in molecule i." );
        
        }
        { //::SireIO::PDB2::nChains
        
            typedef int ( ::SireIO::PDB2::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireIO::PDB2::nChains );
            
            PDB2_exposer.def( 
                "nChains"
                , nChains_function_value
                , "Return the total number of chains." );
        
        }
        { //::SireIO::PDB2::nChains
        
            typedef int ( ::SireIO::PDB2::*nChains_function_type)( int ) const;
            nChains_function_type nChains_function_value( &::SireIO::PDB2::nChains );
            
            PDB2_exposer.def( 
                "nChains"
                , nChains_function_value
                , ( bp::arg("i") )
                , "Return the number of chains in molecule i." );
        
        }
        { //::SireIO::PDB2::nMolecules
        
            typedef int ( ::SireIO::PDB2::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireIO::PDB2::nMolecules );
            
            PDB2_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , "Return the number of models (molecules)." );
        
        }
        { //::SireIO::PDB2::nResidues
        
            typedef int ( ::SireIO::PDB2::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireIO::PDB2::nResidues );
            
            PDB2_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , "Return the total number of residues." );
        
        }
        { //::SireIO::PDB2::nResidues
        
            typedef int ( ::SireIO::PDB2::*nResidues_function_type)( int ) const;
            nResidues_function_type nResidues_function_value( &::SireIO::PDB2::nResidues );
            
            PDB2_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , ( bp::arg("i") )
                , "Return the number of residues in molecule i." );
        
        }
        PDB2_exposer.def( bp::self != bp::self );
        { //::SireIO::PDB2::operator=
        
            typedef ::SireIO::PDB2 & ( ::SireIO::PDB2::*assign_function_type)( ::SireIO::PDB2 const & ) ;
            assign_function_type assign_function_value( &::SireIO::PDB2::operator= );
            
            PDB2_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        PDB2_exposer.def( bp::self == bp::self );
        { //::SireIO::PDB2::toLines
        
            typedef ::QVector< QString > ( ::SireIO::PDB2::*toLines_function_type)( bool ) const;
            toLines_function_type toLines_function_value( &::SireIO::PDB2::toLines );
            
            PDB2_exposer.def( 
                "toLines"
                , toLines_function_value
                , ( bp::arg("is_velocity")=(bool)(false) )
                , "Convert the parsed data to a collection of PDB record lines." );
        
        }
        { //::SireIO::PDB2::toString
        
            typedef ::QString ( ::SireIO::PDB2::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::PDB2::toString );
            
            PDB2_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this parser" );
        
        }
        { //::SireIO::PDB2::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::PDB2::typeName );
            
            PDB2_exposer.def( 
                "typeName"
                , typeName_function_value
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::PDB2::what
        
            typedef char const * ( ::SireIO::PDB2::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::PDB2::what );
            
            PDB2_exposer.def( 
                "what"
                , what_function_value
                , "Return the C++ name for this class" );
        
        }
        { //::SireIO::PDB2::writeVelocityFile
        
            typedef bool ( ::SireIO::PDB2::*writeVelocityFile_function_type)( ::QString const & ) const;
            writeVelocityFile_function_type writeVelocityFile_function_value( &::SireIO::PDB2::writeVelocityFile );
            
            PDB2_exposer.def( 
                "writeVelocityFile"
                , writeVelocityFile_function_value
                , ( bp::arg("filename") )
                , "Write a velocity file in PDB format. This can be used as a restart for NAMD simulations." );
        
        }
        PDB2_exposer.staticmethod( "typeName" );
        PDB2_exposer.def( "__copy__", &__copy__);
        PDB2_exposer.def( "__deepcopy__", &__copy__);
        PDB2_exposer.def( "clone", &__copy__);
        PDB2_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::PDB2 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDB2_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::PDB2 >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        PDB2_exposer.def_pickle(sire_pickle_suite< ::SireIO::PDB2 >());
        PDB2_exposer.def( "__str__", &__str__< ::SireIO::PDB2 > );
        PDB2_exposer.def( "__repr__", &__str__< ::SireIO::PDB2 > );
    }

}
