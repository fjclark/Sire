// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ChemicalBondHunter.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/coordgroup.h"

#include "atom.h"

#include "atomcoords.h"

#include "atomelements.h"

#include "atomselection.h"

#include "bondhunter.h"

#include "connectivity.h"

#include "molecule.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "moleculeview.h"

#include "mover.hpp"

#include "selector.hpp"

#include <QDebug>

#include <QMutex>

#include "bondhunter.h"

SireMol::ChemicalBondHunter __copy__(const SireMol::ChemicalBondHunter &other){ return SireMol::ChemicalBondHunter(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChemicalBondHunter_class(){

    { //::SireMol::ChemicalBondHunter
        typedef bp::class_< SireMol::ChemicalBondHunter, bp::bases< SireMol::CovalentBondHunter, SireMol::BondHunter, SireBase::Property > > ChemicalBondHunter_exposer_t;
        ChemicalBondHunter_exposer_t ChemicalBondHunter_exposer = ChemicalBondHunter_exposer_t( "ChemicalBondHunter", "This is a bond hunter that hunts for bonds using the distance between\natoms (and comparing that distance against the sum of atomic covalent\nradii), but then it runs over each atom and ensures that the atom does\nnot contain too many bonds. If it does, then only the n closest bonds\nare retained.\n\nAuthor: Christopher Woods\n", bp::init< >("Construct with the default tolerance") );
        bp::scope ChemicalBondHunter_scope( ChemicalBondHunter_exposer );
        ChemicalBondHunter_exposer.def( bp::init< double >(( bp::arg("tolerance") ), "Construct with the specified tolerance") );
        ChemicalBondHunter_exposer.def( bp::init< SireMol::ChemicalBondHunter const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::ChemicalBondHunter::operator()
        
            typedef ::SireMol::Connectivity ( ::SireMol::ChemicalBondHunter::*__call___function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            __call___function_type __call___function_value( &::SireMol::ChemicalBondHunter::operator() );
            
            ChemicalBondHunter_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::ChemicalBondHunter::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChemicalBondHunter::typeName );
            
            ChemicalBondHunter_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        ChemicalBondHunter_exposer.staticmethod( "typeName" );
        ChemicalBondHunter_exposer.def( "__copy__", &__copy__);
        ChemicalBondHunter_exposer.def( "__deepcopy__", &__copy__);
        ChemicalBondHunter_exposer.def( "clone", &__copy__);
        ChemicalBondHunter_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChemicalBondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChemicalBondHunter_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChemicalBondHunter >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChemicalBondHunter_exposer.def_pickle(sire_pickle_suite< ::SireMol::ChemicalBondHunter >());
        ChemicalBondHunter_exposer.def( "__str__", &__str__< ::SireMol::ChemicalBondHunter > );
        ChemicalBondHunter_exposer.def( "__repr__", &__str__< ::SireMol::ChemicalBondHunter > );
    }

}
