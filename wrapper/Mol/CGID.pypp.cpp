// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CGID.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "atom.h"

#include "cgid.h"

#include "cgidentifier.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "groupgroupids.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "cgid.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CGID_class(){

    { //::SireMol::CGID
        typedef bp::class_< SireMol::CGID, bp::bases< SireID::ID >, boost::noncopyable > CGID_exposer_t;
        CGID_exposer_t CGID_exposer = CGID_exposer_t( "CGID", "This is the base class of all identifiers that are used\nto identify a CutGroup\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope CGID_scope( CGID_exposer );
        { //::SireMol::CGID::any
        
            typedef ::SireID::MatchAll< SireMol::CGID > ( *any_function_type )(  );
            any_function_type any_function_value( &::SireMol::CGID::any );
            
            CGID_exposer.def( 
                "any"
                , any_function_value
                , bp::release_gil_policy()
                , "Return a match for any cutgroup" );
        
        }
        { //::SireMol::CGID::atom
        
            typedef ::SireMol::AtomsIn< SireMol::CGID > ( ::SireMol::CGID::*atom_function_type)( int ) const;
            atom_function_type atom_function_value( &::SireMol::CGID::atom );
            
            CGID_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "Return a specific atom in the matching residues" );
        
        }
        { //::SireMol::CGID::atoms
        
            typedef ::SireMol::AtomsIn< SireMol::CGID > ( ::SireMol::CGID::*atoms_function_type)(  ) const;
            atoms_function_type atoms_function_value( &::SireMol::CGID::atoms );
            
            CGID_exposer.def( 
                "atoms"
                , atoms_function_value
                , bp::release_gil_policy()
                , "Return the atoms in the matching residues" );
        
        }
        { //::SireMol::CGID::atoms
        
            typedef ::SireMol::AtomsIn< SireMol::CGID > ( ::SireMol::CGID::*atoms_function_type)( int,int ) const;
            atoms_function_type atoms_function_value( &::SireMol::CGID::atoms );
            
            CGID_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "Return a range of atoms in the matching residues" );
        
        }
        { //::SireMol::CGID::fromString
        
            typedef ::SireMol::CGIdentifier ( *fromString_function_type )( ::QString const & );
            fromString_function_type fromString_function_value( &::SireMol::CGID::fromString );
            
            CGID_exposer.def( 
                "fromString"
                , fromString_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return an CGID constructed from the passed string" );
        
        }
        { //::SireMol::CGID::inverse
        
            typedef ::SireID::InvertMatch< SireMol::CGID > ( ::SireMol::CGID::*inverse_function_type)(  ) const;
            inverse_function_type inverse_function_value( &::SireMol::CGID::inverse );
            
            CGID_exposer.def( 
                "inverse"
                , inverse_function_value
                , bp::release_gil_policy()
                , "Return the inverse of this match" );
        
        }
        { //::SireMol::CGID::invert
        
            typedef ::SireID::InvertMatch< SireMol::CGID > ( ::SireMol::CGID::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::CGID::invert );
            
            CGID_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "Return the inverse of this match" );
        
        }
        { //::SireMol::CGID::map
        
            typedef ::QList< SireMol::CGIdx > ( ::SireMol::CGID::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::CGID::map );
            
            CGID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "Map this ID back to the indicies of the CutGroups\nwithin the molecule described by the info in molinfo" );
        
        }
        { //::SireMol::CGID::map
        
            typedef ::QList< SireMol::CGIdx > ( ::SireMol::CGID::*map_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            map_function_type map_function_value( &::SireMol::CGID::map );
            
            CGID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "Map this CGICutGroupD to the CutGroups in the passed molecule view\nThrow: SireMol::missing_cutgroup\nThrow: SireError::invalid_index\n" );
        
        }
        CGID_exposer.def( !bp::self );
        CGID_exposer.def( bp::self & bp::self );
        CGID_exposer.def( bp::self & bp::other< SireMol::AtomID >() );
        CGID_exposer.def( bp::self & bp::other< SireMol::SegID >() );
        CGID_exposer.def( bp::self & bp::other< SireMol::ChainID >() );
        CGID_exposer.def( bp::self & bp::other< SireMol::ResID >() );
        { //::SireMol::CGID::operator()
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__call___function_type)( ::SireBase::Range const & ) const;
            __call___function_type __call___function_value( &::SireMol::CGID::operator() );
            
            CGID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("range") )
                , "" );
        
        }
        { //::SireMol::CGID::operator()
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__call___function_type)( ::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::CGID::operator() );
            
            CGID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::CGID::operator()
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__call___function_type)( ::qint64,::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::CGID::operator() );
            
            CGID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("start"), bp::arg("end") )
                , "" );
        
        }
        { //::SireMol::CGID::operator()
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__call___function_type)( ::qint64,::qint64,::qint64 ) const;
            __call___function_type __call___function_value( &::SireMol::CGID::operator() );
            
            CGID_exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("start"), bp::arg("end"), bp::arg("increment") )
                , "" );
        
        }
        CGID_exposer.def( bp::self * bp::self );
        CGID_exposer.def( bp::self * bp::other< SireMol::AtomID >() );
        CGID_exposer.def( bp::self + bp::self );
        CGID_exposer.def( bp::self + bp::other< SireMol::AtomID >() );
        CGID_exposer.def( bp::self + bp::other< SireMol::SegID >() );
        CGID_exposer.def( bp::self + bp::other< SireMol::ChainID >() );
        CGID_exposer.def( bp::self + bp::other< SireMol::ResID >() );
        CGID_exposer.def( bp::self - bp::self );
        CGID_exposer.def( bp::self - bp::other< SireMol::AtomID >() );
        CGID_exposer.def( bp::self - bp::other< SireMol::SegID >() );
        CGID_exposer.def( bp::self - bp::other< SireMol::ChainID >() );
        CGID_exposer.def( bp::self - bp::other< SireMol::ResID >() );
        CGID_exposer.def( -bp::self );
        { //::SireMol::CGID::operator[]
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__getitem___function_type)( ::qint64 ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGID::operator[] );
            
            CGID_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , "" );
        
        }
        { //::SireMol::CGID::operator[]
        
            typedef ::SireID::Specify< SireMol::CGID > ( ::SireMol::CGID::*__getitem___function_type)( ::SireBase::Range const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::CGID::operator[] );
            
            CGID_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("range") )
                , "" );
        
        }
        CGID_exposer.def( bp::self | bp::self );
        CGID_exposer.def( bp::self | bp::other< SireMol::AtomID >() );
        { //::SireMol::CGID::selectAllFrom
        
            typedef ::SireMol::Selector< SireMol::CutGroup > ( ::SireMol::CGID::*selectAllFrom_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::CGID::selectAllFrom );
            
            CGID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "Select all the CutGroups from the passed view that match this ID\nThrow: SireMol::missing_cutgroup\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::CutGroup > > ( ::SireMol::CGID::*selectAllFrom_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::CGID::selectAllFrom );
            
            CGID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return all of the CutGroups from the molecules that match\nthis ID\nThrow: SireMol::missing_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::CutGroup > > ( ::SireMol::CGID::*selectAllFrom_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::CGID::selectAllFrom );
            
            CGID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the atoms from the molecule group molgroup that match\nthis ID\nThrow: SireMol::missing_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectAllFrom
        
            typedef ::QHash< SireMol::MolNum, SireMol::Selector< SireMol::CutGroup > > ( ::SireMol::CGID::*selectAllFrom_function_type)( ::SireMol::MolGroupsBase const &,::SireBase::PropertyMap const & ) const;
            selectAllFrom_function_type selectAllFrom_function_value( &::SireMol::CGID::selectAllFrom );
            
            CGID_exposer.def( 
                "selectAllFrom"
                , selectAllFrom_function_value
                , ( bp::arg("molgroups"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the set of atoms that match this ID in the molecule groups\nset molgroups\nThrow: SireMol::missing_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectFrom
        
            typedef ::SireMol::CutGroup ( ::SireMol::CGID::*selectFrom_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::CGID::selectFrom );
            
            CGID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molview"), bp::arg("map")=SireBase::PropertyMap() )
                , "Select the CutGroup from the passed view that matches this ID\nThrow: SireMol::missing_cutgroup\nThrow: SireError::invalid_index\nThrow: SireMol::duplicate_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectFrom
        
            typedef ::SireMol::CutGroup ( ::SireMol::CGID::*selectFrom_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::CGID::selectFrom );
            
            CGID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molecules"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the atom from the molecules molecules that matches\nthis ID\nThrow: SireMol::missing_cutgroup\nThrow: SireMol::duplicate_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectFrom
        
            typedef ::SireMol::CutGroup ( ::SireMol::CGID::*selectFrom_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::CGID::selectFrom );
            
            CGID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the atom from the molecule group molgroup that matches\nthis ID\nThrow: SireMol::missing_cutgroup\nThrow: SireMol::duplicate_cutgroup\n" );
        
        }
        { //::SireMol::CGID::selectFrom
        
            typedef ::SireMol::CutGroup ( ::SireMol::CGID::*selectFrom_function_type)( ::SireMol::MolGroupsBase const &,::SireBase::PropertyMap const & ) const;
            selectFrom_function_type selectFrom_function_value( &::SireMol::CGID::selectFrom );
            
            CGID_exposer.def( 
                "selectFrom"
                , selectFrom_function_value
                , ( bp::arg("molgroups"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the atom from the molecule groups molgroups that matches\nthis ID\nThrow: SireMol::missing_cutgroup\nThrow: SireMol::duplicate_cutgroup\n" );
        
        }
        { //::SireMol::CGID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGID::typeName );
            
            CGID_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGID_exposer.staticmethod( "any" );
        CGID_exposer.staticmethod( "fromString" );
        CGID_exposer.staticmethod( "typeName" );
        CGID_exposer.def( "__str__", &__str__< ::SireMol::CGID > );
        CGID_exposer.def( "__repr__", &__str__< ::SireMol::CGID > );
        CGID_exposer.def( "__hash__", &::SireMol::CGID::hash );
    }

}
