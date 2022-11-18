// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Molecule.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "chain.h"

#include "core.h"

#include "cutgroup.h"

#include "evaluator.h"

#include "molecule.h"

#include "moleculeinfo.h"

#include "moleditor.h"

#include "molviewproperty.h"

#include "mover.hpp"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "selectormol.h"

#include <QDebug>

#include "molecule.h"

SireMol::Molecule __copy__(const SireMol::Molecule &other){ return SireMol::Molecule(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Molecule_class(){

    { //::SireMol::Molecule
        typedef bp::class_< SireMol::Molecule, bp::bases< SireMol::MoleculeView, SireBase::Property > > Molecule_exposer_t;
        Molecule_exposer_t Molecule_exposer = Molecule_exposer_t( "Molecule", "A Molecule represents a complete molecule.\n\nMost of the manipulation of a molecule is handled by the orer classes,\ne.g. Mover, Selector, Editer, Evaluator.\n\nThese classes provide additional member functions, thereby allowing me\nto keep the API of Molecule small.\n\nExamples of use include;\n\nmol = mol.move().translate( Vector(1,2,3) )\npoint = mol.evaluate().center()\nmass = mol.evaluate().mass()\n\nmol = mol.edit().rename( ResNum(43)[0], ALA ).commit()\n\nEqually, we can quickly select well-defined subgroups within the\nmolecule, e.g. atom(s), residue(e), chain(s), CutGroup(s) and\nsegment(s), via the select functions, e.g.\n\nala49 = mol.select( ResName(ala) + ResNum(49) );\n\nor if there is more than one residue designate ALA:49\n\nala49_0 = mol.select( (ResName(ala)+ResNum(49))[0] );\n\nor to get all of these residues, do\n\nall_ala49 = mol.selectAll( ResName(ala) + ResNum(49) );\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Molecule_scope( Molecule_exposer );
        Molecule_exposer.def( bp::init< QString const & >(( bp::arg("molname") ), "Construct a new Molecule, called molname") );
        Molecule_exposer.def( bp::init< SireMol::MoleculeData const & >(( bp::arg("moldata") ), "Construct from the passed MoleculeData") );
        Molecule_exposer.def( bp::init< SireMol::Molecule const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Molecule::assertContainsMetadata
        
            typedef void ( ::SireMol::Molecule::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Molecule::assertContainsMetadata );
            
            Molecule_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this molecule contains some metadata at metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::assertContainsMetadata
        
            typedef void ( ::SireMol::Molecule::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Molecule::assertContainsMetadata );
            
            Molecule_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this molecule contains some metadata at metakey metakey\nassociated with the property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::assertContainsProperty
        
            typedef void ( ::SireMol::Molecule::*assertContainsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::Molecule::assertContainsProperty );
            
            Molecule_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Assert that this molecule contains a property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::edit
        
            typedef ::SireMol::MolEditor ( ::SireMol::Molecule::*edit_function_type)(  ) const;
            edit_function_type edit_function_value( &::SireMol::Molecule::edit );
            
            Molecule_exposer.def( 
                "edit"
                , edit_function_value
                , bp::release_gil_policy()
                , "Return an Editor that can edit any part of this molecule" );
        
        }
        { //::SireMol::Molecule::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Molecule::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Molecule::evaluate );
            
            Molecule_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return an Evaluator that evaluates values using\nall of the atoms in this molecule" );
        
        }
        { //::SireMol::Molecule::hasMetadata
        
            typedef bool ( ::SireMol::Molecule::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Molecule::hasMetadata );
            
            Molecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not this molecule posseses metadata with\nmetakey metakey" );
        
        }
        { //::SireMol::Molecule::hasMetadata
        
            typedef bool ( ::SireMol::Molecule::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Molecule::hasMetadata );
            
            Molecule_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not the property of this molecule at\nkey key has metadata at metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::hasProperty
        
            typedef bool ( ::SireMol::Molecule::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Molecule::hasProperty );
            
            Molecule_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return whether or not this molecule posseses a property at key key" );
        
        }
        { //::SireMol::Molecule::info
        
            typedef ::SireMol::MoleculeInfo ( ::SireMol::Molecule::*info_function_type)(  ) const;
            info_function_type info_function_value( &::SireMol::Molecule::info );
            
            Molecule_exposer.def( 
                "info"
                , info_function_value
                , bp::release_gil_policy()
                , "Return the MoleculeInfo object that holds information about the layout\nof the atoms, residues, chains and segments in the molecule" );
        
        }
        { //::SireMol::Molecule::invert
        
            typedef ::SireMol::SelectorMol ( ::SireMol::Molecule::*invert_function_type)(  ) const;
            invert_function_type invert_function_value( &::SireMol::Molecule::invert );
            
            Molecule_exposer.def( 
                "invert"
                , invert_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecule::isEmpty
        
            typedef bool ( ::SireMol::Molecule::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Molecule::isEmpty );
            
            Molecule_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is empty" );
        
        }
        { //::SireMol::Molecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::Molecule::*metadata_function_type)( ::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::Molecule::metadata );
            
            Molecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the metadata for the metakey metakey\nThrow: SireMol::missing_property\n" );
        
        }
        { //::SireMol::Molecule::metadata
        
            typedef ::SireBase::Property const & ( ::SireMol::Molecule::*metadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            metadata_function_type metadata_function_value( &::SireMol::Molecule::metadata );
            
            Molecule_exposer.def( 
                "metadata"
                , metadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the metadata for the metakey metakey for\nthe property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Molecule::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Molecule::metadataKeys );
            
            Molecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "Return the metakeys of all the metadata in this molecule" );
        
        }
        { //::SireMol::Molecule::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Molecule::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Molecule::metadataKeys );
            
            Molecule_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the metakeys for all of the metadata for the property\nat key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Molecule::move
        
            typedef ::SireMol::Mover< SireMol::Molecule > ( ::SireMol::Molecule::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Molecule::move );
            
            Molecule_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return a Mover that moves all of the atoms in\nthis molecule" );
        
        }
        { //::SireMol::Molecule::nAtoms
        
            typedef int ( ::SireMol::Molecule::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Molecule::nAtoms );
            
            Molecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this molecule" );
        
        }
        { //::SireMol::Molecule::nAtoms
        
            typedef int ( ::SireMol::Molecule::*nAtoms_function_type)( ::SireMol::ResID const & ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Molecule::nAtoms );
            
            Molecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return the number of atoms in the identified residue(s)" );
        
        }
        { //::SireMol::Molecule::nAtoms
        
            typedef int ( ::SireMol::Molecule::*nAtoms_function_type)( ::SireMol::ChainID const & ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Molecule::nAtoms );
            
            Molecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return the number of atoms in the identified chain(s)" );
        
        }
        { //::SireMol::Molecule::nAtoms
        
            typedef int ( ::SireMol::Molecule::*nAtoms_function_type)( ::SireMol::SegID const & ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Molecule::nAtoms );
            
            Molecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("seg") )
                , bp::release_gil_policy()
                , "Return the number of atoms in the identified segment(s)" );
        
        }
        { //::SireMol::Molecule::nAtoms
        
            typedef int ( ::SireMol::Molecule::*nAtoms_function_type)( ::SireMol::CGID const & ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Molecule::nAtoms );
            
            Molecule_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return the number of atoms in the identified cutgroups(s)" );
        
        }
        { //::SireMol::Molecule::nChains
        
            typedef int ( ::SireMol::Molecule::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::Molecule::nChains );
            
            Molecule_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "Return the number of chains in this molecule" );
        
        }
        { //::SireMol::Molecule::nCutGroups
        
            typedef int ( ::SireMol::Molecule::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::Molecule::nCutGroups );
            
            Molecule_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "Return the number of CutGroups in this molecule" );
        
        }
        { //::SireMol::Molecule::nResidues
        
            typedef int ( ::SireMol::Molecule::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::Molecule::nResidues );
            
            Molecule_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , bp::release_gil_policy()
                , "Return the number of residues in this molecule" );
        
        }
        { //::SireMol::Molecule::nResidues
        
            typedef int ( ::SireMol::Molecule::*nResidues_function_type)( ::SireMol::ChainID const & ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::Molecule::nResidues );
            
            Molecule_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , ( bp::arg("id") )
                , bp::release_gil_policy()
                , "Return the number of residues in the identified chain(s)" );
        
        }
        { //::SireMol::Molecule::nSegments
        
            typedef int ( ::SireMol::Molecule::*nSegments_function_type)(  ) const;
            nSegments_function_type nSegments_function_value( &::SireMol::Molecule::nSegments );
            
            Molecule_exposer.def( 
                "nSegments"
                , nSegments_function_value
                , bp::release_gil_policy()
                , "Return the number of segments in this molecule" );
        
        }
        { //::SireMol::Molecule::name
        
            typedef ::SireMol::MolName const & ( ::SireMol::Molecule::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireMol::Molecule::name );
            
            Molecule_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the name of this molecule" );
        
        }
        { //::SireMol::Molecule::number
        
            typedef ::SireMol::MolNum ( ::SireMol::Molecule::*number_function_type)(  ) const;
            number_function_type number_function_value( &::SireMol::Molecule::number );
            
            Molecule_exposer.def( 
                "number"
                , number_function_value
                , bp::release_gil_policy()
                , "Return the number of this molecule - this is used\nto identify the molecule" );
        
        }
        Molecule_exposer.def( bp::self != bp::self );
        { //::SireMol::Molecule::operator=
        
            typedef ::SireMol::Molecule & ( ::SireMol::Molecule::*assign_function_type)( ::SireMol::Molecule const & ) ;
            assign_function_type assign_function_value( &::SireMol::Molecule::operator= );
            
            Molecule_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Molecule_exposer.def( bp::self == bp::self );
        { //::SireMol::Molecule::properties
        
            typedef ::SireBase::Properties const & ( ::SireMol::Molecule::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireMol::Molecule::properties );
            
            Molecule_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the properties of this molecule" );
        
        }
        { //::SireMol::Molecule::property
        
            typedef ::SireBase::Property const & ( ::SireMol::Molecule::*property_function_type)( ::SireBase::PropertyName const & ) const;
            property_function_type property_function_value( &::SireMol::Molecule::property );
            
            Molecule_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("key") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the property associated with the key key\nThrow: SireMol::missing_property\n" );
        
        }
        { //::SireMol::Molecule::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Molecule::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Molecule::propertyKeys );
            
            Molecule_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all of the properties in this molecule" );
        
        }
        { //::SireMol::Molecule::selectedAll
        
            typedef bool ( ::SireMol::Molecule::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Molecule::selectedAll );
            
            Molecule_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a complete molecule" );
        
        }
        { //::SireMol::Molecule::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Molecule::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Molecule::selection );
            
            Molecule_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return which atoms are selected in this view" );
        
        }
        { //::SireMol::Molecule::selector
        
            typedef ::SireMol::SelectorMol ( ::SireMol::Molecule::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::Molecule::selector );
            
            Molecule_exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecule::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::Molecule::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMol::Molecule::toSelector );
            
            Molecule_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecule::toString
        
            typedef ::QString ( ::SireMol::Molecule::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Molecule::toString );
            
            Molecule_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this molecule" );
        
        }
        { //::SireMol::Molecule::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Molecule::typeName );
            
            Molecule_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecule::update
        
            typedef void ( ::SireMol::Molecule::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Molecule::update );
            
            Molecule_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update this molecule with the passed molecule data.\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMol::Molecule::version
        
            typedef ::quint64 ( ::SireMol::Molecule::*version_function_type)(  ) const;
            version_function_type version_function_value( &::SireMol::Molecule::version );
            
            Molecule_exposer.def( 
                "version"
                , version_function_value
                , bp::release_gil_policy()
                , "Return the version number of this molecule - all molecules\nwith the same ID number and version number must be identical" );
        
        }
        { //::SireMol::Molecule::version
        
            typedef ::quint64 ( ::SireMol::Molecule::*version_function_type)( ::SireBase::PropertyName const & ) const;
            version_function_type version_function_value( &::SireMol::Molecule::version );
            
            Molecule_exposer.def( 
                "version"
                , version_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the version number of the property at key key.\nAll molecules with the same ID number and same property version\nnumber must have the same value of this property\n(although this says nothing about any metadata associated\nwith this property)\nThrow: SireBase::missing_property\n" );
        
        }
        Molecule_exposer.staticmethod( "typeName" );
        Molecule_exposer.def( "__copy__", &__copy__);
        Molecule_exposer.def( "__deepcopy__", &__copy__);
        Molecule_exposer.def( "clone", &__copy__);
        Molecule_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Molecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Molecule_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Molecule >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Molecule_exposer.def_pickle(sire_pickle_suite< ::SireMol::Molecule >());
        Molecule_exposer.def( "__str__", &__str__< ::SireMol::Molecule > );
        Molecule_exposer.def( "__repr__", &__str__< ::SireMol::Molecule > );
        Molecule_exposer.def( "__len__", &__len_size< ::SireMol::Molecule > );
    }

}
