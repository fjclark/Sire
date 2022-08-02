// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Segment.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "evaluator.h"

#include "groupatomids.h"

#include "molecule.h"

#include "mover.hpp"

#include "mover_metaid.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "segment.h"

#include "segproperty.hpp"

const QString& get_Metadata_SireMol_SegStringProperty_function1(const SireMol::Segment &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_SegStringProperty_function2(const SireMol::Segment &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_SegIntProperty_function1(const SireMol::Segment &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_SegIntProperty_function2(const SireMol::Segment &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_SegFloatProperty_function1(const SireMol::Segment &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_SegFloatProperty_function2(const SireMol::Segment &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_SegVariantProperty_function1(const SireMol::Segment &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_SegVariantProperty_function2(const SireMol::Segment &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_SegPropertyProperty_function1(const SireMol::Segment &atom,
                                   const QString &metakey){ return atom.metadata< SireBase::PropertyPtr >(metakey); }

const SireBase::PropertyPtr& get_Metadata_SireMol_SegPropertyProperty_function2(const SireMol::Segment &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< SireBase::PropertyPtr >(key, metakey); }

SireMol::Segment __copy__(const SireMol::Segment &other){ return SireMol::Segment(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Segment_class(){

    { //::SireMol::Segment
        typedef bp::class_< SireMol::Segment, bp::bases< SireMol::MoleculeView, SireBase::Property > > Segment_exposer_t;
        Segment_exposer_t Segment_exposer = Segment_exposer_t( "Segment", "This is a view of a single segment within a molecule\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope Segment_scope( Segment_exposer );
        Segment_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::SegID const & >(( bp::arg("data"), bp::arg("segid") ), "Construct the Segment at ID cgid in the molecule whose data\nis in moldata\nThrow: SireMol::missing_Segment\nThrow: SireMol::duplicate_Segment\nThrow: SireError::invalid_index\n") );
        Segment_exposer.def( bp::init< SireMol::Segment const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Segment::assertContainsMetadata
        
            typedef void ( ::SireMol::Segment::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Segment::assertContainsMetadata );
            
            Segment_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that this segment has an SegProperty piece of metadata\nat metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Segment::assertContainsMetadata
        
            typedef void ( ::SireMol::Segment::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Segment::assertContainsMetadata );
            
            Segment_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Assert that the property at key key has an SegProperty\npiece of metadata at metakey metakey\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Segment::assertContainsProperty
        
            typedef void ( ::SireMol::Segment::*assertContainsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::Segment::assertContainsProperty );
            
            Segment_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Assert that this segment has an SegProperty at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Segment::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > const & ( ::SireMol::Segment::*atomIdxs_function_type)(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::Segment::atomIdxs );
            
            Segment_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the indicies of the atoms in this segment, in the\norder that they appear in this segment" );
        
        }
        { //::SireMol::Segment::contains
        
            typedef bool ( ::SireMol::Segment::*contains_function_type)( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::Segment::contains );
            
            Segment_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") )
                , bp::release_gil_policy()
                , "Return whether or not this segment contains the atom\nat index atomidx" );
        
        }
        { //::SireMol::Segment::contains
        
            typedef bool ( ::SireMol::Segment::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Segment::contains );
            
            Segment_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this segment contains all of\nthe atoms identified by the ID atomid" );
        
        }
        { //::SireMol::Segment::edit
        
            typedef ::SireMol::SegEditor ( ::SireMol::Segment::*edit_function_type)(  ) const;
            edit_function_type edit_function_value( &::SireMol::Segment::edit );
            
            Segment_exposer.def( 
                "edit"
                , edit_function_value
                , bp::release_gil_policy()
                , "Return an editor that can edit this Segment" );
        
        }
        { //::SireMol::Segment::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Segment::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Segment::evaluate );
            
            Segment_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , bp::release_gil_policy()
                , "Return an evaluator that can evaluate properties\nof this Segment" );
        
        }
        { //::SireMol::Segment::hasMetadata
        
            typedef bool ( ::SireMol::Segment::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Segment::hasMetadata );
            
            Segment_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether or not there is a SegProperty at metakey metakey" );
        
        }
        { //::SireMol::Segment::hasMetadata
        
            typedef bool ( ::SireMol::Segment::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Segment::hasMetadata );
            
            Segment_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") )
                , bp::release_gil_policy()
                , "Return whether the metadata at metakey metakey for the property\nat key key is a SegProperty\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Segment::hasProperty
        
            typedef bool ( ::SireMol::Segment::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Segment::hasProperty );
            
            Segment_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return whether or not there is a SegProperty at key key" );
        
        }
        { //::SireMol::Segment::index
        
            typedef ::SireMol::SegIdx ( ::SireMol::Segment::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireMol::Segment::index );
            
            Segment_exposer.def( 
                "index"
                , index_function_value
                , bp::release_gil_policy()
                , "Return the index of this Segment in the molecule" );
        
        }
        { //::SireMol::Segment::intersects
        
            typedef bool ( ::SireMol::Segment::*intersects_function_type)( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Segment::intersects );
            
            Segment_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return whether or not this segment contains some of\nthe atoms identified by the ID atomid" );
        
        }
        { //::SireMol::Segment::isEmpty
        
            typedef bool ( ::SireMol::Segment::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Segment::isEmpty );
            
            Segment_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether or not this segment is empty" );
        
        }
        { //::SireMol::Segment::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Segment::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Segment::metadataKeys );
            
            Segment_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , bp::release_gil_policy()
                , "Return the metakeys of all SegProperty metadata" );
        
        }
        { //::SireMol::Segment::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Segment::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Segment::metadataKeys );
            
            Segment_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the metakeys of all SegProperty metadata for\nthe property at key key\nThrow: SireBase::missing_property\n" );
        
        }
        { //::SireMol::Segment::move
        
            typedef ::SireMol::Mover< SireMol::Segment > ( ::SireMol::Segment::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Segment::move );
            
            Segment_exposer.def( 
                "move"
                , move_function_value
                , bp::release_gil_policy()
                , "Return an object that can move a copy of this Segment" );
        
        }
        { //::SireMol::Segment::nAtoms
        
            typedef int ( ::SireMol::Segment::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Segment::nAtoms );
            
            Segment_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this Segment" );
        
        }
        { //::SireMol::Segment::name
        
            typedef ::SireMol::SegName const & ( ::SireMol::Segment::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireMol::Segment::name );
            
            Segment_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the name of this Segment" );
        
        }
        { //::SireMol::Segment::number
        
            typedef ::SireMol::SegIdx ( ::SireMol::Segment::*number_function_type)(  ) const;
            number_function_type number_function_value( &::SireMol::Segment::number );
            
            Segment_exposer.def( 
                "number"
                , number_function_value
                , bp::release_gil_policy()
                , "Return the number of this segment (same as its index)" );
        
        }
        Segment_exposer.def( bp::self != bp::self );
        { //::SireMol::Segment::operator=
        
            typedef ::SireMol::Segment & ( ::SireMol::Segment::*assign_function_type)( ::SireMol::Segment const & ) ;
            assign_function_type assign_function_value( &::SireMol::Segment::operator= );
            
            Segment_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Segment_exposer.def( bp::self == bp::self );
        { //::SireMol::Segment::propertyAsProperty
        
            typedef ::SireBase::PropertyPtr ( ::SireMol::Segment::*propertyAsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            propertyAsProperty_function_type propertyAsProperty_function_value( &::SireMol::Segment::propertyAsProperty );
            
            Segment_exposer.def( 
                "propertyAsProperty"
                , propertyAsProperty_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the specified property as a PropertyPtr" );
        
        }
        { //::SireMol::Segment::propertyAsVariant
        
            typedef ::QVariant ( ::SireMol::Segment::*propertyAsVariant_function_type)( ::SireBase::PropertyName const & ) const;
            propertyAsVariant_function_type propertyAsVariant_function_value( &::SireMol::Segment::propertyAsVariant );
            
            Segment_exposer.def( 
                "propertyAsVariant"
                , propertyAsVariant_function_value
                , ( bp::arg("key") )
                , bp::release_gil_policy()
                , "Return the specified property as a QVariant" );
        
        }
        { //::SireMol::Segment::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Segment::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Segment::propertyKeys );
            
            Segment_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value
                , bp::release_gil_policy()
                , "Return the keys of all SegProperty properties" );
        
        }
        { //::SireMol::Segment::selectedAll
        
            typedef bool ( ::SireMol::Segment::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Segment::selectedAll );
            
            Segment_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not this segment contains the entire molecule" );
        
        }
        { //::SireMol::Segment::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Segment::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Segment::selection );
            
            Segment_exposer.def( 
                "selection"
                , selection_function_value
                , bp::release_gil_policy()
                , "Return the atoms that are in this Segment" );
        
        }
        { //::SireMol::Segment::selector
        
            typedef ::SireMol::Selector< SireMol::Segment > ( ::SireMol::Segment::*selector_function_type)(  ) const;
            selector_function_type selector_function_value( &::SireMol::Segment::selector );
            
            Segment_exposer.def( 
                "selector"
                , selector_function_value
                , bp::release_gil_policy()
                , "Return a selector that can be used to change the selection\nof segments from the molecule" );
        
        }
        { //::SireMol::Segment::toSelector
        
            typedef ::SireMol::MolViewPtr ( ::SireMol::Segment::*toSelector_function_type)(  ) const;
            toSelector_function_type toSelector_function_value( &::SireMol::Segment::toSelector );
            
            Segment_exposer.def( 
                "toSelector"
                , toSelector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Segment::toString
        
            typedef ::QString ( ::SireMol::Segment::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Segment::toString );
            
            Segment_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this segment" );
        
        }
        { //::SireMol::Segment::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Segment::typeName );
            
            Segment_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Segment::update
        
            typedef void ( ::SireMol::Segment::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Segment::update );
            
            Segment_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update this segment with the passed molecule data.\nThrow: SireError::incompatible_error\n" );
        
        }
        Segment_exposer.staticmethod( "typeName" );
        Segment_exposer.def( "_get_property_SireMol_SegStringProperty", &SireMol::Segment::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegStringProperty", get_Metadata_SireMol_SegStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegStringProperty", &get_Metadata_SireMol_SegStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_property_SireMol_SegIntProperty", &SireMol::Segment::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegIntProperty", get_Metadata_SireMol_SegIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegIntProperty", &get_Metadata_SireMol_SegIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_property_SireMol_SegFloatProperty", &SireMol::Segment::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegFloatProperty", get_Metadata_SireMol_SegFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegFloatProperty", &get_Metadata_SireMol_SegFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_property_SireMol_SegVariantProperty", &SireMol::Segment::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegVariantProperty", get_Metadata_SireMol_SegVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegVariantProperty", &get_Metadata_SireMol_SegVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_property_SireMol_SegPropertyProperty", &SireMol::Segment::property< SireBase::PropertyPtr >, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegPropertyProperty", get_Metadata_SireMol_SegPropertyProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "_get_metadata_SireMol_SegPropertyProperty", &get_Metadata_SireMol_SegPropertyProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Segment_exposer.def( "__copy__", &__copy__);
        Segment_exposer.def( "__deepcopy__", &__copy__);
        Segment_exposer.def( "clone", &__copy__);
        Segment_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Segment >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Segment_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Segment >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Segment_exposer.def_pickle(sire_pickle_suite< ::SireMol::Segment >());
        Segment_exposer.def( "__str__", &__str__< ::SireMol::Segment > );
        Segment_exposer.def( "__repr__", &__str__< ::SireMol::Segment > );
        Segment_exposer.def( "__len__", &__len_size< ::SireMol::Segment > );
    }

}
