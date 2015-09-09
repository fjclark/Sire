// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Bead.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "bead.h"

#include "beadeditor.h"

#include "beads.h"

#include "mover.hpp"

#include "tostring.h"

#include "bead.h"

#include "beadproperty.hpp"

const QString& get_Metadata_SireMol_BeadStringProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< QString >(metakey); }

const QString& get_Metadata_SireMol_BeadStringProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QString >(key, metakey); }

const qint64& get_Metadata_SireMol_BeadIntProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< qint64 >(metakey); }

const qint64& get_Metadata_SireMol_BeadIntProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< qint64 >(key, metakey); }

const double& get_Metadata_SireMol_BeadFloatProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< double >(metakey); }

const double& get_Metadata_SireMol_BeadFloatProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< double >(key, metakey); }

const QVariant& get_Metadata_SireMol_BeadVariantProperty_function1(const SireMol::Bead &atom,
                                   const QString &metakey){ return atom.metadata< QVariant >(metakey); }

const QVariant& get_Metadata_SireMol_BeadVariantProperty_function2(const SireMol::Bead &atom,
                                   const QString &key, const QString &metakey){
                                        return atom.metadata< QVariant >(key, metakey); }

SireMol::Bead __copy__(const SireMol::Bead &other){ return SireMol::Bead(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_Bead_class(){

    { //::SireMol::Bead
        typedef bp::class_< SireMol::Bead, bp::bases< SireMol::MoleculeView, SireBase::Property > > Bead_exposer_t;
        Bead_exposer_t Bead_exposer = Bead_exposer_t( "Bead", bp::init< >() );
        bp::scope Bead_scope( Bead_exposer );
        Bead_exposer.def( bp::init< SireMol::MoleculeData const &, SireMol::BeadIdx const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("moldata"), bp::arg("bead"), bp::arg("map")=SireBase::PropertyMap() )) );
        Bead_exposer.def( bp::init< SireMol::Bead const & >(( bp::arg("other") )) );
        { //::SireMol::Bead::assertContainsMetadata
        
            typedef void ( ::SireMol::Bead::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Bead::assertContainsMetadata );
            
            Bead_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::Bead::assertContainsMetadata
        
            typedef void ( ::SireMol::Bead::*assertContainsMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            assertContainsMetadata_function_type assertContainsMetadata_function_value( &::SireMol::Bead::assertContainsMetadata );
            
            Bead_exposer.def( 
                "assertContainsMetadata"
                , assertContainsMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::Bead::assertContainsProperty
        
            typedef void ( ::SireMol::Bead::*assertContainsProperty_function_type)( ::SireBase::PropertyName const & ) const;
            assertContainsProperty_function_type assertContainsProperty_function_value( &::SireMol::Bead::assertContainsProperty );
            
            Bead_exposer.def( 
                "assertContainsProperty"
                , assertContainsProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Bead::at
        
            typedef ::SireMol::Atom ( ::SireMol::Bead::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::Bead::at );
            
            Bead_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Bead::atom
        
            typedef ::SireMol::Atom ( ::SireMol::Bead::*atom_function_type)( int ) const;
            atom_function_type atom_function_value( &::SireMol::Bead::atom );
            
            Bead_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Bead::atomIdxs
        
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::Bead::*atomIdxs_function_type)(  ) const;
            atomIdxs_function_type atomIdxs_function_value( &::SireMol::Bead::atomIdxs );
            
            Bead_exposer.def( 
                "atomIdxs"
                , atomIdxs_function_value );
        
        }
        { //::SireMol::Bead::beading
        
            typedef ::SireMol::Beading const & ( ::SireMol::Bead::*beading_function_type)(  ) const;
            beading_function_type beading_function_value( &::SireMol::Bead::beading );
            
            Bead_exposer.def( 
                "beading"
                , beading_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::Bead::beads
        
            typedef ::SireMol::Beads ( ::SireMol::Bead::*beads_function_type)(  ) const;
            beads_function_type beads_function_value( &::SireMol::Bead::beads );
            
            Bead_exposer.def( 
                "beads"
                , beads_function_value );
        
        }
        { //::SireMol::Bead::contains
        
            typedef bool ( ::SireMol::Bead::*contains_function_type)( ::SireMol::AtomIdx ) const;
            contains_function_type contains_function_value( &::SireMol::Bead::contains );
            
            Bead_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomidx") ) );
        
        }
        { //::SireMol::Bead::contains
        
            typedef bool ( ::SireMol::Bead::*contains_function_type)( ::SireMol::AtomID const & ) const;
            contains_function_type contains_function_value( &::SireMol::Bead::contains );
            
            Bead_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::Bead::count
        
            typedef int ( ::SireMol::Bead::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::Bead::count );
            
            Bead_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::Bead::edit
        
            typedef ::SireMol::BeadEditor ( ::SireMol::Bead::*edit_function_type)(  ) const;
            edit_function_type edit_function_value( &::SireMol::Bead::edit );
            
            Bead_exposer.def( 
                "edit"
                , edit_function_value );
        
        }
        { //::SireMol::Bead::evaluate
        
            typedef ::SireMol::Evaluator ( ::SireMol::Bead::*evaluate_function_type)(  ) const;
            evaluate_function_type evaluate_function_value( &::SireMol::Bead::evaluate );
            
            Bead_exposer.def( 
                "evaluate"
                , evaluate_function_value );
        
        }
        { //::SireMol::Bead::hasMetadata
        
            typedef bool ( ::SireMol::Bead::*hasMetadata_function_type)( ::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Bead::hasMetadata );
            
            Bead_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("metakey") ) );
        
        }
        { //::SireMol::Bead::hasMetadata
        
            typedef bool ( ::SireMol::Bead::*hasMetadata_function_type)( ::SireBase::PropertyName const &,::SireBase::PropertyName const & ) const;
            hasMetadata_function_type hasMetadata_function_value( &::SireMol::Bead::hasMetadata );
            
            Bead_exposer.def( 
                "hasMetadata"
                , hasMetadata_function_value
                , ( bp::arg("key"), bp::arg("metakey") ) );
        
        }
        { //::SireMol::Bead::hasProperty
        
            typedef bool ( ::SireMol::Bead::*hasProperty_function_type)( ::SireBase::PropertyName const & ) const;
            hasProperty_function_type hasProperty_function_value( &::SireMol::Bead::hasProperty );
            
            Bead_exposer.def( 
                "hasProperty"
                , hasProperty_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Bead::index
        
            typedef ::SireMol::BeadIdx ( ::SireMol::Bead::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireMol::Bead::index );
            
            Bead_exposer.def( 
                "index"
                , index_function_value );
        
        }
        { //::SireMol::Bead::intersects
        
            typedef bool ( ::SireMol::Bead::*intersects_function_type)( ::SireMol::AtomID const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Bead::intersects );
            
            Bead_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("atomid") ) );
        
        }
        { //::SireMol::Bead::isEmpty
        
            typedef bool ( ::SireMol::Bead::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Bead::isEmpty );
            
            Bead_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::Bead::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*metadataKeys_function_type)(  ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Bead::metadataKeys );
            
            Bead_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value );
        
        }
        { //::SireMol::Bead::metadataKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*metadataKeys_function_type)( ::SireBase::PropertyName const & ) const;
            metadataKeys_function_type metadataKeys_function_value( &::SireMol::Bead::metadataKeys );
            
            Bead_exposer.def( 
                "metadataKeys"
                , metadataKeys_function_value
                , ( bp::arg("key") ) );
        
        }
        { //::SireMol::Bead::move
        
            typedef ::SireMol::Mover< SireMol::Bead > ( ::SireMol::Bead::*move_function_type)(  ) const;
            move_function_type move_function_value( &::SireMol::Bead::move );
            
            Bead_exposer.def( 
                "move"
                , move_function_value );
        
        }
        { //::SireMol::Bead::nAtoms
        
            typedef int ( ::SireMol::Bead::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::Bead::nAtoms );
            
            Bead_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        Bead_exposer.def( bp::self != bp::self );
        { //::SireMol::Bead::operator=
        
            typedef ::SireMol::Bead & ( ::SireMol::Bead::*assign_function_type)( ::SireMol::Bead const & ) ;
            assign_function_type assign_function_value( &::SireMol::Bead::operator= );
            
            Bead_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Bead_exposer.def( bp::self == bp::self );
        { //::SireMol::Bead::operator[]
        
            typedef ::SireMol::Atom ( ::SireMol::Bead::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Bead::operator[] );
            
            Bead_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireMol::Bead::propertyKeys
        
            typedef ::QStringList ( ::SireMol::Bead::*propertyKeys_function_type)(  ) const;
            propertyKeys_function_type propertyKeys_function_value( &::SireMol::Bead::propertyKeys );
            
            Bead_exposer.def( 
                "propertyKeys"
                , propertyKeys_function_value );
        
        }
        { //::SireMol::Bead::selectedAll
        
            typedef bool ( ::SireMol::Bead::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::Bead::selectedAll );
            
            Bead_exposer.def( 
                "selectedAll"
                , selectedAll_function_value );
        
        }
        { //::SireMol::Bead::selection
        
            typedef ::SireMol::AtomSelection ( ::SireMol::Bead::*selection_function_type)(  ) const;
            selection_function_type selection_function_value( &::SireMol::Bead::selection );
            
            Bead_exposer.def( 
                "selection"
                , selection_function_value );
        
        }
        { //::SireMol::Bead::size
        
            typedef int ( ::SireMol::Bead::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::Bead::size );
            
            Bead_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::Bead::toString
        
            typedef ::QString ( ::SireMol::Bead::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Bead::toString );
            
            Bead_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::Bead::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Bead::typeName );
            
            Bead_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::Bead::update
        
            typedef void ( ::SireMol::Bead::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Bead::update );
            
            Bead_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") ) );
        
        }
        Bead_exposer.staticmethod( "typeName" );
        Bead_exposer.def( "_get_property_SireMol_BeadStringProperty", &SireMol::Bead::property< QString >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadStringProperty", get_Metadata_SireMol_BeadStringProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadStringProperty", &get_Metadata_SireMol_BeadStringProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadIntProperty", &SireMol::Bead::property< qint64 >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadIntProperty", get_Metadata_SireMol_BeadIntProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadIntProperty", &get_Metadata_SireMol_BeadIntProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadFloatProperty", &SireMol::Bead::property< double >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadFloatProperty", get_Metadata_SireMol_BeadFloatProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadFloatProperty", &get_Metadata_SireMol_BeadFloatProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_property_SireMol_BeadVariantProperty", &SireMol::Bead::property< QVariant >, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadVariantProperty", get_Metadata_SireMol_BeadVariantProperty_function1, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "_get_metadata_SireMol_BeadVariantProperty", &get_Metadata_SireMol_BeadVariantProperty_function2, bp::return_value_policy<bp::copy_const_reference>());
        Bead_exposer.def( "__copy__", &__copy__);
        Bead_exposer.def( "__deepcopy__", &__copy__);
        Bead_exposer.def( "clone", &__copy__);
        Bead_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Bead >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bead_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Bead >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Bead_exposer.def( "__str__", &__str__< ::SireMol::Bead > );
        Bead_exposer.def( "__repr__", &__str__< ::SireMol::Bead > );
        Bead_exposer.def( "__len__", &__len_size< ::SireMol::Bead > );
    }

}
