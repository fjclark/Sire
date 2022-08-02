// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomPropertyList.pypp.hpp"

namespace bp = boost::python;

#include "atompropertylist.h"

#include "atompropertylist.h"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

SireMol::AtomProperty<SireBase::PropertyList> __copy__(const SireMol::AtomProperty<SireBase::PropertyList> &other){ return SireMol::AtomProperty<SireBase::PropertyList>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_AtomPropertyList_class(){

    { //::SireMol::AtomProperty< SireBase::PropertyList >
        typedef bp::class_< SireMol::AtomProperty< SireBase::PropertyList >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomPropertyList_exposer_t;
        AtomPropertyList_exposer_t AtomPropertyList_exposer = AtomPropertyList_exposer_t( "AtomPropertyList", "", bp::init< >("") );
        bp::scope AtomPropertyList_scope( AtomPropertyList_exposer );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeInfo const &, SireBase::PropertyList const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeView const &, SireBase::PropertyList const & >(( bp::arg("molview"), bp::arg("default_value") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireBase::PropertyList const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireBase::PropertyList const & >(( bp::arg("value") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireBase::PackedArray2D< SireBase::PropertyList > const & >(( bp::arg("values") ), "") );
        AtomPropertyList_exposer.def( bp::init< SireMol::AtomProperty< SireBase::PropertyList > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomProperty< SireBase::PropertyList >::array
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::PropertyList > const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::array );
            
            AtomPropertyList_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::PropertyList >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::assertCanConvert );
            
            AtomPropertyList_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::assignFrom
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::PropertyList >::*assignFrom_function_type)( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::assignFrom );
            
            AtomPropertyList_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::at
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::PropertyList >::Array const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*at_function_type)( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::at );
            
            AtomPropertyList_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::at
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyList const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*at_function_type)( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::at );
            
            AtomPropertyList_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::canConvert
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::PropertyList >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::canConvert );
            
            AtomPropertyList_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::copyFrom
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::PropertyList >::*copyFrom_function_type)( ::QVector< SireBase::PropertyList > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::copyFrom );
            
            AtomPropertyList_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::copyFrom
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireBase::PropertyList >::*copyFrom_function_type)( ::QVector< SireBase::PropertyList > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::copyFrom );
            
            AtomPropertyList_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::count
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::PropertyList >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::count );
            
            AtomPropertyList_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::divide
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::PropertyList >::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::divide );
            
            AtomPropertyList_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::divideByResidue
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::PropertyList >::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::divideByResidue );
            
            AtomPropertyList_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::fromVariant
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::PropertyList > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::fromVariant );
            
            AtomPropertyList_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::get
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::PropertyList >::Array const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::get );
            
            AtomPropertyList_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::get
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyList const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::get );
            
            AtomPropertyList_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::getAsProperty
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::PropertyList >::*getAsProperty_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::getAsProperty );
            
            AtomPropertyList_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::getAsVariant
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::QVariant ( ::SireMol::AtomProperty< SireBase::PropertyList >::*getAsVariant_function_type)( ::SireMol::CGAtomIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::getAsVariant );
            
            AtomPropertyList_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("cgatomidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::PropertyList >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::isCompatibleWith );
            
            AtomPropertyList_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::PropertyList >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfo const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::isCompatibleWith );
            
            AtomPropertyList_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::isEmpty
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireBase::PropertyList >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::isEmpty );
            
            AtomPropertyList_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::matchToSelection
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::PropertyList > ( ::SireMol::AtomProperty< SireBase::PropertyList >::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::matchToSelection );
            
            AtomPropertyList_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::merge
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireBase::PropertyList >::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::merge );
            
            AtomPropertyList_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::nAtoms
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::PropertyList >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::nAtoms );
            
            AtomPropertyList_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::nAtoms
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::PropertyList >::*nAtoms_function_type)( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::nAtoms );
            
            AtomPropertyList_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::nCutGroups
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::PropertyList >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::nCutGroups );
            
            AtomPropertyList_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomPropertyList_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireBase::PropertyList >::operator=
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::PropertyList > & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*assign_function_type)( ::SireMol::AtomProperty< SireBase::PropertyList > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::operator= );
            
            AtomPropertyList_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomPropertyList_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireBase::PropertyList >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireBase::PropertyList >::Array const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*__getitem___function_type)( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::operator[] );
            
            AtomPropertyList_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::operator[]
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireBase::PropertyList const & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*__getitem___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::operator[] );
            
            AtomPropertyList_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::set
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::PropertyList > & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireBase::PropertyList const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::set );
            
            AtomPropertyList_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::set
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< SireBase::PropertyList > & ( ::SireMol::AtomProperty< SireBase::PropertyList >::*set_function_type)( ::SireMol::CGIdx,::QVector< SireBase::PropertyList > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::set );
            
            AtomPropertyList_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::size
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireBase::PropertyList >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::size );
            
            AtomPropertyList_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::toString
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireBase::PropertyList >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::toString );
            
            AtomPropertyList_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::toVariant
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireBase::PropertyList >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::toVariant );
            
            AtomPropertyList_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::toVector
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::QVector< SireBase::PropertyList > ( ::SireMol::AtomProperty< SireBase::PropertyList >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::toVector );
            
            AtomPropertyList_exposer.def( 
                "toVector"
                , toVector_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::toVector
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef ::QVector< SireBase::PropertyList > ( ::SireMol::AtomProperty< SireBase::PropertyList >::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::toVector );
            
            AtomPropertyList_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireBase::PropertyList >::typeName
        
            typedef SireMol::AtomProperty< SireBase::PropertyList > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireBase::PropertyList >::typeName );
            
            AtomPropertyList_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        AtomPropertyList_exposer.staticmethod( "fromVariant" );
        AtomPropertyList_exposer.staticmethod( "typeName" );
        AtomPropertyList_exposer.def( "__copy__", &__copy__);
        AtomPropertyList_exposer.def( "__deepcopy__", &__copy__);
        AtomPropertyList_exposer.def( "clone", &__copy__);
        AtomPropertyList_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireBase::PropertyList> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPropertyList_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireBase::PropertyList> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomPropertyList_exposer.def_pickle(sire_pickle_suite< ::SireMol::AtomProperty<SireBase::PropertyList> >());
        AtomPropertyList_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireBase::PropertyList> > );
        AtomPropertyList_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireBase::PropertyList> > );
        AtomPropertyList_exposer.def( "__len__", &__len_size< ::SireMol::AtomProperty<SireBase::PropertyList> > );
    }

}
