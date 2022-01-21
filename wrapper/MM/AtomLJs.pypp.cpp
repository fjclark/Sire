// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomLJs.pypp.hpp"

namespace bp = boost::python;

#include "atomljs.h"

#include "atomljs.h"

#include "SireMol/moleculeview.h"

SireMol::AtomProperty<SireMM::LJParameter> __copy__(const SireMol::AtomProperty<SireMM::LJParameter> &other){ return SireMol::AtomProperty<SireMM::LJParameter>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_AtomLJs_class(){

    { //::SireMol::AtomProperty< SireMM::LJParameter >
        typedef bp::class_< SireMol::AtomProperty< SireMM::LJParameter >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomLJs_exposer_t;
        AtomLJs_exposer_t AtomLJs_exposer = AtomLJs_exposer_t( "AtomLJs", "", bp::init< >("") );
        bp::scope AtomLJs_scope( AtomLJs_exposer );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeInfo const & >(( bp::arg("molinfo") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeInfo const &, SireMM::LJParameter const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeView const &, SireMM::LJParameter const & >(( bp::arg("molview"), bp::arg("default_value") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireMM::LJParameter const & >(( bp::arg("molinfo"), bp::arg("default_value") ), "") );
        AtomLJs_exposer.def( bp::init< SireMM::LJParameter const & >(( bp::arg("value") ), "") );
        AtomLJs_exposer.def( bp::init< SireBase::PackedArray2D< SireMM::LJParameter > const & >(( bp::arg("values") ), "") );
        AtomLJs_exposer.def( bp::init< SireMol::AtomProperty< SireMM::LJParameter > const & >(( bp::arg("other") ), "") );
        { //::SireMol::AtomProperty< SireMM::LJParameter >::array
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMM::LJParameter > const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::array );
            
            AtomLJs_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMM::LJParameter >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::assertCanConvert );
            
            AtomLJs_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::assignFrom
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMM::LJParameter >::*assignFrom_function_type)( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::assignFrom );
            
            AtomLJs_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::at
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMM::LJParameter >::Array const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*at_function_type)( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::at );
            
            AtomLJs_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::at
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMM::LJParameter const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*at_function_type)( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::at );
            
            AtomLJs_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::canConvert
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMM::LJParameter >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::canConvert );
            
            AtomLJs_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::copyFrom
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMM::LJParameter >::*copyFrom_function_type)( ::QVector< SireMM::LJParameter > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::copyFrom );
            
            AtomLJs_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::copyFrom
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireMM::LJParameter >::*copyFrom_function_type)( ::QVector< SireMM::LJParameter > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::copyFrom );
            
            AtomLJs_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::count
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMM::LJParameter >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::count );
            
            AtomLJs_exposer.def( 
                "count"
                , count_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::divide
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMM::LJParameter >::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::divide );
            
            AtomLJs_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::divideByResidue
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMM::LJParameter >::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::divideByResidue );
            
            AtomLJs_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::fromVariant
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMM::LJParameter > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::fromVariant );
            
            AtomLJs_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::get
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMM::LJParameter >::Array const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::get );
            
            AtomLJs_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::get
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMM::LJParameter const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::get );
            
            AtomLJs_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMM::LJParameter >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::isCompatibleWith );
            
            AtomLJs_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMM::LJParameter >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfo const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::isCompatibleWith );
            
            AtomLJs_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::isEmpty
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireMM::LJParameter >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::isEmpty );
            
            AtomLJs_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::matchToSelection
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMM::LJParameter > ( ::SireMol::AtomProperty< SireMM::LJParameter >::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::matchToSelection );
            
            AtomLJs_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::merge
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::AtomProperty< SireMM::LJParameter >::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::merge );
            
            AtomLJs_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::nAtoms
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMM::LJParameter >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::nAtoms );
            
            AtomLJs_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::nAtoms
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMM::LJParameter >::*nAtoms_function_type)( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::nAtoms );
            
            AtomLJs_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::nCutGroups
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMM::LJParameter >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::nCutGroups );
            
            AtomLJs_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , "" );
        
        }
        AtomLJs_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireMM::LJParameter >::operator=
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMM::LJParameter > & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*assign_function_type)( ::SireMol::AtomProperty< SireMM::LJParameter > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::operator= );
            
            AtomLJs_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AtomLJs_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireMM::LJParameter >::operator[]
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireMM::LJParameter >::Array const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*__getitem___function_type)( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::operator[] );
            
            AtomLJs_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::operator[]
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMM::LJParameter const & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*__getitem___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::operator[] );
            
            AtomLJs_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::set
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMM::LJParameter > & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireMM::LJParameter const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::set );
            
            AtomLJs_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::set
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< SireMM::LJParameter > & ( ::SireMol::AtomProperty< SireMM::LJParameter >::*set_function_type)( ::SireMol::CGIdx,::QVector< SireMM::LJParameter > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::set );
            
            AtomLJs_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::size
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireMM::LJParameter >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::size );
            
            AtomLJs_exposer.def( 
                "size"
                , size_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::toString
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireMM::LJParameter >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::toString );
            
            AtomLJs_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::toVariant
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireMM::LJParameter >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::toVariant );
            
            AtomLJs_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::toVector
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::QVector< SireMM::LJParameter > ( ::SireMol::AtomProperty< SireMM::LJParameter >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::toVector );
            
            AtomLJs_exposer.def( 
                "toVector"
                , toVector_function_value
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::toVector
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef ::QVector< SireMM::LJParameter > ( ::SireMol::AtomProperty< SireMM::LJParameter >::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::toVector );
            
            AtomLJs_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") )
                , "" );
        
        }
        { //::SireMol::AtomProperty< SireMM::LJParameter >::typeName
        
            typedef SireMol::AtomProperty< SireMM::LJParameter > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireMM::LJParameter >::typeName );
            
            AtomLJs_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        AtomLJs_exposer.staticmethod( "fromVariant" );
        AtomLJs_exposer.staticmethod( "typeName" );
        AtomLJs_exposer.def( "__copy__", &__copy__);
        AtomLJs_exposer.def( "__deepcopy__", &__copy__);
        AtomLJs_exposer.def( "clone", &__copy__);
        AtomLJs_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireMM::LJParameter> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomLJs_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireMM::LJParameter> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomLJs_exposer.def( "__getstate_manages_dict__", true);
        AtomLJs_exposer.def( "__safe_for_unpickling__", true);
        AtomLJs_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::AtomProperty<SireMM::LJParameter> > );
        AtomLJs_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::AtomProperty<SireMM::LJParameter> > );
        AtomLJs_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireMM::LJParameter> > );
        AtomLJs_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireMM::LJParameter> > );
        AtomLJs_exposer.def( "__len__", &__len_size< ::SireMol::AtomProperty<SireMM::LJParameter> > );
    }

}
