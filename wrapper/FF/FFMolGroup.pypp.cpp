// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "FFMolGroup.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/propertymap.h"

#include "SireError/errors.h"

#include "SireMol/mgname.h"

#include "SireMol/mgnum.h"

#include "SireMol/molecule.h"

#include "SireMol/molecules.h"

#include "SireMol/molnum.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireMol/viewsofmol.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "ffmolgroup.h"

#include <QDebug>

#include "ffmolgroup.h"

SireFF::FFMolGroup __copy__(const SireFF::FFMolGroup &other){ return SireFF::FFMolGroup(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_FFMolGroup_class(){

    { //::SireFF::FFMolGroup
        typedef bp::class_< SireFF::FFMolGroup, bp::bases< SireMol::MoleculeGroup, SireBase::Property > > FFMolGroup_exposer_t;
        FFMolGroup_exposer_t FFMolGroup_exposer = FFMolGroup_exposer_t( "FFMolGroup", "This class holds a molecule group that is part of a forcefield.\nThis actually holds a copy of the forcefield that contains\nthis molecule group, so that changes to this group change\nthe actual forcefield to which this group belongs. This is\nthe publicly available class that corresponds to the private\nSireFF::detail::FFMolGroupPvt class.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor - this creates a useless FFMolGroup. You can\nonly construct a valid FFMolGroup by getting a reference to\na valid FFMolGroupPvt from a forcefield") );
        bp::scope FFMolGroup_scope( FFMolGroup_exposer );
        FFMolGroup_exposer.def( bp::init< SireFF::detail::FFMolGroupPvt const & >(( bp::arg("ffmolgroup") ), "Construct from an FFMolGroupPvt - this grabs a copy of the\nforcefield that contains the FFMolGroupPvt") );
        FFMolGroup_exposer.def( bp::init< SireMol::MoleculeGroup const & >(( bp::arg("other") ), "Construct from a MoleculeGroup - you can only construct from a FFMolGroup\nor from an FFMolGroupPvt\nThrow: SireError::invalid_arg\n") );
        FFMolGroup_exposer.def( bp::init< SireFF::FFMolGroup const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireFF::FFMolGroup::accept
        
            typedef void ( ::SireFF::FFMolGroup::*accept_function_type)(  ) ;
            accept_function_type accept_function_value( &::SireFF::FFMolGroup::accept );
            
            FFMolGroup_exposer.def( 
                "accept"
                , accept_function_value
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::MoleculeView const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molview") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::ViewsOfMol const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molviews") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::Molecules const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::MoleculeGroup const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molgroup") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molview"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::ViewsOfMol const &,::SireBase::PropertyMap const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molviews"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molecules"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::add
        
            typedef void ( ::SireFF::FFMolGroup::*add_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            add_function_type add_function_value( &::SireFF::FFMolGroup::add );
            
            FFMolGroup_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molgroup"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef bool ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::MoleculeView const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molview") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::SireMol::ViewsOfMol ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::ViewsOfMol const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molviews") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::Molecules const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::MoleculeGroup const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molgroup") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef bool ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molview"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::SireMol::ViewsOfMol ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::ViewsOfMol const &,::SireBase::PropertyMap const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molviews"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molecules"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::addIfUnique
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*addIfUnique_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireFF::FFMolGroup::addIfUnique );
            
            FFMolGroup_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molgroup"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::forceField
        
            typedef ::SireFF::FF const & ( ::SireFF::FFMolGroup::*forceField_function_type)(  ) const;
            forceField_function_type forceField_function_value( &::SireFF::FFMolGroup::forceField );
            
            FFMolGroup_exposer.def( 
                "forceField"
                , forceField_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the forcefield that contains this molecule group" );
        
        }
        { //::SireFF::FFMolGroup::index
        
            typedef ::SireMol::MGIdx ( ::SireFF::FFMolGroup::*index_function_type)(  ) const;
            index_function_type index_function_value( &::SireFF::FFMolGroup::index );
            
            FFMolGroup_exposer.def( 
                "index"
                , index_function_value
                , "Return the index of this group in the parent forcefield" );
        
        }
        { //::SireFF::FFMolGroup::needsAccepting
        
            typedef bool ( ::SireFF::FFMolGroup::*needsAccepting_function_type)(  ) const;
            needsAccepting_function_type needsAccepting_function_value( &::SireFF::FFMolGroup::needsAccepting );
            
            FFMolGroup_exposer.def( 
                "needsAccepting"
                , needsAccepting_function_value
                , "" );
        
        }
        { //::SireFF::FFMolGroup::operator=
        
            typedef ::SireFF::FFMolGroup & ( ::SireFF::FFMolGroup::*assign_function_type)( ::SireFF::FFMolGroup const & ) ;
            assign_function_type assign_function_value( &::SireFF::FFMolGroup::operator= );
            
            FFMolGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::FFMolGroup::operator=
        
            typedef ::SireFF::FFMolGroup & ( ::SireFF::FFMolGroup::*assign_function_type)( ::SireMol::MoleculeGroup const & ) ;
            assign_function_type assign_function_value( &::SireFF::FFMolGroup::operator= );
            
            FFMolGroup_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef bool ( ::SireFF::FFMolGroup::*remove_function_type)( ::SireMol::MoleculeView const & ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molview") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef ::SireMol::ViewsOfMol ( ::SireFF::FFMolGroup::*remove_function_type)( ::SireMol::ViewsOfMol const & ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molviews") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*remove_function_type)( ::SireMol::Molecules const & ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*remove_function_type)( ::SireMol::MoleculeGroup const & ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molgroup") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef ::SireMol::ViewsOfMol ( ::SireFF::FFMolGroup::*remove_function_type)( ::SireMol::MolNum ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molnum") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::remove
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*remove_function_type)( ::QSet< SireMol::MolNum > const & ) ;
            remove_function_type remove_function_value( &::SireFF::FFMolGroup::remove );
            
            FFMolGroup_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molnums") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::removeAll
        
            typedef bool ( ::SireFF::FFMolGroup::*removeAll_function_type)( ::SireMol::MoleculeView const & ) ;
            removeAll_function_type removeAll_function_value( &::SireFF::FFMolGroup::removeAll );
            
            FFMolGroup_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molview") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::removeAll
        
            typedef ::SireMol::ViewsOfMol ( ::SireFF::FFMolGroup::*removeAll_function_type)( ::SireMol::ViewsOfMol const & ) ;
            removeAll_function_type removeAll_function_value( &::SireFF::FFMolGroup::removeAll );
            
            FFMolGroup_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molviews") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::removeAll
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*removeAll_function_type)( ::SireMol::Molecules const & ) ;
            removeAll_function_type removeAll_function_value( &::SireFF::FFMolGroup::removeAll );
            
            FFMolGroup_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::removeAll
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireFF::FFMolGroup::*removeAll_function_type)( ::SireMol::MoleculeGroup const & ) ;
            removeAll_function_type removeAll_function_value( &::SireFF::FFMolGroup::removeAll );
            
            FFMolGroup_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molgroup") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::removeAll
        
            typedef void ( ::SireFF::FFMolGroup::*removeAll_function_type)(  ) ;
            removeAll_function_type removeAll_function_value( &::SireFF::FFMolGroup::removeAll );
            
            FFMolGroup_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::MoleculeView const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molview") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::ViewsOfMol const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molviews") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::Molecules const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molecules") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::MoleculeGroup const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molgroup") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molview"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::ViewsOfMol const &,::SireBase::PropertyMap const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molviews"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::Molecules const &,::SireBase::PropertyMap const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molecules"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setContents
        
            typedef bool ( ::SireFF::FFMolGroup::*setContents_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            setContents_function_type setContents_function_value( &::SireFF::FFMolGroup::setContents );
            
            FFMolGroup_exposer.def( 
                "setContents"
                , setContents_function_value
                , ( bp::arg("molgroup"), bp::arg("map") )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::setName
        
            typedef void ( ::SireFF::FFMolGroup::*setName_function_type)( ::QString const & ) ;
            setName_function_type setName_function_value( &::SireFF::FFMolGroup::setName );
            
            FFMolGroup_exposer.def( 
                "setName"
                , setName_function_value
                , ( bp::arg("name") )
                , "Set the name of this molecule group" );
        
        }
        { //::SireFF::FFMolGroup::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::FFMolGroup::typeName );
            
            FFMolGroup_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::FFMolGroup::update
        
            typedef bool ( ::SireFF::FFMolGroup::*update_function_type)( ::SireMol::MoleculeData const &,bool ) ;
            update_function_type update_function_value( &::SireFF::FFMolGroup::update );
            
            FFMolGroup_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata"), bp::arg("auto_commit")=(bool)(true) )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::update
        
            typedef ::QList< SireMol::Molecule > ( ::SireFF::FFMolGroup::*update_function_type)( ::SireMol::Molecules const &,bool ) ;
            update_function_type update_function_value( &::SireFF::FFMolGroup::update );
            
            FFMolGroup_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules"), bp::arg("auto_commit")=(bool)(true) )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::update
        
            typedef ::QList< SireMol::Molecule > ( ::SireFF::FFMolGroup::*update_function_type)( ::SireMol::MoleculeGroup const &,bool ) ;
            update_function_type update_function_value( &::SireFF::FFMolGroup::update );
            
            FFMolGroup_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup"), bp::arg("auto_commit")=(bool)(true) )
                , "" );
        
        }
        { //::SireFF::FFMolGroup::what
        
            typedef char const * ( ::SireFF::FFMolGroup::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::FFMolGroup::what );
            
            FFMolGroup_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        FFMolGroup_exposer.staticmethod( "typeName" );
        FFMolGroup_exposer.def( "__copy__", &__copy__);
        FFMolGroup_exposer.def( "__deepcopy__", &__copy__);
        FFMolGroup_exposer.def( "clone", &__copy__);
        FFMolGroup_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::FFMolGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFMolGroup_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::FFMolGroup >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFMolGroup_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::FFMolGroup > );
        FFMolGroup_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::FFMolGroup > );
        FFMolGroup_exposer.def( "__str__", &__str__< ::SireFF::FFMolGroup > );
        FFMolGroup_exposer.def( "__repr__", &__str__< ::SireFF::FFMolGroup > );
    }

}
