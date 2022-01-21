// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ImproperID.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/property.h"

#include "SireMaths/torsion.h"

#include "SireMaths/vector.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atomcoords.h"

#include "improperid.h"

#include "moleculedata.h"

#include "moleculeinfodata.h"

#include "improperid.h"

SireMol::ImproperID __copy__(const SireMol::ImproperID &other){ return SireMol::ImproperID(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ImproperID_class(){

    { //::SireMol::ImproperID
        typedef bp::class_< SireMol::ImproperID, bp::bases< SireID::ID > > ImproperID_exposer_t;
        ImproperID_exposer_t ImproperID_exposer = ImproperID_exposer_t( "ImproperID", "This class provides a generic ID for an\nimproper angle between four atoms. The improper\nangle is the angle 1-2-3-4, which has the effect of measuring\nthe angle between the plane formed by the atoms 1,3,4 and the\nplane formed by the atoms 2,3,4. This measures by how much\natom 1 lies out of the plane formed by the atoms 2,3,4.\n\n1\n|\n2\n\\n3       4\n\nThe internal move Monte Carlo move changes an improper\nby splitting the molecule about the 1-2 bond, and then\nrotating the atom 1 group, and the atom 2-3-4 group about\nthe vector 3-4, about point 2.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ImproperID_scope( ImproperID_exposer );
        ImproperID_exposer.def( bp::init< SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const &, SireMol::AtomID const & >(( bp::arg("atom0"), bp::arg("atom1"), bp::arg("atom2"), bp::arg("atom3") ), "Construct a improper between the two specified atoms. The order\nis important, as this improper may be between two different\nmolecules") );
        ImproperID_exposer.def( bp::init< SireMol::ImproperID const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::ImproperID::atom0
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ImproperID::*atom0_function_type)(  ) const;
            atom0_function_type atom0_function_value( &::SireMol::ImproperID::atom0 );
            
            ImproperID_exposer.def( 
                "atom0"
                , atom0_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ID of the first atom of the improper" );
        
        }
        { //::SireMol::ImproperID::atom1
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ImproperID::*atom1_function_type)(  ) const;
            atom1_function_type atom1_function_value( &::SireMol::ImproperID::atom1 );
            
            ImproperID_exposer.def( 
                "atom1"
                , atom1_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ID of the second atom of the improper" );
        
        }
        { //::SireMol::ImproperID::atom2
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ImproperID::*atom2_function_type)(  ) const;
            atom2_function_type atom2_function_value( &::SireMol::ImproperID::atom2 );
            
            ImproperID_exposer.def( 
                "atom2"
                , atom2_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ID of the third atom of the improper" );
        
        }
        { //::SireMol::ImproperID::atom3
        
            typedef ::SireMol::AtomID const & ( ::SireMol::ImproperID::*atom3_function_type)(  ) const;
            atom3_function_type atom3_function_value( &::SireMol::ImproperID::atom3 );
            
            ImproperID_exposer.def( 
                "atom3"
                , atom3_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ID of the fourth atom of the improper" );
        
        }
        { //::SireMol::ImproperID::equivalent
        
            typedef bool ( ::SireMol::ImproperID::*equivalent_function_type)( ::SireMol::ImproperID const & ) const;
            equivalent_function_type equivalent_function_value( &::SireMol::ImproperID::equivalent );
            
            ImproperID_exposer.def( 
                "equivalent"
                , equivalent_function_value
                , ( bp::arg("other") )
                , "Are these impropers generally equivalent, i.e. do they contain the same\natom indices. This is useful since the ordering of improper atoms is\ninconsistent between different molecular topology formats.\n" );
        
        }
        { //::SireMol::ImproperID::hash
        
            typedef ::uint ( ::SireMol::ImproperID::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ImproperID::hash );
            
            ImproperID_exposer.def( 
                "hash"
                , hash_function_value
                , "Return a hash for this ID" );
        
        }
        { //::SireMol::ImproperID::isNull
        
            typedef bool ( ::SireMol::ImproperID::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ImproperID::isNull );
            
            ImproperID_exposer.def( 
                "isNull"
                , isNull_function_value
                , "Return whether this is a null ID" );
        
        }
        { //::SireMol::ImproperID::map
        
            typedef ::boost::tuples::tuple< SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMol::ImproperID::*map_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            map_function_type map_function_value( &::SireMol::ImproperID::map );
            
            ImproperID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "Return the indicies of the four atoms in this improper - this returns\nthem in the order\ntuple(improper.atom0(),improper.atom1(),improper.atom2(),improper.atom3())\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::map
        
            typedef ::boost::tuples::tuple< SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, SireMol::AtomIdx, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > ( ::SireMol::ImproperID::*map_function_type)( ::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const &,::SireMol::MoleculeInfoData const & ) const;
            map_function_type map_function_value( &::SireMol::ImproperID::map );
            
            ImproperID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("mol0info"), bp::arg("mol1info"), bp::arg("mol2info"), bp::arg("mol3info") )
                , "Return the indicies of the four atoms of this improper, between the\ntwo molecules whose data is in mol0info (containing improper.atom0()),\nmol1info (containing improper.atom1()), mol2info (containing\nimproper.atom2()) and mol3info (containing improper.atom3())\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        ImproperID_exposer.def( bp::self != bp::self );
        { //::SireMol::ImproperID::operator=
        
            typedef ::SireMol::ImproperID & ( ::SireMol::ImproperID::*assign_function_type)( ::SireMol::ImproperID const & ) ;
            assign_function_type assign_function_value( &::SireMol::ImproperID::operator= );
            
            ImproperID_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ImproperID_exposer.def( bp::self == bp::other< SireID::ID >() );
        ImproperID_exposer.def( bp::self == bp::self );
        { //::SireMol::ImproperID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::ImproperID::*size_function_type)( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::ImproperID::size );
            
            ImproperID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("moldata"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the size of this improper in the molecule whose data\nis in moldata, using map to find the coordinates property\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::ImproperID::*size_function_type)( ::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::ImproperID::size );
            
            ImproperID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("mol0data"), bp::arg("mol1data"), bp::arg("mol2data"), bp::arg("mol3data"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the size of the improper between atom0() in the\nmolecule whose data is in mol0data, atom1() in the\nmolecule whose data is in mol1data, atom2() in\nthe molecule whose data is in mol2data, and\natom3() in the molecule whose data is in mol3data,\nusing map to find the coordinates property of the\nmolecules\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::size
        
            typedef ::SireUnits::Dimension::Angle ( ::SireMol::ImproperID::*size_function_type)( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            size_function_type size_function_value( &::SireMol::ImproperID::size );
            
            ImproperID_exposer.def( 
                "size"
                , size_function_value
                , ( bp::arg("mol0data"), bp::arg("map0"), bp::arg("mol1data"), bp::arg("map1"), bp::arg("mol2data"), bp::arg("map2"), bp::arg("mol3data"), bp::arg("map3") )
                , "Return the size of the improper between atom0() in the\nmolecule whose data is in mol0data, atom1() in the\nmolecule whose data is in mol1data, atom2() in\nthe molecule whose data is in mol2data, and\natom3() in the molecule whose data is in mol3data, using map0\nto the find the coordinates property of mol0,\nmap1 to find the coordinates property of mol1,\nmap2 to find the coordinates property of mol2 and\nmap3 to find the coordinates property of mol3\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::toString
        
            typedef ::QString ( ::SireMol::ImproperID::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ImproperID::toString );
            
            ImproperID_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this ID" );
        
        }
        { //::SireMol::ImproperID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::ImproperID::*torsion_function_type)( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::ImproperID::torsion );
            
            ImproperID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("moldata"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the geometric torsion formed by the four atoms\nof this improper in the molecule whose data is in moldata,\nusing map to find the coordinates property.\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::ImproperID::*torsion_function_type)( ::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::ImproperID::torsion );
            
            ImproperID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("mol0data"), bp::arg("mol1data"), bp::arg("mol2data"), bp::arg("mol3data"), bp::arg("map")=SireBase::PropertyMap() )
                , "Return the geometric torsion formed by the four atoms,\natom0() in the molecule whose data is in mol0data,\natom1() from mol1data, atom2() from mol2data, and\natom3() from mol3data,\nusing map to find the coordinates property of\nthe molecules\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::torsion
        
            typedef ::SireMaths::Torsion ( ::SireMol::ImproperID::*torsion_function_type)( ::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const &,::SireMol::MoleculeData const &,::SireBase::PropertyMap const & ) const;
            torsion_function_type torsion_function_value( &::SireMol::ImproperID::torsion );
            
            ImproperID_exposer.def( 
                "torsion"
                , torsion_function_value
                , ( bp::arg("mol0data"), bp::arg("map0"), bp::arg("mol1data"), bp::arg("map1"), bp::arg("mol2data"), bp::arg("map2"), bp::arg("mol3data"), bp::arg("map3") )
                , "Return the geometric torsion formed by the four atoms,\natom0() in the molecule whose data is in mol0data,\natom1() from mol1data, atom2() from mol2data, and\natom3() from mol3data,\nusing map0 to find the coordinates property of mol0,\nmap1 to find the coordinates property of mol1,\nmap2 to find the coordinates property of mol2 and\nmap3 to find the coordinates property of mol3.\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::ImproperID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ImproperID::typeName );
            
            ImproperID_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ImproperID::what
        
            typedef char const * ( ::SireMol::ImproperID::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ImproperID::what );
            
            ImproperID_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ImproperID_exposer.staticmethod( "typeName" );
        ImproperID_exposer.def( "__copy__", &__copy__);
        ImproperID_exposer.def( "__deepcopy__", &__copy__);
        ImproperID_exposer.def( "clone", &__copy__);
        ImproperID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ImproperID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ImproperID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ImproperID >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ImproperID_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::ImproperID > );
        ImproperID_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::ImproperID > );
        ImproperID_exposer.def( "__str__", &__str__< ::SireMol::ImproperID > );
        ImproperID_exposer.def( "__repr__", &__str__< ::SireMol::ImproperID > );
        ImproperID_exposer.def( "__hash__", &::SireMol::ImproperID::hash );
    }

}
