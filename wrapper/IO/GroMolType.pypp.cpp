// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "GroMolType.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/booleanproperty.h"

#include "SireBase/parallel.h"

#include "SireBase/stringproperty.h"

#include "SireError/errors.h"

#include "SireIO/errors.h"

#include "SireIO/grotop.h"

#include "SireMM/atomljs.h"

#include "SireMM/cljnbpairs.h"

#include "SireMM/fouratomfunctions.h"

#include "SireMM/internalff.h"

#include "SireMM/threeatomfunctions.h"

#include "SireMM/twoatomfunctions.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/connectivity.h"

#include "SireMol/errors.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleditor.h"

#include "SireMol/select.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "grotop.h"

#include <QDateTime>

#include <QDir>

#include <QElapsedTimer>

#include <QFileInfo>

#include <QRegularExpression>

#include "grotop.h"

SireIO::GroMolType __copy__(const SireIO::GroMolType &other){ return SireIO::GroMolType(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_GroMolType_class(){

    { //::SireIO::GroMolType
        typedef bp::class_< SireIO::GroMolType > GroMolType_exposer_t;
        GroMolType_exposer_t GroMolType_exposer = GroMolType_exposer_t( "GroMolType", "This class is used by GroTop to hold an intermediate representation of a\nGromacs moleculetype. This provides metadata about the molecule that is\nneeded to construct the whole molecule.\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope GroMolType_scope( GroMolType_exposer );
        GroMolType_exposer.def( bp::init< SireMol::Molecule const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("mol"), bp::arg("map")=SireBase::PropertyMap() ), "Construct from the passed molecule") );
        GroMolType_exposer.def( bp::init< SireIO::GroMolType const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::GroMolType::addAngle
        
            typedef void ( ::SireIO::GroMolType::*addAngle_function_type)( ::SireMol::AngleID const &,::SireMM::GromacsAngle const &,bool ) ;
            addAngle_function_type addAngle_function_value( &::SireIO::GroMolType::addAngle );
            
            GroMolType_exposer.def( 
                "addAngle"
                , addAngle_function_value
                , ( bp::arg("angle"), bp::arg("parm"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed angle to the molecule" );
        
        }
        { //::SireIO::GroMolType::addAngles
        
            typedef void ( ::SireIO::GroMolType::*addAngles_function_type)( ::QMultiHash< SireMol::AngleID, SireMM::GromacsAngle > const &,bool ) ;
            addAngles_function_type addAngles_function_value( &::SireIO::GroMolType::addAngles );
            
            GroMolType_exposer.def( 
                "addAngles"
                , addAngles_function_value
                , ( bp::arg("angles"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed angles to the molecule" );
        
        }
        { //::SireIO::GroMolType::addAtom
        
            typedef void ( ::SireIO::GroMolType::*addAtom_function_type)( ::SireIO::GroAtom const &,bool ) ;
            addAtom_function_type addAtom_function_value( &::SireIO::GroMolType::addAtom );
            
            GroMolType_exposer.def( 
                "addAtom"
                , addAtom_function_value
                , ( bp::arg("atom"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add an atom to this moleculetype, with specified atom type, residue number,\nresidue name, atom name, charge group, charge and mass" );
        
        }
        { //::SireIO::GroMolType::addBond
        
            typedef void ( ::SireIO::GroMolType::*addBond_function_type)( ::SireMol::BondID const &,::SireMM::GromacsBond const &,bool ) ;
            addBond_function_type addBond_function_value( &::SireIO::GroMolType::addBond );
            
            GroMolType_exposer.def( 
                "addBond"
                , addBond_function_value
                , ( bp::arg("bond"), bp::arg("parm"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed bond to the molecule" );
        
        }
        { //::SireIO::GroMolType::addBonds
        
            typedef void ( ::SireIO::GroMolType::*addBonds_function_type)( ::QMultiHash< SireMol::BondID, SireMM::GromacsBond > const &,bool ) ;
            addBonds_function_type addBonds_function_value( &::SireIO::GroMolType::addBonds );
            
            GroMolType_exposer.def( 
                "addBonds"
                , addBonds_function_value
                , ( bp::arg("bonds"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed bonds to the molecule" );
        
        }
        { //::SireIO::GroMolType::addDihedral
        
            typedef void ( ::SireIO::GroMolType::*addDihedral_function_type)( ::SireMol::DihedralID const &,::SireMM::GromacsDihedral const &,bool ) ;
            addDihedral_function_type addDihedral_function_value( &::SireIO::GroMolType::addDihedral );
            
            GroMolType_exposer.def( 
                "addDihedral"
                , addDihedral_function_value
                , ( bp::arg("dihedral"), bp::arg("parm"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed dihedral to the molecule" );
        
        }
        { //::SireIO::GroMolType::addDihedrals
        
            typedef void ( ::SireIO::GroMolType::*addDihedrals_function_type)( ::QMultiHash< SireMol::DihedralID, SireMM::GromacsDihedral > const &,bool ) ;
            addDihedrals_function_type addDihedrals_function_value( &::SireIO::GroMolType::addDihedrals );
            
            GroMolType_exposer.def( 
                "addDihedrals"
                , addDihedrals_function_value
                , ( bp::arg("dihedrals"), bp::arg("is_lambda1")=(bool)(false) )
                , "Add the passed dihedrals to the molecule" );
        
        }
        { //::SireIO::GroMolType::addWarning
        
            typedef void ( ::SireIO::GroMolType::*addWarning_function_type)( ::QString const & ) ;
            addWarning_function_type addWarning_function_value( &::SireIO::GroMolType::addWarning );
            
            GroMolType_exposer.def( 
                "addWarning"
                , addWarning_function_value
                , ( bp::arg("warning") )
                , bp::release_gil_policy()
                , "Add a warning that has been generated while parsing or creatig this object" );
        
        }
        { //::SireIO::GroMolType::angles
        
            typedef ::QMultiHash< SireMol::AngleID, SireMM::GromacsAngle > ( ::SireIO::GroMolType::*angles_function_type)( bool ) const;
            angles_function_type angles_function_value( &::SireIO::GroMolType::angles );
            
            GroMolType_exposer.def( 
                "angles"
                , angles_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the angles" );
        
        }
        { //::SireIO::GroMolType::atom
        
            typedef ::SireIO::GroAtom ( ::SireIO::GroMolType::*atom_function_type)( ::SireMol::AtomIdx const &,bool ) const;
            atom_function_type atom_function_value( &::SireIO::GroMolType::atom );
            
            GroMolType_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomidx"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return the atom at index atomidx" );
        
        }
        { //::SireIO::GroMolType::atom
        
            typedef ::SireIO::GroAtom ( ::SireIO::GroMolType::*atom_function_type)( ::SireMol::AtomNum const &,bool ) const;
            atom_function_type atom_function_value( &::SireIO::GroMolType::atom );
            
            GroMolType_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomnum"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return the atom with number atomnum" );
        
        }
        { //::SireIO::GroMolType::atom
        
            typedef ::SireIO::GroAtom ( ::SireIO::GroMolType::*atom_function_type)( ::SireMol::AtomName const &,bool ) const;
            atom_function_type atom_function_value( &::SireIO::GroMolType::atom );
            
            GroMolType_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomnam"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return the first atom with name atomnam. If you want all atoms\nwith this name then call atoms(AtomName atomname)" );
        
        }
        { //::SireIO::GroMolType::atoms
        
            typedef ::QVector< SireIO::GroAtom > ( ::SireIO::GroMolType::*atoms_function_type)( bool ) const;
            atoms_function_type atoms_function_value( &::SireIO::GroMolType::atoms );
            
            GroMolType_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the atoms in this molecule" );
        
        }
        { //::SireIO::GroMolType::atoms
        
            typedef ::QVector< SireIO::GroAtom > ( ::SireIO::GroMolType::*atoms_function_type)( ::SireMol::AtomName const &,bool ) const;
            atoms_function_type atoms_function_value( &::SireIO::GroMolType::atoms );
            
            GroMolType_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("atomnam"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return all atoms that have the passed name. Returns an empty\nlist if there are no atoms with this name" );
        
        }
        { //::SireIO::GroMolType::atoms
        
            typedef ::QVector< SireIO::GroAtom > ( ::SireIO::GroMolType::*atoms_function_type)( ::SireMol::ResIdx const &,bool ) const;
            atoms_function_type atoms_function_value( &::SireIO::GroMolType::atoms );
            
            GroMolType_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("residx"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the atoms in the specified residue" );
        
        }
        { //::SireIO::GroMolType::atoms
        
            typedef ::QVector< SireIO::GroAtom > ( ::SireIO::GroMolType::*atoms_function_type)( ::SireMol::ResNum const &,bool ) const;
            atoms_function_type atoms_function_value( &::SireIO::GroMolType::atoms );
            
            GroMolType_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("resnum"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the atoms in the specified residue(s)" );
        
        }
        { //::SireIO::GroMolType::atoms
        
            typedef ::QVector< SireIO::GroAtom > ( ::SireIO::GroMolType::*atoms_function_type)( ::SireMol::ResName const &,bool ) const;
            atoms_function_type atoms_function_value( &::SireIO::GroMolType::atoms );
            
            GroMolType_exposer.def( 
                "atoms"
                , atoms_function_value
                , ( bp::arg("resnam"), bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the atoms in the specified residue(s)" );
        
        }
        { //::SireIO::GroMolType::bonds
        
            typedef ::QMultiHash< SireMol::BondID, SireMM::GromacsBond > ( ::SireIO::GroMolType::*bonds_function_type)( bool ) const;
            bonds_function_type bonds_function_value( &::SireIO::GroMolType::bonds );
            
            GroMolType_exposer.def( 
                "bonds"
                , bonds_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the bonds" );
        
        }
        { //::SireIO::GroMolType::dihedrals
        
            typedef ::QMultiHash< SireMol::DihedralID, SireMM::GromacsDihedral > ( ::SireIO::GroMolType::*dihedrals_function_type)( bool ) const;
            dihedrals_function_type dihedrals_function_value( &::SireIO::GroMolType::dihedrals );
            
            GroMolType_exposer.def( 
                "dihedrals"
                , dihedrals_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return all of the dihedrals" );
        
        }
        { //::SireIO::GroMolType::forcefield
        
            typedef ::SireMM::MMDetail ( ::SireIO::GroMolType::*forcefield_function_type)( bool ) const;
            forcefield_function_type forcefield_function_value( &::SireIO::GroMolType::forcefield );
            
            GroMolType_exposer.def( 
                "forcefield"
                , forcefield_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return the guessed forcefield for this molecule type" );
        
        }
        { //::SireIO::GroMolType::isNull
        
            typedef bool ( ::SireIO::GroMolType::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireIO::GroMolType::isNull );
            
            GroMolType_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "Return whether or not this object is null" );
        
        }
        { //::SireIO::GroMolType::isPerturbable
        
            typedef bool ( ::SireIO::GroMolType::*isPerturbable_function_type)(  ) const;
            isPerturbable_function_type isPerturbable_function_value( &::SireIO::GroMolType::isPerturbable );
            
            GroMolType_exposer.def( 
                "isPerturbable"
                , isPerturbable_function_value
                , bp::release_gil_policy()
                , "Return whether or not this molecule is perturbable." );
        
        }
        { //::SireIO::GroMolType::isWater
        
            typedef bool ( ::SireIO::GroMolType::*isWater_function_type)( bool ) const;
            isWater_function_type isWater_function_value( &::SireIO::GroMolType::isWater );
            
            GroMolType_exposer.def( 
                "isWater"
                , isWater_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return whether or not this is a topology for water. This should\nreturn true for all water models (including TIP4P and TIP5P)" );
        
        }
        { //::SireIO::GroMolType::nAtoms
        
            typedef int ( ::SireIO::GroMolType::*nAtoms_function_type)( bool ) const;
            nAtoms_function_type nAtoms_function_value( &::SireIO::GroMolType::nAtoms );
            
            GroMolType_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return the number of atoms in this molecule" );
        
        }
        { //::SireIO::GroMolType::nExcludedAtoms
        
            typedef ::qint64 ( ::SireIO::GroMolType::*nExcludedAtoms_function_type)( bool ) const;
            nExcludedAtoms_function_type nExcludedAtoms_function_value( &::SireIO::GroMolType::nExcludedAtoms );
            
            GroMolType_exposer.def( 
                "nExcludedAtoms"
                , nExcludedAtoms_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return the number of excluded atoms" );
        
        }
        { //::SireIO::GroMolType::nResidues
        
            typedef int ( ::SireIO::GroMolType::*nResidues_function_type)( bool ) const;
            nResidues_function_type nResidues_function_value( &::SireIO::GroMolType::nResidues );
            
            GroMolType_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return the number of residues in this molecule" );
        
        }
        { //::SireIO::GroMolType::name
        
            typedef ::QString ( ::SireIO::GroMolType::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireIO::GroMolType::name );
            
            GroMolType_exposer.def( 
                "name"
                , name_function_value
                , bp::release_gil_policy()
                , "Return the name of this moleculetype" );
        
        }
        { //::SireIO::GroMolType::needsSanitising
        
            typedef bool ( ::SireIO::GroMolType::*needsSanitising_function_type)( bool ) const;
            needsSanitising_function_type needsSanitising_function_value( &::SireIO::GroMolType::needsSanitising );
            
            GroMolType_exposer.def( 
                "needsSanitising"
                , needsSanitising_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return whether or not this molecule needs sanitising" );
        
        }
        GroMolType_exposer.def( bp::self != bp::self );
        { //::SireIO::GroMolType::operator=
        
            typedef ::SireIO::GroMolType & ( ::SireIO::GroMolType::*assign_function_type)( ::SireIO::GroMolType const & ) ;
            assign_function_type assign_function_value( &::SireIO::GroMolType::operator= );
            
            GroMolType_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        GroMolType_exposer.def( bp::self == bp::self );
        { //::SireIO::GroMolType::sanitise
        
            typedef void ( ::SireIO::GroMolType::*sanitise_function_type)( ::QString,::QString,::QString,double,double,bool ) ;
            sanitise_function_type sanitise_function_value( &::SireIO::GroMolType::sanitise );
            
            GroMolType_exposer.def( 
                "sanitise"
                , sanitise_function_value
                , ( bp::arg("elecstyle"), bp::arg("vdwstyle"), bp::arg("combrule"), bp::arg("elec14"), bp::arg("vdw14"), bp::arg("is_lambda1")=(bool)(false) )
                , "Sanitise this moleculetype. This assumes that the moleculetype has\nbeen fully specified, so it collects everything together and checks that the\nmolecule makes sense. Any warnings generated can be retrieved using the\nwarnings function. It also uses the passed defaults from the top file,\ntogether with the information in the molecule to guess the forcefield for\nthe molecule" );
        
        }
        { //::SireIO::GroMolType::setNExcludedAtoms
        
            typedef void ( ::SireIO::GroMolType::*setNExcludedAtoms_function_type)( ::qint64,bool ) ;
            setNExcludedAtoms_function_type setNExcludedAtoms_function_value( &::SireIO::GroMolType::setNExcludedAtoms );
            
            GroMolType_exposer.def( 
                "setNExcludedAtoms"
                , setNExcludedAtoms_function_value
                , ( bp::arg("nexcl"), bp::arg("is_lambda1")=(bool)(false) )
                , "Set the number of excluded atoms" );
        
        }
        { //::SireIO::GroMolType::setName
        
            typedef void ( ::SireIO::GroMolType::*setName_function_type)( ::QString const & ) ;
            setName_function_type setName_function_value( &::SireIO::GroMolType::setName );
            
            GroMolType_exposer.def( 
                "setName"
                , setName_function_value
                , ( bp::arg("name") )
                , bp::release_gil_policy()
                , "Set the name of this moleculetype" );
        
        }
        { //::SireIO::GroMolType::settlesLines
        
            typedef ::QStringList ( ::SireIO::GroMolType::*settlesLines_function_type)( bool ) const;
            settlesLines_function_type settlesLines_function_value( &::SireIO::GroMolType::settlesLines );
            
            GroMolType_exposer.def( 
                "settlesLines"
                , settlesLines_function_value
                , ( bp::arg("is_lambda1")=(bool)(false) )
                , "Return the settles lines for this molecule. This currently only returns\nsettles lines for water molecules. These lines are used to constrain the\nbondsangles of the water molecule" );
        
        }
        { //::SireIO::GroMolType::toString
        
            typedef ::QString ( ::SireIO::GroMolType::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireIO::GroMolType::toString );
            
            GroMolType_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string form for this object" );
        
        }
        { //::SireIO::GroMolType::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::GroMolType::typeName );
            
            GroMolType_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::GroMolType::warnings
        
            typedef ::QStringList ( ::SireIO::GroMolType::*warnings_function_type)(  ) const;
            warnings_function_type warnings_function_value( &::SireIO::GroMolType::warnings );
            
            GroMolType_exposer.def( 
                "warnings"
                , warnings_function_value
                , bp::release_gil_policy()
                , "Return any warnings associated with this moleculetype" );
        
        }
        { //::SireIO::GroMolType::what
        
            typedef char const * ( ::SireIO::GroMolType::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireIO::GroMolType::what );
            
            GroMolType_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        GroMolType_exposer.staticmethod( "typeName" );
        GroMolType_exposer.def( "__copy__", &__copy__);
        GroMolType_exposer.def( "__deepcopy__", &__copy__);
        GroMolType_exposer.def( "clone", &__copy__);
        GroMolType_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::GroMolType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroMolType_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::GroMolType >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        GroMolType_exposer.def_pickle(sire_pickle_suite< ::SireIO::GroMolType >());
        GroMolType_exposer.def( "__str__", &__str__< ::SireIO::GroMolType > );
        GroMolType_exposer.def( "__repr__", &__str__< ::SireIO::GroMolType > );
    }

}
