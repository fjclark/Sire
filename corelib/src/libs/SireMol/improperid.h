/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef SIREMOL_IMPROPERID_H
#define SIREMOL_IMPROPERID_H

#include "atomidentifier.h"

#include "SireBase/propertymap.h"
#include "SireUnits/dimensions.h"

#include <boost/tuple/tuple.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class ImproperID;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::ImproperID&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::ImproperID&);

namespace SireMaths
{
class Vector;
class Torsion;
}

namespace SireMol
{

class MoleculeData;
class MoleculeInfoData;
class AtomIdx;

using SireMaths::Vector;
using SireMaths::Torsion;

using SireBase::PropertyMap;

using boost::tuple;

/** This class provides a generic ID for an
    improper angle between four atoms. The improper
    angle is the angle 1-2-3-4, which has the effect of measuring
    the angle between the plane formed by the atoms 1,3,4 and the
    plane formed by the atoms 2,3,4. This measures by how much
    atom 1 lies out of the plane formed by the atoms 2,3,4.

            1
            |
            2
          /   \
        3       4

    The internal move Monte Carlo move changes an improper
    by splitting the molecule about the 1-2 bond, and then
    rotating the atom 1 group, and the atom 2-3-4 group about
    the vector 3-4, about point 2.

    @author Christopher Woods
*/
class SIREMOL_EXPORT ImproperID : public SireID::ID
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const ImproperID&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, ImproperID&);

public:
    ImproperID();
    ImproperID(const AtomID &atom0, const AtomID &atom1,
               const AtomID &atom2, const AtomID &atom3);

    ImproperID(const ImproperID &other);

    ~ImproperID();

    static const char* typeName();

    const char* what() const
    {
        return ImproperID::typeName();
    }

    ImproperID* clone() const;

    uint hash() const;

    QString toString() const;

    bool isNull() const;

    ImproperID& operator=(const ImproperID &other);

    bool operator==(const SireID::ID &other) const;

    bool operator==(const ImproperID &other) const;
    bool operator!=(const ImproperID &other) const;

    bool equivalent(const ImproperID &other) const;

    const AtomID& operator[](int i) const;

    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>
    map(const MoleculeInfoData &molinfo) const;

    tuple<AtomIdx,AtomIdx,AtomIdx,AtomIdx>
    map(const MoleculeInfoData &mol0info,
        const MoleculeInfoData &mol1info,
        const MoleculeInfoData &mol2info,
        const MoleculeInfoData &mol3info) const;

    Torsion torsion(const MoleculeData &moldata,
                    const PropertyMap &map = PropertyMap()) const;

    Torsion torsion(const MoleculeData &mol0data,
                    const MoleculeData &mol1data,
                    const MoleculeData &mol2data,
                    const MoleculeData &mol3data,
                    const PropertyMap &map = PropertyMap()) const;

    Torsion torsion(const MoleculeData &mol0data,
                    const PropertyMap &map0,
                    const MoleculeData &mol1data,
                    const PropertyMap &map1,
                    const MoleculeData &mol2data,
                    const PropertyMap &map2,
                    const MoleculeData &mol3data,
                    const PropertyMap &map3) const;

    SireUnits::Dimension::Angle size(const MoleculeData &moldata,
                                     const PropertyMap &map = PropertyMap()) const;

    SireUnits::Dimension::Angle size(const MoleculeData &mol0data,
                                     const MoleculeData &mol1data,
                                     const MoleculeData &mol2data,
                                     const MoleculeData &mol3data,
                                     const PropertyMap &map = PropertyMap()) const;

    SireUnits::Dimension::Angle size(const MoleculeData &mol0data,
                                     const PropertyMap &map0,
                                     const MoleculeData &mol1data,
                                     const PropertyMap &map1,
                                     const MoleculeData &mol2data,
                                     const PropertyMap &map2,
                                     const MoleculeData &mol3data,
                                     const PropertyMap &map3) const;

    const AtomID& atom0() const;
    const AtomID& atom1() const;
    const AtomID& atom2() const;
    const AtomID& atom3() const;

private:
    /** The identifiers of the four atoms */
    AtomIdentifier atm0,atm1,atm2,atm3;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

SIRE_ALWAYS_INLINE uint qHash(const ImproperID &improperid)
{
    return improperid.hash();
}

#endif // SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::ImproperID);

SIRE_EXPOSE_CLASS( SireMol::ImproperID )

SIRE_END_HEADER

#endif
