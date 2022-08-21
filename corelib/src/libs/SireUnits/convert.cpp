/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <QString>
#include <QStringList>
#include <QMap>
#include <QPair>

#include <boost/shared_ptr.hpp>

#include "convert.h"
#include "dimensions.h"

#include <QDebug>

namespace SireUnits
{

namespace Dimension
{

static void appendString(int M, QString rep, QStringList &pos, QStringList &neg)
{
    if (M > 1)
    {
        pos.append( QString("%1^%2").arg(rep).arg(M) );
    }
    else if (M == 1)
    {
        pos.append(rep);
    }
    else if (M < 0)
    {
        neg.append( QString("%1%2").arg(rep).arg(M) );
    }
}

static QString getGenericUnitString(int M, int L, int T, int C,
                                    int t, int Q, int A)
{
    QStringList pos;
    QStringList neg;

    appendString(M, "M", pos, neg);
    appendString(L, "L", pos, neg);
    appendString(T, "T", pos, neg);
    appendString(C, "C", pos, neg);
    appendString(t, "t", pos, neg);
    appendString(Q, "Q", pos, neg);
    appendString(A, "A", pos, neg);

    if (pos.isEmpty() and neg.isEmpty())
        return "";
    else if (neg.isEmpty())
        return pos.join(" ");
    else if (pos.isEmpty())
        return neg.join(" ");
    else
        return QString("%1 %2").arg(pos.join(" "), neg.join(" "));
}

class DimensionKey
{
public:
    template<int M, int L, int T, int C, int t, int Q, int A>
    DimensionKey(const PhysUnit<M,L,T,C,t,Q,A> &unit)
          : M_(M), L_(L), T_(T), C_(C), t_(t), Q_(Q), A_(A)
    {}

    DimensionKey(int M, int L, int T, int C, int t, int Q, int A)
          : M_(M), L_(L), T_(T), C_(C), t_(t), Q_(Q), A_(A)
    {}

    DimensionKey(const DimensionKey &other)
          : M_(other.M_), L_(other.L_), T_(other.T_),
            C_(other.C_), t_(other.t_), Q_(other.Q_), A_(other.A_)
    {}

    ~DimensionKey()
    {}

    int M_, L_, T_, C_, t_, Q_, A_;

    bool operator==(const DimensionKey &other) const
    {
        return M_ == other.M_ and L_ == other.L_ and T_ == other.T_ and
               C_ == other.C_ and t_ == other.t_ and Q_ == other.Q_ and
               A_ == other.A_;
    }

    bool operator!=(const DimensionKey &other) const
    {
        return not this->operator==(other);
    }

    bool operator<(const DimensionKey &other) const
    {
        return (M_ < other.M_) or

               (M_ == other.M_ and L_ < other.L_) or

               (M_ == other.M_ and L_ == other.L_ and T_ < other.T_) or

               (M_ == other.M_ and L_ == other.L_ and T_ == other.T_ and
                C_ < other.C_) or

               (M_ == other.M_ and L_ == other.L_ and T_ == other.T_ and
                C_ == other.C_ and t_ < other.t_) or

               (M_ == other.M_ and L_ == other.L_ and T_ == other.T_ and
                C_ == other.C_ and t_ == other.t_ and Q_ < other.Q_) or

               (M_ == other.M_ and L_ == other.L_ and T_ == other.T_ and
                C_ == other.C_ and t_ == other.t_ and Q_ == other.Q_ and
                A_ < other.A_);
    }

    bool operator<=(const DimensionKey &other) const
    {
        return this->operator==(other) or this->operator<(other);
    }

    bool operator>=(const DimensionKey &other) const
    {
        return not this->operator<(other);
    }

    bool operator>(const DimensionKey &other) const
    {
        return not this->operator<=(other);
    }
};

boost::shared_ptr< QMap< DimensionKey,QPair<double,QString> > > default_strings;

/** Return a string representing the unit with specified dimensions */
QString getUnitString(double value,
                      int M, int L, int T, int C,
                      int t, int Q, int A)
{
    if (default_strings == 0)
    {
        boost::shared_ptr< QMap< DimensionKey, QPair<double,QString> > > strings(
                new QMap< DimensionKey,QPair<double,QString> >() );

        strings->insert( DimensionKey(kcal_per_mol),
                         QPair<double,QString>( kcal_per_mol, " kcal mol-1" ) );

        strings->insert( DimensionKey(kelvin),
                         QPair<double,QString>( kelvin, "°K"));

        strings->insert( DimensionKey(degree),
                         QPair<double,QString>( degree, "°" ) );

        strings->insert( DimensionKey(1/angstrom),
                         QPair<double,QString>( 1/angstrom, " Å-1"));

        strings->insert( DimensionKey(1/(angstrom*angstrom)),
                         QPair<double,QString>( 1/(angstrom*angstrom), " Å-2"));

        strings->insert( DimensionKey(angstrom),
                         QPair<double,QString>( angstrom, " Å" ) );

        strings->insert( DimensionKey(angstrom2),
                         QPair<double,QString>( angstrom2, " Å^2" ) );

        strings->insert( DimensionKey(angstrom3),
                         QPair<double,QString>( angstrom3, " Å^3" ) );

        strings->insert( DimensionKey(g_per_mol),
                         QPair<double,QString>( g_per_mol, " g mol-1" ) );

        strings->insert( DimensionKey(mole),
                         QPair<double,QString>( mole, " mol"));

        strings->insert( DimensionKey(mod_electron),
                         QPair<double,QString>( mod_electron, " |e|" ) );

        strings->insert( DimensionKey(picosecond),
                         QPair<double,QString>( picosecond, " ps" ) );

        strings->insert( DimensionKey(atm),
                         QPair<double,QString>( atm, " atm" ) );

        strings->insert( DimensionKey(angstrom/femtosecond),
                         QPair<double,QString>( angstrom/femtosecond, " Å fs-1" ) );

        strings->insert( DimensionKey(angstrom/(femtosecond*femtosecond)),
                         QPair<double,QString>( angstrom/(femtosecond*femtosecond), " Å fs-2" ) );

        default_strings = strings;
    }

    QMap< DimensionKey,QPair<double,QString> >::const_iterator
                 it = default_strings->constFind(DimensionKey(M,L,T,C,t,Q,A));

    if (it != default_strings->constEnd())
    {
        return QString("%1%2").arg( value / it->first )
                              .arg( it->second );
    }
    else
        return QString("%1 %2").arg(value).arg(getGenericUnitString(M,L,T,C,t,Q,A));
}

}

}

