/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2012  Christopher Woods
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

#include "gridff.h"
#include "cljpotential.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomcharges.h"
#include "atomljs.h"

#include "SireMaths/constants.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>
#include <QTime>
#include <QElapsedTimer>

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

using namespace SireMM;
using namespace SireStream;
using namespace SireVol;
using namespace SireMol;
using namespace SireMaths;
using namespace SireBase;

static const RegisterMetaType<GridFF> r_gridff;

QDataStream &operator<<(QDataStream &ds, const GridFF &gridff)
{
    writeHeader(ds, r_gridff, 5);

    SharedDataStream sds(ds);

    sds << gridff.gridbox << gridff.dimx << gridff.dimy << gridff.dimz
        << gridff.gridpot
        << gridff.buffer_size << gridff.grid_spacing
        << gridff.coul_cutoff << gridff.lj_cutoff
        << gridff.fixedatoms_coords
        << gridff.fixedatoms_params
        << gridff.closemols_coords
        << gridff.closemols_params;

    sds << quint32( gridff.oldnrgs.count() );

    for (QHash<MolNum,CLJEnergy>::const_iterator it = gridff.oldnrgs.constBegin();
         it != gridff.oldnrgs.constEnd();
         ++it)
    {
        sds << it.key() << it.value().coulomb() << it.value().lj();
    }

    //collect together the used LJParameters and write those to the stream.
    //This will let us know if we need to update the LJIDs...
    QHash<quint32,LJParameter> used_ljs;

    LJParameterDB::lock();
    for (int i=0; i<gridff.fixedatoms_params.count(); ++i)
    {
        const SireMM::detail::CLJParameter &param = gridff.fixedatoms_params.constData()[i];
        used_ljs.insert(param.ljid, LJParameterDB::_locked_getLJParameter(param.ljid));
    }
    for (int i=0; i<gridff.closemols_params.count(); ++i)
    {
        const SireMM::detail::CLJParameter &param = gridff.closemols_params.constData()[i];
        used_ljs.insert(param.ljid, LJParameterDB::_locked_getLJParameter(param.ljid));
    }
    LJParameterDB::unlock();

    sds << used_ljs;

    sds << static_cast<const InterGroupCLJFF&>(gridff);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, GridFF &gridff)
{
    VersionID v = readHeader(ds, r_gridff);

    if (v == 5)
    {
        SharedDataStream sds(ds);

        gridff = GridFF();

        gridff.fixedatoms_coords.clear();
        gridff.fixedatoms_params.clear();

        QVector<SireMM::detail::CLJParameter> fixedatoms_params;
        QVector<SireMM::detail::CLJParameter> closemols_params;
        QHash<quint32,LJParameter> used_ljs;

        sds >> gridff.gridbox >> gridff.dimx >> gridff.dimy >> gridff.dimz
            >> gridff.gridpot
            >> gridff.buffer_size >> gridff.grid_spacing
            >> gridff.coul_cutoff >> gridff.lj_cutoff
            >> gridff.fixedatoms_coords
            >> fixedatoms_params
            >> gridff.closemols_coords
            >> closemols_params;

        quint32 noldnrgs;
        sds >> noldnrgs;

        gridff.oldnrgs.clear();
        gridff.oldnrgs.reserve(noldnrgs);

        for (quint32 i=0; i<noldnrgs; ++i)
        {
            MolNum molnum;
            double cnrg, ljnrg;

            sds >> molnum >> cnrg >> ljnrg;
            gridff.oldnrgs.insert(molnum,CLJEnergy(cnrg,ljnrg));
        }

        sds >> used_ljs;

        QHash<quint32,quint32> ljidmap;

        LJParameterDB::lock();
        for (QHash<quint32,LJParameter>::const_iterator it = used_ljs.constBegin();
             it != used_ljs.constEnd();
             ++it)
        {
            quint32 newid = LJParameterDB::_locked_addLJParameter(it.value());

            if (it.key() != newid)
                ljidmap.insert(it.key(), newid);
        }
        LJParameterDB::unlock();

        if (not ljidmap.isEmpty())
        {
            //some of the LJIDs have changed
            for (int i=0; i<fixedatoms_params.count(); ++i)
            {
                SireMM::detail::CLJParameter &param = fixedatoms_params[i];
                param.ljid = ljidmap.value(param.ljid, param.ljid);
            }
            for (int i=0; i<closemols_params.count(); ++i)
            {
                SireMM::detail::CLJParameter &param = closemols_params[i];
                param.ljid = ljidmap.value(param.ljid, param.ljid);
            }
        }

        gridff.fixedatoms_params = fixedatoms_params;
        gridff.closemols_params = closemols_params;

        sds >> static_cast<InterGroupCLJFF&>(gridff);

        gridff.need_update_ljpairs = true;
    }
    else if (v == 4)
    {
        SharedDataStream sds(ds);

        gridff = GridFF();

        gridff.fixedatoms_coords.clear();
        gridff.fixedatoms_params.clear();

        QVector<SireMM::detail::CLJParameter> fixedatoms_params;
        QHash<quint32,LJParameter> used_ljs;

        sds >> gridff.buffer_size >> gridff.grid_spacing
            >> gridff.coul_cutoff >> gridff.lj_cutoff
            >> gridff.fixedatoms_coords
            >> fixedatoms_params
            >> used_ljs;

        QHash<quint32,quint32> ljidmap;

        LJParameterDB::lock();
        for (QHash<quint32,LJParameter>::const_iterator it = used_ljs.constBegin();
             it != used_ljs.constEnd();
             ++it)
        {
            quint32 newid = LJParameterDB::_locked_addLJParameter(it.value());

            if (it.key() != newid)
                ljidmap.insert(it.key(), newid);
        }
        LJParameterDB::unlock();

        if (not ljidmap.isEmpty())
        {
            //some of the LJIDs have changed
            for (int i=0; i<fixedatoms_params.count(); ++i)
            {
                SireMM::detail::CLJParameter &param = fixedatoms_params[i];
                param.ljid = ljidmap.value(param.ljid, param.ljid);
            }
        }

        gridff.fixedatoms_params = fixedatoms_params;

        sds >> static_cast<InterGroupCLJFF&>(gridff);

        gridff.need_update_ljpairs = true;
    }
    else if (v == 3)
    {
        SharedDataStream sds(ds);

        gridff = GridFF();

        gridff.fixedatoms_coords.clear();
        gridff.fixedatoms_params.clear();

        sds >> gridff.buffer_size >> gridff.grid_spacing
            >> gridff.coul_cutoff >> gridff.lj_cutoff
            >> gridff.fixedatoms_coords;

        gridff.fixedatoms_params.resize(gridff.fixedatoms_coords.count());

        LJParameterDB::lock();

        try
        {
            for (int i=0; i<gridff.fixedatoms_coords.count(); ++i)
            {
                double reduced_chg;
                LJParameter ljparam;

                sds >> reduced_chg >> ljparam;

                gridff.fixedatoms_params[i].reduced_charge = reduced_chg;
                gridff.fixedatoms_params[i].ljid = LJParameterDB::_locked_addLJParameter(ljparam);

                LJParameterDB::unlock();
            }

            gridff.need_update_ljpairs = true;
        }
        catch(...)
        {
            LJParameterDB::unlock();
            throw;
        }

        sds >> static_cast<InterGroupCLJFF&>(gridff);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);

        gridff = GridFF();

        gridff.fixedatoms_coords.clear();
        gridff.fixedatoms_params.clear();

        sds >> gridff.buffer_size >> gridff.grid_spacing
            >> gridff.coul_cutoff >> gridff.lj_cutoff
            >> static_cast<InterGroupCLJFF&>(gridff);
    }
    else
        throw version_error(v, "2,3,4,5", r_gridff, CODELOC);

    return ds;
}

/** Empty constructor */
GridFF::GridFF()
       : ConcreteProperty<GridFF,InterGroupCLJFF>(),
         buffer_size(2.5), grid_spacing(1.0),
         coul_cutoff(50), lj_cutoff(7.5)
{
    this->setSwitchingFunction(NoCutoff());
}

/** Construct a grid forcefield with a specified name */
GridFF::GridFF(const QString &name)
       : ConcreteProperty<GridFF,InterGroupCLJFF>(name),
         buffer_size(2.5), grid_spacing(1.0),
         coul_cutoff(50), lj_cutoff(7.5)
{
    this->setSwitchingFunction(NoCutoff());
}

/** Copy constructor */
GridFF::GridFF(const GridFF &other)
       : ConcreteProperty<GridFF,InterGroupCLJFF>(other),
         gridbox(other.gridbox),
         buffer_size(other.buffer_size), grid_spacing(other.grid_spacing),
         coul_cutoff(other.coul_cutoff), lj_cutoff(other.lj_cutoff),
         dimx(other.dimx), dimy(other.dimy), dimz(other.dimz),
         gridpot(other.gridpot),
         fixedatoms_coords(other.fixedatoms_coords),
         fixedatoms_params(other.fixedatoms_params),
         closemols_coords(other.closemols_coords),
         closemols_params(other.closemols_params),
         oldnrgs(other.oldnrgs)
{
    this->setSwitchingFunction(NoCutoff());
}

/** Destructor */
GridFF::~GridFF()
{}

const char* GridFF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<GridFF>() );
}

const char* GridFF::what() const
{
    return GridFF::typeName();
}

/** Copy assignment operator */
GridFF& GridFF::operator=(const GridFF &other)
{
    if (this != &other)
    {
        gridbox = other.gridbox;
        buffer_size = other.buffer_size;
        grid_spacing = other.grid_spacing;
        coul_cutoff = other.coul_cutoff;
        lj_cutoff = other.lj_cutoff;
        dimx = other.dimx;
        dimy = other.dimy;
        dimz = other.dimz;
        gridpot = other.gridpot;
        fixedatoms_coords = other.fixedatoms_coords;
        fixedatoms_params = other.fixedatoms_params;
        closemols_coords = other.closemols_coords;
        closemols_params = other.closemols_params;
        oldnrgs = other.oldnrgs;

        InterGroupCLJFF::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool GridFF::operator==(const GridFF &other) const
{
    return buffer_size == other.buffer_size and
           grid_spacing == other.grid_spacing and
           InterGroupCLJFF::operator==(other);
}

/** Comparison operator */
bool GridFF::operator!=(const GridFF &other) const
{
    return not GridFF::operator==(other);
}

GridFF* GridFF::clone() const
{
    return new GridFF(*this);
}

/** Add fixed atoms to the grid. This will copy the fixed atoms from
    the passed GridFF. This allows multiple grid forcefields to share the
    memory cost of a shared set of fixed atoms */
void GridFF::addFixedAtoms(const GridFF &other)
{
    if (other.fixedatoms_coords.isEmpty())
        return;

    if (fixedatoms_coords.isEmpty())
    {
        fixedatoms_coords = other.fixedatoms_coords;
        fixedatoms_params = other.fixedatoms_params;
    }
    else
    {
        fixedatoms_coords += other.fixedatoms_coords;
        fixedatoms_params += other.fixedatoms_params;
    }

    need_update_ljpairs = true;

    this->mustNowRecalculateFromScratch();
}

/** Add fixed atoms to the grid. These are atoms that will never change
    position or charge during the simulation, and that you wish to be
    included in the energy expression. The atoms can be placed here, and
    then do not need to be added to the simulation System. This is useful
    if you are simulating a small cutout of the system and do not want to
    have all of the atoms loaded into the system during the simulation */
void GridFF::addFixedAtoms(const MoleculeView &fixed_atoms, const PropertyMap &map)
{
    const PropertyName coords_property = map["coordinates"];
    const PropertyName chg_property = map["charge"];
    const PropertyName lj_property = map["LJ"];

    const QVector<Vector> coords = fixed_atoms.molecule().property(coords_property)
                                      .asA<AtomCoords>().toVector(fixed_atoms.selection());

    const QVector<SireUnits::Dimension::Charge> charges =
                            fixed_atoms.molecule().property(chg_property)
                                      .asA<AtomCharges>().toVector(fixed_atoms.selection());

    const QVector<LJParameter> ljs = fixed_atoms.molecule().property(lj_property)
                                      .asA<AtomLJs>().toVector(fixed_atoms.selection());

    int nats = coords.count();

    fixedatoms_coords.reserve( fixedatoms_coords.count() + nats );
    fixedatoms_params.reserve( fixedatoms_params.count() + nats );

    fixedatoms_coords += coords;

    LJParameterDB::lock();

    const double sqrt_one_over_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);

    for (int i=0; i<nats; ++i)
    {
        SireMM::detail::CLJParameter cljparam;

        cljparam.reduced_charge = charges[i] * sqrt_one_over_4pieps0;
        cljparam.ljid = LJParameterDB::_locked_addLJParameter(ljs[i]);

        fixedatoms_params.append(cljparam);
    }

    LJParameterDB::unlock();

    fixedatoms_coords.squeeze();
    fixedatoms_params.squeeze();

    need_update_ljpairs = true;

    this->mustNowRecalculateFromScratch();
}

/** Add fixed atoms to the grid. These are atoms that will never change
    position or charge during the simulation, and that you wish to be
    included in the energy expression. The atoms can be placed here, and
    then do not need to be added to the simulation System. This is useful
    if you are simulating a small cutout of the system and do not want to
    have all of the atoms loaded into the system during the simulation */
void GridFF::addFixedAtoms(const SireMol::Molecules &fixed_atoms, const PropertyMap &map)
{
    const PropertyName coords_property = map["coordinates"];
    const PropertyName chg_property = map["charge"];
    const PropertyName lj_property = map["LJ"];

    const double sqrt_one_over_4pieps0 = std::sqrt(SireUnits::one_over_four_pi_eps0);

    LJParameterDB::lock();

    for (SireMol::Molecules::const_iterator it = fixed_atoms.constBegin();
         it != fixed_atoms.constEnd();
         ++it)
    {
        SireMol::Molecule mol = it->molecule();
        AtomSelection selection = it->selection();

        const QVector<Vector> coords = mol.property(coords_property)
                                        .asA<AtomCoords>().toVector(selection);

        const QVector<SireUnits::Dimension::Charge> charges = mol.property(chg_property)
                                        .asA<AtomCharges>().toVector(selection);

        const QVector<LJParameter> ljs = mol.property(lj_property)
                                            .asA<AtomLJs>().toVector(selection);

        int nats = coords.count();

        fixedatoms_coords.reserve( fixedatoms_coords.count() + nats );
        fixedatoms_params.reserve( fixedatoms_params.count() + nats );

        fixedatoms_coords += coords;

        for (int i=0; i<nats; ++i)
        {
            SireMM::detail::CLJParameter cljparam;

            cljparam.reduced_charge = charges[i] * sqrt_one_over_4pieps0;
            cljparam.ljid = LJParameterDB::_locked_addLJParameter(ljs[i]);

            fixedatoms_params.append(cljparam);
        }
    }

    LJParameterDB::unlock();

    need_update_ljpairs = true;

    fixedatoms_coords.squeeze();
    fixedatoms_params.squeeze();

    this->mustNowRecalculateFromScratch();
}

/** Add all of the atoms in the molecules in the passed molecule group to the set
    of fixed atoms */
void GridFF::addFixedAtoms(const MoleculeGroup &group, const PropertyMap &map)
{
    this->addFixedAtoms(group.molecules(), map);
}

/** Set the buffer when building the grid. This adds a buffer space
    around the grid when it is built, to try to reduce the number of
    times it needs to be rebuilt */
void GridFF::setBuffer(SireUnits::Dimension::Length buffer)
{
    buffer_size = buffer.value();
}

/** Set the grid spacing (the distance between grid points). The
    smaller the spacing the more memory is required to hold the grid,
    but the more accurate the energy */
void GridFF::setGridSpacing(SireUnits::Dimension::Length spacing)
{
    if (grid_spacing != spacing.value())
    {
        grid_spacing = spacing.value();
        this->mustNowRecalculateFromScratch();
    }
}

/** Set the cutoff for the coulomb energy - this can be greater
    than the box size as multiple periodic images can be used */
void GridFF::setCoulombCutoff(SireUnits::Dimension::Length cutoff)
{
    if (coul_cutoff != cutoff.value())
    {
        coul_cutoff = cutoff.value();
        this->mustNowRecalculateFromScratch();
    }
}

/** Set the cutoff for the LJ energy - this can be greater than
    the box size as multiple periodic images can be used */
void GridFF::setLJCutoff(SireUnits::Dimension::Length cutoff)
{
    if (lj_cutoff != cutoff.value())
    {
        lj_cutoff = cutoff.value();
        this->mustNowRecalculateFromScratch();
    }
}

/** Turn on or off use of the force shifted potential */
bool GridFF::setShiftElectrostatics(bool on)
{
    if (InterGroupCLJFF::setShiftElectrostatics(on))
    {
        this->mustNowRecalculateFromScratch();
        return true;
    }
    else
        return false;
}

/** Turn on or off the use of the reaction field */
bool GridFF::setUseReactionField(bool on)
{
    if (InterGroupCLJFF::setUseReactionField(on))
    {
        this->mustNowRecalculateFromScratch();
        return true;
    }
    else
        return false;
}

/** Set the dielectric constant to use with the reaction field potential */
bool GridFF::setReactionFieldDielectric(double dielectric)
{
    if (InterGroupCLJFF::setReactionFieldDielectric(dielectric))
    {
        this->mustNowRecalculateFromScratch();
        return true;
    }

    return false;
}

/** Return the buffer size used when building grids */
SireUnits::Dimension::Length GridFF::buffer() const
{
    return SireUnits::Dimension::Length(buffer_size);
}

/** Return the spacing between grid points */
SireUnits::Dimension::Length GridFF::spacing() const
{
    return SireUnits::Dimension::Length(grid_spacing);
}

/** Return the cutoff for the coulomb energy */
SireUnits::Dimension::Length GridFF::coulombCutoff() const
{
    return SireUnits::Dimension::Length(coul_cutoff);
}

/** Return the cutoff for the LJ energy */
SireUnits::Dimension::Length GridFF::ljCutoff() const
{
    return SireUnits::Dimension::Length(lj_cutoff);
}

SIRE_ALWAYS_INLINE GridFF::Vector4::Vector4(const Vector &v, double chg)
              : x(v.x()), y(v.y()), z(v.z()), q(chg)
{}

void GridFF::appendTo(QVector<GridFF::Vector4> &coords_and_charges,
                      const Vector *coords, const SireMM::detail::CLJParameter *params,
                      int nats)
{
    if (nats > 0)
    {
        coords_and_charges.reserve( coords_and_charges.count() + nats );

        for (int i=0; i<nats; ++i)
        {
            double q = params[i].reduced_charge;

            if (q != 0)
            {
                coords_and_charges.append( Vector4(coords[i],q) );
            }
        }
    }
}

#ifdef SIRE_USE_SSE
    SIRE_ALWAYS_INLINE QString toString(const __m128d &sseval)
    {
        return QString("{ %1, %2 }").arg(*((const double*)&sseval))
                                    .arg(*( ((const double*)&sseval) + 1));
    }
#endif

void GridFF::addToGrid(const QVector<GridFF::Vector4> &coords_and_charges)
{
    double *pot = gridpot.data();
    const Vector4 *array = coords_and_charges.constData();

    Vector minpoint = gridbox.minCoords();

    //loop over all grid points
    const int npts = dimx*dimy*dimz;
    const int nats = coords_and_charges.count();

    const double Rc = coul_cutoff;

    QElapsedTimer t;
    t.start();

    if (shiftElectrostatics())
    {
        const double one_over_Rc = double(1) / Rc;
        const double one_over_Rc2 = double(1) / (Rc*Rc);

        #ifdef SIRE_USE_SSE
        {
            const int remainder = npts % 2;

            quint32 i0=0;
            quint32 j0=0;
            quint32 k0=0;

            quint32 i1=0;
            quint32 j1=0;
            quint32 k1=1;

            double gx0 = minpoint.x();
            double gy0 = minpoint.y();
            double gz0 = minpoint.z();

            double gx1 = minpoint.x();
            double gy1 = minpoint.y();
            double gz1 = minpoint.z()+grid_spacing;

            const __m128d sse_one = { 1.0, 1.0 };

            const __m128d sse_Rc = { Rc, Rc };
            const __m128d sse_one_over_Rc = { one_over_Rc, one_over_Rc };
            const __m128d sse_one_over_Rc2 = { one_over_Rc2, one_over_Rc2 };

            for (int ipt=0; ipt<(npts-1); ipt+=2)
            {
                //set the sse values - note that _mm_set_pd is backwards,
                //so the lower value is gx0, not gx1
                const __m128d sse_gx = _mm_set_pd(gx1,gx0);
                const __m128d sse_gy = _mm_set_pd(gy1,gy0);
                const __m128d sse_gz = _mm_set_pd(gz1,gz0);

                __m128d sse_cnrg = { 0.0, 0.0 };

                //loop through each atom
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const __m128d sse_dx = _mm_sub_pd(_mm_set1_pd(c.x), sse_gx);
                    const __m128d sse_dy = _mm_sub_pd(_mm_set1_pd(c.y), sse_gy);
                    const __m128d sse_dz = _mm_sub_pd(_mm_set1_pd(c.z), sse_gz);

                    const __m128d sse_r = _mm_sqrt_pd(
                                        _mm_add_pd( _mm_mul_pd(sse_dx,sse_dx),
                                          _mm_add_pd( _mm_mul_pd(sse_dy,sse_dy),
                                                      _mm_mul_pd(sse_dz,sse_dz) ) ) );

                    const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                    __m128d nrg = _mm_sub_pd(sse_r, sse_Rc);
                    nrg = _mm_mul_pd(nrg, sse_one_over_Rc2);
                    nrg = _mm_add_pd(nrg, sse_one_over_r);
                    nrg = _mm_sub_pd(nrg, sse_one_over_Rc);
                    nrg = _mm_mul_pd(nrg, _mm_set1_pd(c.q));

                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);
                    nrg = _mm_and_pd(nrg, sse_in_cutoff);

                    sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                }

                pot[ipt] += *((const double*)&sse_cnrg);
                pot[ipt+1] += *( ((const double*)&sse_cnrg) + 1 );

                for (int ii=0; ii<2; ++ii)
                {
                    //advance the indicies (twice, as two grid points per iteration)
                    k0 += 1;
                    gz0 += grid_spacing;

                    if (k0 == dimz)
                    {
                        k0 = 0;
                        gz0 = minpoint.z();

                        j0 += 1;
                        gy0 += grid_spacing;

                        if (j0 == dimy)
                        {
                            j0 = 0;
                            gy0 = minpoint.y();

                            i0 += 1;
                            gx0 += grid_spacing;

                            if (i0 == dimx)
                            {
                                i0 = 0;
                                gx0 = minpoint.x();
                            }
                        }
                    }

                    k1 += 1;
                    gz1 += grid_spacing;

                    if (k1 == dimz)
                    {
                        k1 = 0;
                        gz1 = minpoint.z();

                        j1 += 1;
                        gy1 += grid_spacing;

                        if (j1 == dimy)
                        {
                            j1 = 0;
                            gy1 = minpoint.y();

                            i1 += 1;
                            gx1 += grid_spacing;

                            if (i1 == dimx)
                            {
                                i1 = 0;
                                gx1 = minpoint.x();
                            }
                        }
                    }
                }
            }

            if (remainder)
            {
                //we need to process the last grid point
                const double gx = minpoint.x() + (dimx-1)*grid_spacing;
                const double gy = minpoint.y() + (dimy-1)*grid_spacing;
                const double gz = minpoint.z() + (dimz-1)*grid_spacing;

                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double one_over_r = double(1) / r;

                    const double in_cutoff = (r < Rc);
                    total += in_cutoff * c.q * (one_over_r - one_over_Rc +
                                                   one_over_Rc2*(r - Rc));
                }

                pot[npts-1] += total;
            }
        }
        #else
        {
            int i=0;
            int j=0;
            int k=0;

            double gx = minpoint.x();
            double gy = minpoint.y();
            double gz = minpoint.z();

            for (int ipt=0; ipt<npts; ++ipt)
            {
                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double one_over_r = double(1) / r;

                    const double in_cutoff = (r < Rc);

                    total += in_cutoff * c.q * (one_over_r - one_over_Rc +
                                                one_over_Rc2*(r - Rc));
                }

                pot[ipt] += total;

                k += 1;
                gz += grid_spacing;

                if (k == dimz)
                {
                    k = 0;
                    gz = minpoint.z();

                    j += 1;
                    gy += grid_spacing;

                    if (j == dimy)
                    {
                        j = 0;
                        gy = minpoint.y();

                        i += 1;
                        gx += grid_spacing;

                        if (i == dimx)
                        {
                            i = 0;
                            gx = minpoint.x();
                        }
                    }
                }
            }
        }
        #endif
    }
    else if (useReactionField())
    {
        //use the reaction field potential
        double rf_dielectric = reactionFieldDielectric();
        const double k_rf = (1.0 / pow_3(Rc)) * ( (rf_dielectric-1) /
                                                  (2*rf_dielectric + 1) );
        const double c_rf = (1.0 / Rc) * ( (3*rf_dielectric) /
                                           (2*rf_dielectric + 1) );

        #ifdef SIRE_USE_SSE
        {
            const int remainder = npts % 2;

            quint32 i0=0;
            quint32 j0=0;
            quint32 k0=0;

            quint32 i1=0;
            quint32 j1=0;
            quint32 k1=1;

            double gx0 = minpoint.x();
            double gy0 = minpoint.y();
            double gz0 = minpoint.z();

            double gx1 = minpoint.x();
            double gy1 = minpoint.y();
            double gz1 = minpoint.z()+grid_spacing;

            const __m128d sse_one = { 1.0, 1.0 };

            const __m128d sse_Rc = { Rc, Rc };
            const __m128d sse_k_rf = { k_rf, k_rf };
            const __m128d sse_c_rf = { c_rf, c_rf };

            for (int ipt=0; ipt<(npts-1); ipt+=2)
            {
                //set the sse values - note that _mm_set_pd is backwards,
                //so the lower value is gx0, not gx1
                const __m128d sse_gx = _mm_set_pd(gx1,gx0);
                const __m128d sse_gy = _mm_set_pd(gy1,gy0);
                const __m128d sse_gz = _mm_set_pd(gz1,gz0);

                __m128d sse_cnrg = { 0.0, 0.0 };

                //loop through each atom
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const __m128d sse_dx = _mm_sub_pd(_mm_set1_pd(c.x), sse_gx);
                    const __m128d sse_dy = _mm_sub_pd(_mm_set1_pd(c.y), sse_gy);
                    const __m128d sse_dz = _mm_sub_pd(_mm_set1_pd(c.z), sse_gz);

                    const __m128d sse_r = _mm_sqrt_pd(
                                        _mm_add_pd( _mm_mul_pd(sse_dx,sse_dx),
                                          _mm_add_pd( _mm_mul_pd(sse_dy,sse_dy),
                                                      _mm_mul_pd(sse_dz,sse_dz) ) ) );

                    const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                    __m128d nrg = _mm_mul_pd(sse_r, sse_r);
                    nrg = _mm_mul_pd(nrg, sse_k_rf);
                    nrg = _mm_sub_pd(nrg, sse_c_rf);
                    nrg = _mm_add_pd(nrg, sse_one_over_r);
                    nrg = _mm_mul_pd(nrg, _mm_set1_pd(c.q));

                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);
                    nrg = _mm_and_pd(nrg, sse_in_cutoff);

                    sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                }

                pot[ipt] += *((const double*)&sse_cnrg);
                pot[ipt+1] += *( ((const double*)&sse_cnrg) + 1 );

                for (int ii=0; ii<2; ++ii)
                {
                    //advance the indicies (twice, as two grid points per iteration)
                    k0 += 1;
                    gz0 += grid_spacing;

                    if (k0 == dimz)
                    {
                        k0 = 0;
                        gz0 = minpoint.z();

                        j0 += 1;
                        gy0 += grid_spacing;

                        if (j0 == dimy)
                        {
                            j0 = 0;
                            gy0 = minpoint.y();

                            i0 += 1;
                            gx0 += grid_spacing;

                            if (i0 == dimx)
                            {
                                i0 = 0;
                                gx0 = minpoint.x();
                            }
                        }
                    }

                    k1 += 1;
                    gz1 += grid_spacing;

                    if (k1 == dimz)
                    {
                        k1 = 0;
                        gz1 = minpoint.z();

                        j1 += 1;
                        gy1 += grid_spacing;

                        if (j1 == dimy)
                        {
                            j1 = 0;
                            gy1 = minpoint.y();

                            i1 += 1;
                            gx1 += grid_spacing;

                            if (i1 == dimx)
                            {
                                i1 = 0;
                                gx1 = minpoint.x();
                            }
                        }
                    }
                }
            }

            if (remainder)
            {
                //we need to process the last grid point
                const double gx = minpoint.x() + (dimx-1)*grid_spacing;
                const double gy = minpoint.y() + (dimy-1)*grid_spacing;
                const double gz = minpoint.z() + (dimz-1)*grid_spacing;

                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double one_over_r = double(1) / r;

                    const double in_cutoff = (r < Rc);
                    total += in_cutoff * c.q * (one_over_r + k_rf*r*r - c_rf);
                }

                pot[npts-1] += total;
            }
        }
        #else
        {
            int i=0;
            int j=0;
            int k=0;

            double gx = minpoint.x();
            double gy = minpoint.y();
            double gz = minpoint.z();

            for (int ipt=0; ipt<npts; ++ipt)
            {
                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double one_over_r = double(1) / r;

                    const double in_cutoff = (r < Rc);

                    total += in_cutoff * c.q * (one_over_r + k_rf*r*r - c_rf);
                }

                pot[ipt] += total;

                k += 1;
                gz += grid_spacing;

                if (k == dimz)
                {
                    k = 0;
                    gz = minpoint.z();

                    j += 1;
                    gy += grid_spacing;

                    if (j == dimy)
                    {
                        j = 0;
                        gy = minpoint.y();

                        i += 1;
                        gx += grid_spacing;

                        if (i == dimx)
                        {
                            i = 0;
                            gx = minpoint.x();
                        }
                    }
                }
            }
        }
        #endif
    }
    else
    {
        //we use a simple atom-based cutoff
        #ifdef SIRE_USE_SSE
        {
            const int remainder = npts % 2;

            quint32 i0=0;
            quint32 j0=0;
            quint32 k0=0;

            quint32 i1=0;
            quint32 j1=0;
            quint32 k1=1;

            double gx0 = minpoint.x();
            double gy0 = minpoint.y();
            double gz0 = minpoint.z();

            double gx1 = minpoint.x();
            double gy1 = minpoint.y();
            double gz1 = minpoint.z()+grid_spacing;

            const __m128d sse_one = { 1.0, 1.0 };

            const __m128d sse_Rc = { Rc, Rc };

            for (int ipt=0; ipt<(npts-1); ipt+=2)
            {
                //set the sse values - note that _mm_set_pd is backwards,
                //so the lower value is gx0, not gx1
                const __m128d sse_gx = _mm_set_pd(gx1,gx0);
                const __m128d sse_gy = _mm_set_pd(gy1,gy0);
                const __m128d sse_gz = _mm_set_pd(gz1,gz0);

                __m128d sse_cnrg = { 0.0, 0.0 };

                //loop through each atom
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const __m128d sse_dx = _mm_sub_pd(_mm_set1_pd(c.x), sse_gx);
                    const __m128d sse_dy = _mm_sub_pd(_mm_set1_pd(c.y), sse_gy);
                    const __m128d sse_dz = _mm_sub_pd(_mm_set1_pd(c.z), sse_gz);

                    const __m128d sse_r = _mm_sqrt_pd(
                                        _mm_add_pd( _mm_mul_pd(sse_dx,sse_dx),
                                          _mm_add_pd( _mm_mul_pd(sse_dy,sse_dy),
                                                      _mm_mul_pd(sse_dz,sse_dz) ) ) );

                    const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                    __m128d nrg = _mm_mul_pd(sse_one_over_r, _mm_set1_pd(c.q));

                    const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);
                    nrg = _mm_and_pd(nrg, sse_in_cutoff);

                    sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                }

                pot[ipt] += *((const double*)&sse_cnrg);
                pot[ipt+1] += *( ((const double*)&sse_cnrg) + 1 );

                for (int ii=0; ii<2; ++ii)
                {
                    //advance the indicies (twice, as two grid points per iteration)
                    k0 += 1;
                    gz0 += grid_spacing;

                    if (k0 == dimz)
                    {
                        k0 = 0;
                        gz0 = minpoint.z();

                        j0 += 1;
                        gy0 += grid_spacing;

                        if (j0 == dimy)
                        {
                            j0 = 0;
                            gy0 = minpoint.y();

                            i0 += 1;
                            gx0 += grid_spacing;

                            if (i0 == dimx)
                            {
                                i0 = 0;
                                gx0 = minpoint.x();
                            }
                        }
                    }

                    k1 += 1;
                    gz1 += grid_spacing;

                    if (k1 == dimz)
                    {
                        k1 = 0;
                        gz1 = minpoint.z();

                        j1 += 1;
                        gy1 += grid_spacing;

                        if (j1 == dimy)
                        {
                            j1 = 0;
                            gy1 = minpoint.y();

                            i1 += 1;
                            gx1 += grid_spacing;

                            if (i1 == dimx)
                            {
                                i1 = 0;
                                gx1 = minpoint.x();
                            }
                        }
                    }
                }
            }

            if (remainder)
            {
                //we need to process the last grid point
                const double gx = minpoint.x() + (dimx-1)*grid_spacing;
                const double gy = minpoint.y() + (dimy-1)*grid_spacing;
                const double gz = minpoint.z() + (dimz-1)*grid_spacing;

                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);
                    const double one_over_r = double(1) / r;

                    const double in_cutoff = (r < Rc);
                    total += in_cutoff * c.q * one_over_r;
                }

                pot[npts-1] += total;
            }
        }
        #else
        {
            int i=0;
            int j=0;
            int k=0;

            double gx = minpoint.x();
            double gy = minpoint.y();
            double gz = minpoint.z();

            for (int ipt=0; ipt<npts; ++ipt)
            {
                double total = 0;

                //loop through each atom...
                for (int iat=0; iat<nats; ++iat)
                {
                    const Vector4 &c = array[iat];

                    const double dx = c.x - gx;
                    const double dy = c.y - gy;
                    const double dz = c.z - gz;

                    const double r = std::sqrt(dx*dx + dy*dy + dz*dz);

                    if (r < Rc)
                    {
                        total += c.q / r;
                    }
                }

                pot[ipt] += total;

                k += 1;
                gz += grid_spacing;

                if (k == dimz)
                {
                    k = 0;
                    gz = minpoint.z();

                    j += 1;
                    gy += grid_spacing;

                    if (j == dimy)
                    {
                        j = 0;
                        gy = minpoint.y();

                        i += 1;
                        gx += grid_spacing;

                        if (i == dimx)
                        {
                            i = 0;
                            gx = minpoint.x();
                        }
                    }
                }
            }
        }
        #endif
    }

    //int ms = t.elapsed();

    //qDebug() << "Added" << nats << "more atoms to" << npts << "grid points in"
    //         << ms << "ms";
}

SIRE_ALWAYS_INLINE double getDist(double p, double minp, double maxp)
{
    if (p < minp)
        return minp - p;
    else if (p > maxp)
        return p - maxp;
    else
        return std::min(p - minp, maxp - p);
}

/** Return the minimum distance between a point and the passed AABox. This
    returns 0 if the point is inside the box */
double minimumDistanceToGrid(const Vector &coords, const AABox &gridbox)
{
    const Vector mincoords = gridbox.minCoords();
    const Vector maxcoords = gridbox.maxCoords();

    int ix(0), iy(0), iz(0);
    double dx(0), dy(0), dz(0);

    if (coords.x() < mincoords.x())
    {
        ix = -1;
        dx = mincoords.x() - coords.x();
    }
    else if (coords.x() > maxcoords.x())
    {
        ix = 1;
        dx = coords.x() - maxcoords.x();
    }

    if (coords.y() < mincoords.y())
    {
        iy = -1;
        dy = mincoords.y() - coords.y();
    }
    else if (coords.y() > maxcoords.y())
    {
        iy = 1;
        dy = coords.y() - maxcoords.y();
    }

    if (coords.z() < mincoords.z())
    {
        iz = -1;
        dz = mincoords.z() - coords.z();
    }
    else if (coords.z() > maxcoords.z())
    {
        iz = 1;
        dz = coords.z() - maxcoords.z();
    }

    if (ix == 0 and iy == 0 and iz == 0)
        //the box contains the point
        return 0;
    else if (ix == 0 and iy == 0)
    {
        //the point is above or below the xy face
        return dz;
    }
    else if (ix == 0 and iz == 0)
    {
        //the point is above or below the xz face
        return dy;
    }
    else if (iy == 0 and iz == 0)
    {
        //the point is above or below the yz face
        return dx;
    }
    else if (ix == 0)
    {
        //the point is next to the yz corner
        return std::sqrt(dy*dy + dz*dz);
    }
    else if (iy == 0)
    {
        //the point is next to the xz corner
        return std::sqrt(dx*dx + dz*dz);
    }
    else if (iz == 0)
    {
        //the point is next to the xy corner
        return std::sqrt(dx*dx + dy*dy);
    }
    else
    {
        //the point is off the 3D corner
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
}

/** Internal function used to rebuild the coulomb potential grid */
void GridFF::rebuildGrid()
{
    //QTime t;

    if (mols[0].moleculesByIndex().isEmpty())
        return;

    //find the bounding box that contains all of the atoms from group0
    AABox group0_box;

    for (ChunkedVector<CLJMolecule>::const_iterator
                                it = mols[0].moleculesByIndex().constBegin();
         it != mols[0].moleculesByIndex().constEnd();
         ++it)
    {
        const CoordGroupArray &coords = (*it).coordinates();

        for (int i=0; i<coords.count(); ++i)
        {
            if (group0_box.isNull())
            {
                group0_box = coords.constData()[i].aaBox();
            }
            else
            {
                group0_box += coords.constData()[i].aaBox();
            }
        }
    }

    //now add the buffer region onto the box - this extends the grid by
    //buffer_size in all dimensions
    if (buffer_size > 0)
        gridbox = AABox::from( group0_box.minCoords() - Vector(buffer_size),
                               group0_box.maxCoords() + Vector(buffer_size) );

    //now work out how many gridpoints are needed in each dimension
    if (grid_spacing <= 0)
        grid_spacing = 0.5;

    Vector boxsize = gridbox.maxCoords() - gridbox.minCoords();

    dimx = 2 + int(boxsize.x() / grid_spacing);
    dimy = 2 + int(boxsize.y() / grid_spacing);
    dimz = 2 + int(boxsize.z() / grid_spacing);

    const quint32 MAX_DIM = 250;

    while ( dimx > MAX_DIM or dimy > MAX_DIM or dimz > MAX_DIM )
    {
        double grid_x = boxsize.x() / MAX_DIM;
        double grid_y = boxsize.y() / MAX_DIM;
        double grid_z = boxsize.z() / MAX_DIM;

        grid_spacing = qMax( grid_x, qMax(grid_y,grid_z) );

        dimx = 1 + int(boxsize.x() / grid_spacing);
        dimy = 1 + int(boxsize.y() / grid_spacing);
        dimz = 1 + int(boxsize.z() / grid_spacing);
    }

    Vector maxpoint = gridbox.minCoords() + Vector( (dimx-1) * grid_spacing,
                                                    (dimy-1) * grid_spacing,
                                                    (dimz-1) * grid_spacing );

    gridbox = AABox::from(gridbox.minCoords(), maxpoint);
    const Vector &grid_center = gridbox.center();

    //create space for the grid
    gridpot = QVector<double>(dimx*dimy*dimz, 0.0);
    gridpot.squeeze();

    closemols_coords.clear();
    closemols_params.clear();

    //now build the grid - we take atoms that are within the LJ cutoff
    //and calculate those manually (saving them in the closemols list), and
    //only add coulomb potentials of CoordGroups that are beyond the LJ cutoff
    const ChunkedVector<CLJMolecule> &cljmols = mols[1].moleculesByIndex();

    const Space &spce = this->space();

    QVector<Vector4> far_mols;

    int atomcount = 0;
    int gridcount = 0;

    closemols_coords.reserve(cljmols.count() + (fixedatoms_coords.count() / 2));
    closemols_params.reserve(cljmols.count() + (fixedatoms_coords.count() / 2));

    if (fixedatoms_coords.count() > 0)
    {
        for (int i=0; i<fixedatoms_coords.count(); ++i)
        {
            Vector coords = fixedatoms_coords.constData()[i];

            if (spce.isPeriodic())
            {
                coords = spce.getMinimumImage(coords, grid_center);
            }

            const SireMM::detail::CLJParameter &params = fixedatoms_params.constData()[i];

            if (params.reduced_charge == 0 and params.ljid == 0)
                continue;

            //calculate the closest distance between this point and the grid
            double dist = ::minimumDistanceToGrid(coords, gridbox);

            //only explicitly evaluate points within the LJ cutoff of the grid
            if (dist < lj_cutoff)
            {
                atomcount += 1;
                closemols_coords.append(coords);
                closemols_params.append(params);
            }
            //all other points are evaluated using the grid
            else
            {
                if (params.reduced_charge != 0)
                {
                    far_mols.append( Vector4(coords, params.reduced_charge) );

                    if (far_mols.count() > 1023)
                    {
                        addToGrid(far_mols);
                        gridcount += far_mols.count();
                        far_mols.clear();
                        //qDebug() << "Added" << i+1 << "of" << fixedatoms_coords.count()
                        //         << "fixed atoms to the grid...";
                    }
                }
            }
        }

        addToGrid(far_mols);
        gridcount += far_mols.count();
        far_mols.clear();
    }

    if (not cljmols.isEmpty())
    {
        int nmols = 0;

        for (ChunkedVector<CLJMolecule>::const_iterator it = cljmols.constBegin();
             it != cljmols.constEnd();
             ++it)
        {
            nmols += 1;

            const CLJMolecule &cljmol = *it;

            //loop through each CutGroup of this molecule
            const int ngroups = cljmol.coordinates().count();

            const CoordGroup *groups_array = cljmol.coordinates().constData();

            const CLJParameters::Array *params_array
                                    = cljmol.parameters().atomicParameters().constData();

            for (int igroup=0; igroup<ngroups; ++igroup)
            {
                CoordGroup coordgroup = groups_array[igroup];
                const CLJParameters::Array &paramsgroup = params_array[igroup];

                for (int i=0; i<coordgroup.count(); ++i)
                {
                    Vector coords = coordgroup.constData()[i];

                    if (spce.isPeriodic())
                    {
                        coords = spce.getMinimumImage(coords, grid_center);
                    }

                    const SireMM::detail::CLJParameter &params = paramsgroup.constData()[i];

                    if (params.reduced_charge == 0 and params.ljid == 0)
                        continue;

                    //calculate the closest distance between this point and the grid
                    double dist = ::minimumDistanceToGrid(coords, gridbox);

                    //only explicitly evaluate points within the LJ cutoff of the grid
                    if (dist < lj_cutoff)
                    {
                        atomcount += 1;
                        closemols_coords.append(coords);
                        closemols_params.append(params);
                    }
                    //all other points are evaluated using the grid
                    else
                    {
                        if (params.reduced_charge != 0)
                        {
                            far_mols.append( Vector4(coords, params.reduced_charge) );

                            if (far_mols.count() > 1023)
                            {
                                addToGrid(far_mols);
                                gridcount += far_mols.count();
                                far_mols.clear();
                            }
                        }
                    }
                }
            }
        }

        addToGrid(far_mols);
        gridcount += far_mols.count();
        far_mols.clear();
    }

    closemols_coords.squeeze();
    closemols_params.squeeze();
}

void GridFF::calculateEnergy(const CoordGroup &coords0,
                             const GridFF::CLJParameters::Array &params0,
                             double &cnrg, double &ljnrg)
{
    double icnrg = 0;
    double iljnrg = 0;

    const int nats0 = coords0.count();

    const Vector *coords0_array = coords0.constData();
    const detail::CLJParameter *params0_array = params0.constData();

    const Vector *coords1_array = closemols_coords.constData();
    const detail::CLJParameter *params1_array = closemols_params.constData();

    const int nats1 = closemols_coords.count();

    BOOST_ASSERT( closemols_coords.count() == closemols_params.count() );

    if (nats1 > 0)
    {
        if (shiftElectrostatics())
        {
            //use the force-shifted cutoff with an alpha value of 0
            const double Rc = coul_cutoff;
            const double Rlj = lj_cutoff;

            const double one_over_Rc = double(1) / Rc;
            const double one_over_Rc2 = double(1) / (Rc*Rc);

            #ifdef SIRE_USE_SSE
            {
                const int remainder = nats1 % 2;

                __m128d sse_cnrg = { 0, 0 };
                __m128d sse_ljnrg = { 0, 0 };

                const __m128d sse_one = { 1.0, 1.0 };

                const __m128d sse_Rc = { Rc, Rc };
                const __m128d sse_Rlj = { Rlj, Rlj };

                const __m128d sse_one_over_Rc = { one_over_Rc, one_over_Rc };
                const __m128d sse_one_over_Rc2 = { one_over_Rc2, one_over_Rc2 };

                for (int j=0; j<nats1-1; j += 2)
                {
                    const Parameter &param10 = params1_array[j];
                    const Parameter &param11 = params1_array[j+1];

                    const Vector &c10 = coords1_array[j];
                    const Vector &c11 = coords1_array[j+1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];

                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge,
                                                       param0.reduced_charge );

                        const Vector &c0 = coords0_array[i];

                        __m128d sse_r;
                        {
                            __m128d dx = _mm_set_pd( c0.x() - c10.x(), c0.x() - c11.x() );
                            __m128d dy = _mm_set_pd( c0.y() - c10.y(), c0.y() - c11.y() );
                            __m128d dz = _mm_set_pd( c0.z() - c10.z(), c0.z() - c11.z() );

                            __m128d dx2 = _mm_mul_pd(dx, dx);
                            __m128d dy2 = _mm_mul_pd(dy, dy);
                            __m128d dz2 = _mm_mul_pd(dz, dz);

                            __m128d sse_r2 = _mm_add_pd(dx2,dy2);
                            sse_r2 = _mm_add_pd(sse_r2, dz2);

                            sse_r = _mm_sqrt_pd(sse_r2);
                        }

                        const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                        //calculate the coulomb energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);
                            __m128d nrg = _mm_sub_pd(sse_r, sse_Rc);
                            nrg = _mm_mul_pd(nrg, sse_one_over_Rc2);
                            nrg = _mm_add_pd(nrg, sse_one_over_r);
                            nrg = _mm_sub_pd(nrg, sse_one_over_Rc);

                            __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                          param11.reduced_charge );

                            sse_chg = _mm_mul_pd(sse_chg, sse_chg0);

                            nrg = _mm_mul_pd(sse_chg, nrg);

                            nrg = _mm_and_pd(nrg, sse_in_cutoff);

                            sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                        }

                        //calculate the LJ energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param10.ljid)];

                            const LJPair &ljpair1 = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param11.ljid)];

                            const __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                            const __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                                ljpair1.epsilon() );

                            //calculate (sigma/r)^6 and (sigma/r)^12
                            __m128d sse_sig_over_r2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                            sse_sig_over_r2 = _mm_mul_pd( sse_sig_over_r2,
                                                          sse_sig_over_r2 );

                            __m128d sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r2,
                                                                 sse_sig_over_r2);

                            sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r6,
                                                         sse_sig_over_r2);

                            __m128d sse_sig_over_r12 = _mm_mul_pd(sse_sig_over_r6,
                                                                  sse_sig_over_r6);

                            //calculate LJ energy (the factor of 4 is added later)
                            __m128d tmp = _mm_sub_pd(sse_sig_over_r12,
                                                     sse_sig_over_r6);

                            tmp = _mm_mul_pd(tmp, sse_eps);
                            tmp = _mm_and_pd(tmp, sse_in_cutoff);
                            sse_ljnrg = _mm_add_pd(sse_ljnrg, tmp);
                        }
                    }
                }

                if (remainder == 1)
                {
                    const Vector &c1 = coords1_array[nats1-1];
                    const Parameter &param1 = params1_array[nats1-1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Vector &c0 = coords0_array[i];
                        const Parameter &param0 = params0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += param0.reduced_charge * param1.reduced_charge
                                        * (one_over_r - one_over_Rc + one_over_Rc2*(r-Rc));
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }

                icnrg += *((const double*)&sse_cnrg) +
                         *( ((const double*)&sse_cnrg) + 1 );

                iljnrg += *((const double*)&sse_ljnrg) +
                          *( ((const double*)&sse_ljnrg) + 1 );
            }
            #else
            {
                //we use the force-shifted cutoff potential described in equation 18
                //in Fennell and Gezelter, J. Chem. Phys., 124, 234104, 2006
                //We use alpha=0, as I have seen that a 25 A cutoff gives stable results
                //with alpha=0, and this way we avoid changing the hamiltonian significantly
                //by having an erfc function
                for (int j=0; j<nats1; ++j)
                {
                    const Parameter &param1 = params1_array[j];
                    const Vector &c1 = coords1_array[j];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];
                        const Vector &c0 = coords0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += (param0.reduced_charge * param1.reduced_charge) *
                                        (one_over_r - one_over_Rc + one_over_Rc2 * (r-Rc));
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }
            }
            #endif
        }
        else if (useReactionField())
        {
            //use the reaction field potential
            const double Rc = coul_cutoff;
            const double Rlj = lj_cutoff;

            const double rf_dielectric = reactionFieldDielectric();
            const double k_rf = (1.0 / pow_3(Rc)) * ( (rf_dielectric-1) /
                                                      (2*rf_dielectric + 1) );
            const double c_rf = (1.0 / Rc) * ( (3*rf_dielectric) /
                                               (2*rf_dielectric + 1) );

            #ifdef SIRE_USE_SSE
            {
                const int remainder = nats1 % 2;

                __m128d sse_cnrg = { 0, 0 };
                __m128d sse_ljnrg = { 0, 0 };

                const __m128d sse_one = { 1.0, 1.0 };

                const __m128d sse_Rc = { Rc, Rc };
                const __m128d sse_Rlj = { Rlj, Rlj };

                const __m128d sse_k_rf = { k_rf, k_rf };
                const __m128d sse_c_rf = { c_rf, c_rf };

                for (int j=0; j<nats1-1; j += 2)
                {
                    const Parameter &param10 = params1_array[j];
                    const Parameter &param11 = params1_array[j+1];

                    const Vector &c10 = coords1_array[j];
                    const Vector &c11 = coords1_array[j+1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];

                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge,
                                                       param0.reduced_charge );

                        const Vector &c0 = coords0_array[i];

                        __m128d sse_r;
                        {
                            __m128d dx = _mm_set_pd( c0.x() - c10.x(), c0.x() - c11.x() );
                            __m128d dy = _mm_set_pd( c0.y() - c10.y(), c0.y() - c11.y() );
                            __m128d dz = _mm_set_pd( c0.z() - c10.z(), c0.z() - c11.z() );

                            __m128d dx2 = _mm_mul_pd(dx, dx);
                            __m128d dy2 = _mm_mul_pd(dy, dy);
                            __m128d dz2 = _mm_mul_pd(dz, dz);

                            __m128d sse_r2 = _mm_add_pd(dx2,dy2);
                            sse_r2 = _mm_add_pd(sse_r2, dz2);

                            sse_r = _mm_sqrt_pd(sse_r2);
                        }

                        const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                        //calculate the coulomb energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);

                            __m128d nrg = _mm_mul_pd(sse_r, sse_r);
                            nrg = _mm_mul_pd(nrg, sse_k_rf);
                            nrg = _mm_sub_pd(nrg, sse_c_rf);
                            nrg = _mm_add_pd(nrg, sse_one_over_r);

                            __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                          param11.reduced_charge );

                            sse_chg = _mm_mul_pd(sse_chg, sse_chg0);

                            nrg = _mm_mul_pd(sse_chg, nrg);

                            nrg = _mm_and_pd(nrg, sse_in_cutoff);

                            sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                        }

                        //calculate the LJ energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);
                            const LJPair &ljpair0 = ljpairs.constData()[
                                                ljpairs.map(param0.ljid,
                                                            param10.ljid)];

                            const LJPair &ljpair1 = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param11.ljid)];

                            const __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                            const __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                                ljpair1.epsilon() );

                            //calculate (sigma/r)^6 and (sigma/r)^12
                            __m128d sse_sig_over_r2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                            sse_sig_over_r2 = _mm_mul_pd( sse_sig_over_r2,
                                                          sse_sig_over_r2 );

                            __m128d sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r2,
                                                                 sse_sig_over_r2);

                            sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r6,
                                                         sse_sig_over_r2);

                            __m128d sse_sig_over_r12 = _mm_mul_pd(sse_sig_over_r6,
                                                                  sse_sig_over_r6);

                            //calculate LJ energy (the factor of 4 is added later)
                            __m128d tmp = _mm_sub_pd(sse_sig_over_r12,
                                                     sse_sig_over_r6);

                            tmp = _mm_mul_pd(tmp, sse_eps);
                            tmp = _mm_and_pd(tmp, sse_in_cutoff);
                            sse_ljnrg = _mm_add_pd(sse_ljnrg, tmp);
                        }
                    }
                }

                if (remainder == 1)
                {
                    const Vector &c1 = coords1_array[nats1-1];
                    const Parameter &param1 = params1_array[nats1-1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Vector &c0 = coords0_array[i];
                        const Parameter &param0 = params0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += param0.reduced_charge * param1.reduced_charge
                                        * (one_over_r + k_rf*r*r - c_rf);
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }

                icnrg += *((const double*)&sse_cnrg) +
                         *( ((const double*)&sse_cnrg) + 1 );

                iljnrg += *((const double*)&sse_ljnrg) +
                          *( ((const double*)&sse_ljnrg) + 1 );
            }
            #else
            {
                //we use the force-shifted cutoff potential described in equation 18
                //in Fennell and Gezelter, J. Chem. Phys., 124, 234104, 2006
                //We use alpha=0, as I have seen that a 25 A cutoff gives stable results
                //with alpha=0, and this way we avoid changing the hamiltonian significantly
                //by having an erfc function
                for (int j=0; j<nats1; ++j)
                {
                    const Parameter &param1 = params1_array[j];
                    const Vector &c1 = coords1_array[j];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];
                        const Vector &c0 = coords0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += (param0.reduced_charge * param1.reduced_charge) *
                                     (one_over_r + k_rf*r*r - c_rf);
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }
            }
            #endif
        }
        else
        {
            //use a simple atom-based cutoff
            const double Rc = coul_cutoff;
            const double Rlj = lj_cutoff;

            #ifdef SIRE_USE_SSE
            {
                const int remainder = nats1 % 2;

                __m128d sse_cnrg = { 0, 0 };
                __m128d sse_ljnrg = { 0, 0 };

                const __m128d sse_one = { 1.0, 1.0 };
                const __m128d sse_Rc = { Rc, Rc };
                const __m128d sse_Rlj = { Rlj, Rlj };

                for (int j=0; j<nats1-1; j += 2)
                {
                    const Parameter &param10 = params1_array[j];
                    const Parameter &param11 = params1_array[j+1];

                    const Vector &c10 = coords1_array[j];
                    const Vector &c11 = coords1_array[j+1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];

                        __m128d sse_chg0 = _mm_set_pd( param0.reduced_charge,
                                                       param0.reduced_charge );

                        const Vector &c0 = coords0_array[i];

                        __m128d sse_r;
                        {
                            __m128d dx = _mm_set_pd( c0.x() - c10.x(), c0.x() - c11.x() );
                            __m128d dy = _mm_set_pd( c0.y() - c10.y(), c0.y() - c11.y() );
                            __m128d dz = _mm_set_pd( c0.z() - c10.z(), c0.z() - c11.z() );

                            __m128d dx2 = _mm_mul_pd(dx, dx);
                            __m128d dy2 = _mm_mul_pd(dy, dy);
                            __m128d dz2 = _mm_mul_pd(dz, dz);

                            __m128d sse_r2 = _mm_add_pd(dx2,dy2);
                            sse_r2 = _mm_add_pd(sse_r2, dz2);

                            sse_r = _mm_sqrt_pd(sse_r2);
                        }

                        const __m128d sse_one_over_r = _mm_div_pd(sse_one, sse_r);

                        //calculate the coulomb energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rc);

                            __m128d sse_chg = _mm_set_pd( param10.reduced_charge,
                                                          param11.reduced_charge );

                            sse_chg = _mm_mul_pd(sse_chg, sse_chg0);

                            __m128d nrg = _mm_mul_pd(sse_chg, sse_one_over_r);

                            nrg = _mm_and_pd(nrg, sse_in_cutoff);

                            sse_cnrg = _mm_add_pd(sse_cnrg, nrg);
                        }

                        //calculate the LJ energy
                        {
                            const __m128d sse_in_cutoff = _mm_cmplt_pd(sse_r, sse_Rlj);

                            const LJPair &ljpair0 = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param10.ljid)];

                            const LJPair &ljpair1 = ljpairs.constData()[
                                                        ljpairs.map(param0.ljid,
                                                                    param11.ljid)];

                            const __m128d sse_sig = _mm_set_pd( ljpair0.sigma(), ljpair1.sigma() );
                            const __m128d sse_eps = _mm_set_pd( ljpair0.epsilon(),
                                                                ljpair1.epsilon() );

                            //calculate (sigma/r)^6 and (sigma/r)^12
                            __m128d sse_sig_over_r2 = _mm_mul_pd(sse_sig, sse_one_over_r);
                            sse_sig_over_r2 = _mm_mul_pd( sse_sig_over_r2,
                                                          sse_sig_over_r2 );

                            __m128d sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r2,
                                                                 sse_sig_over_r2);

                            sse_sig_over_r6 = _mm_mul_pd(sse_sig_over_r6,
                                                         sse_sig_over_r2);

                            __m128d sse_sig_over_r12 = _mm_mul_pd(sse_sig_over_r6,
                                                                  sse_sig_over_r6);

                            //calculate LJ energy (the factor of 4 is added later)
                            __m128d tmp = _mm_sub_pd(sse_sig_over_r12,
                                                     sse_sig_over_r6);

                            tmp = _mm_mul_pd(tmp, sse_eps);
                            tmp = _mm_and_pd(tmp, sse_in_cutoff);
                            sse_ljnrg = _mm_add_pd(sse_ljnrg, tmp);
                        }
                    }
                }

                if (remainder == 1)
                {
                    const Vector &c1 = coords1_array[nats1-1];
                    const Parameter &param1 = params1_array[nats1-1];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Vector &c0 = coords0_array[i];
                        const Parameter &param0 = params0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += param0.reduced_charge * param1.reduced_charge
                                                * one_over_r;
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }

                icnrg += *((const double*)&sse_cnrg) +
                         *( ((const double*)&sse_cnrg) + 1 );

                iljnrg += *((const double*)&sse_ljnrg) +
                          *( ((const double*)&sse_ljnrg) + 1 );
            }
            #else
            {
                // use an atomistic cutoff
                for (int j=0; j<nats1; ++j)
                {
                    const Parameter &param1 = params1_array[j];
                    const Vector &c1 = coords1_array[j];

                    for (int i=0; i<nats0; ++i)
                    {
                        const Parameter &param0 = params0_array[i];
                        const Vector &c0 = coords0_array[i];

                        const double r = Vector::distance(c0,c1);
                        const double one_over_r = double(1) / r;

                        if (r < Rc)
                        {
                            icnrg += (param0.reduced_charge * param1.reduced_charge) *
                                        one_over_r;
                        }

                        if (r < Rlj)
                        {
                            const LJPair &ljpair = ljpairs.constData()[
                                                    ljpairs.map(param0.ljid,
                                                                param1.ljid)];

                            double sig_over_dist6 = pow_6(ljpair.sigma()*one_over_r);
                            double sig_over_dist12 = pow_2(sig_over_dist6);

                            iljnrg += ljpair.epsilon() * (sig_over_dist12 -
                                                          sig_over_dist6);
                        }
                    }
                }
            }
            #endif
        }
    }

    double gridnrg = 0;
    const double *gridpot_array = gridpot.constData();

    //now calculate the energy in the grid
    for (int i0=0; i0<nats0; ++i0)
    {
        const Vector &c0 = coords0[i0];
        const detail::CLJParameter &p0 = params0[i0];

        Vector grid_coords = c0 - gridbox.minCoords();

        int i_0 = int(grid_coords.x() / grid_spacing);
        int j_0 = int(grid_coords.y() / grid_spacing);
        int k_0 = int(grid_coords.z() / grid_spacing);

        if (i_0 < 0 or i_0 >= int(dimx-1) or
            j_0 < 0 or j_0 >= int(dimy-1) or
            k_0 < 0 or k_0 >= int(dimz-1))
        {
            qDebug() << "POINT" << c0.toString() << "LIES OUTSIDE OF "
                     << "THE GRID?" << gridbox.toString();
        }
        else
        {
            //use tri-linear interpolation to get the potential at the atom
            //
            // This is described in
            //
            // Davis, Madura and McCammon, Comp. Phys. Comm., 62, 187-197, 1991
            //
            // phi(x,y,z) = phi(i  ,j  ,k  )*(1-R)(1-S)(1-T) +
            //              phi(i+1,j  ,k  )*(  R)(1-S)(1-T) +
            //              phi(i  ,j+1,k  )*(1-R)(  S)(1-T) +
            //              phi(i  ,j  ,k+1)*(1-R)(1-S)(  T) +
            //              phi(i+1,j+1,k  )*(  R)(  S)(1-T) +
            //              phi(i+1,j  ,k+1)*(  R)(1-S)(  T) +
            //              phi(i  ,j+1,k+1)*(1-R)(  S)(  T) +
            //              phi(i+1,j+1,k+1)*(  R)(  S)(  T) +
            //
            // where R, S and T are the coordinates of the atom in
            // fractional grid coordinates from the point (i,j,k), e.g.
            // (0,0,0) is (i,j,k) and (1,1,1) is (i+1,j+1,k+1)
            //
            const Vector c000 = gridbox.minCoords() +
                                    Vector( i_0 * grid_spacing,
                                            j_0 * grid_spacing,
                                            k_0 * grid_spacing );

            const Vector RST = (c0 - c000) / grid_spacing;
            const double R = RST.x();
            const double S = RST.y();
            const double T = RST.z();

            int i000 = (i_0  ) * (dimy*dimz) + (j_0  )*dimz + (k_0  );
            int i001 = (i_0  ) * (dimy*dimz) + (j_0  )*dimz + (k_0+1);
            int i010 = (i_0  ) * (dimy*dimz) + (j_0+1)*dimz + (k_0  );
            int i100 = (i_0+1) * (dimy*dimz) + (j_0  )*dimz + (k_0  );
            int i011 = (i_0  ) * (dimy*dimz) + (j_0+1)*dimz + (k_0+1);
            int i101 = (i_0+1) * (dimy*dimz) + (j_0  )*dimz + (k_0+1);
            int i110 = (i_0+1) * (dimy*dimz) + (j_0+1)*dimz + (k_0  );
            int i111 = (i_0+1) * (dimy*dimz) + (j_0+1)*dimz + (k_0+1);

            double phi = (gridpot_array[i000] * (1-R)*(1-S)*(1-T)) +
                         (gridpot_array[i001] * (1-R)*(1-S)*(  T)) +
                         (gridpot_array[i010] * (1-R)*(  S)*(1-T)) +
                         (gridpot_array[i100] * (  R)*(1-S)*(1-T)) +
                         (gridpot_array[i011] * (1-R)*(  S)*(  T)) +
                         (gridpot_array[i101] * (  R)*(1-S)*(  T)) +
                         (gridpot_array[i110] * (  R)*(  S)*(1-T)) +
                         (gridpot_array[i111] * (  R)*(  S)*(  T));

            gridnrg += phi * p0.reduced_charge;
        }
    }

    cnrg = icnrg + gridnrg;
    ljnrg = 4.0*iljnrg;  // 4 epsilon (....)
}

/** Ensure that the next energy evaluation is from scratch */
void GridFF::mustNowRecalculateFromScratch()
{
    gridpot.clear();
    closemols_coords.clear();
    closemols_params.clear();
    oldnrgs.clear();

    InterGroupCLJFF::mustNowRecalculateFromScratch();
}

/** Any additions mean that the forcefield must be recalculated from scratch */
void GridFF::_pvt_added(quint32 groupid, const PartialMolecule &mol,
                        const PropertyMap &map)
{
    this->mustNowRecalculateFromScratch();
    InterGroupCLJFF::_pvt_added(groupid, mol, map);
}

/** Any removals mean that the forcefield must be recalculated from scratch */
void GridFF::_pvt_removed(quint32 groupid, const PartialMolecule &mol)
{
    this->mustNowRecalculateFromScratch();
    InterGroupCLJFF::_pvt_removed(groupid, mol);
}

/** Any changes to group 1 mean that the forcefield must be recalculated from scratch */
void GridFF::_pvt_changed(quint32 groupid, const SireMol::Molecule &molecule, bool auto_commit)
{
    InterGroupCLJFF::_pvt_changed(groupid, molecule, auto_commit);
}

/** Any changes to group 1 mean that the forcefield must be recalculated from scratch */
void GridFF::_pvt_changed(quint32 groupid, const QList<SireMol::Molecule> &molecules,
                          bool auto_commit)
{
    InterGroupCLJFF::_pvt_changed(groupid, molecules, auto_commit);
}

/** Any removals mean that the forcefield must be recalculated from scratch */
void GridFF::_pvt_removedAll(quint32 groupid)
{
    this->mustNowRecalculateFromScratch();
    InterGroupCLJFF::_pvt_removedAll(groupid);
}

/** Recalculate the total energy */
void GridFF::recalculateEnergy()
{
    CLJPotential::startEvaluation();

    if (mols[0].isEmpty() or (mols[1].isEmpty() and fixedatoms_coords.isEmpty()))
    {
        //one of the two groups is empty, so the energy must be zero
        this->components().setEnergy(*this, CLJEnergy(0,0));
        CLJPotential::finishedEvaluation();
        this->setClean();
        return;
    }

    bool must_recalculate = false;

    //we need to recalculate if either something has changed in group 1,
    //or nothing at all has changed (in which case we assume that changes
    //are not being recorded), or if the grid hasn't been created,
    //or if molecules have been added or removed to group 0.
    if (gridpot.isEmpty())
    {
        //the grid is empty
        must_recalculate = true;
    }
    else if (not changed_mols[1].isEmpty())
    {
        //check to see if the parts in this forcefield have changed
        for (QHash<MolNum,ChangedMolecule>::const_iterator
                                      it = this->changed_mols[1].constBegin();
             it != this->changed_mols[1].constEnd();
             ++it)
        {
            if (not (it->newParts().isEmpty() and it->oldParts().isEmpty()))
            {
                //a part of this molecule in the forcefield has changed
                must_recalculate = true;
                break;
            }
        }

        if (not must_recalculate)
        {
            changed_mols[1].clear();

            if (changed_mols[0].isEmpty())
            {
                //there is nothing to do :-)
                changed_mols[1].clear();
                this->setClean();
                return;
            }

            //otherwise a part of the molecule in group 0 has changed,
            //so evaluate the change
        }
    }
    else if (changed_mols[0].isEmpty())
    {
        //probably not recording changes - assume everything has changed
        must_recalculate = true;
    }

    if (must_recalculate)
    {
        QElapsedTimer t;
        t.start();

        this->mustNowRecalculateFromScratch();
        this->rebuildGrid();

        qint64 ns = t.nsecsElapsed();
        t.restart();

        double total_cnrg(0);
        double total_ljnrg(0);

        //calculate all of the energies from scratch and store
        //them in 'oldnrgs'
        oldnrgs.reserve( mols[0].count() );

        //loop through all of the molecules and calculate the energies.
        //First, calculate the CLJ energies for the closemols
        for (ChunkedVector<CLJMolecule>::const_iterator
                    it = mols[0].moleculesByIndex().constBegin();
             it != mols[0].moleculesByIndex().constEnd();
             ++it)
        {
            const CLJMolecule &cljmol = *it;

            double cnrg(0);
            double ljnrg(0);

            //loop through each CutGroup of this molecule
            const int ngroups = cljmol.coordinates().count();

            const CoordGroup *groups_array = cljmol.coordinates().constData();

            const CLJParameters::Array *params_array
                                = cljmol.parameters().atomicParameters().constData();

            for (int igroup=0; igroup<ngroups; ++igroup)
            {
                const CoordGroup &group = groups_array[igroup];

                if (not gridbox.contains(group.aaBox()))
                {
                    //this group lies outside the grid - we need to recalculate
                    //the grid
                    this->mustNowRecalculateFromScratch();
                    this->recalculateEnergy();
                    return;
                }

                double icnrg, iljnrg;
                calculateEnergy(group, params_array[igroup],
                                icnrg, iljnrg);

                cnrg += icnrg;
                ljnrg += iljnrg;

            }

            oldnrgs.insert(cljmol.number(), CLJEnergy(cnrg,ljnrg));

            total_cnrg += cnrg;
            total_ljnrg += ljnrg;
        }

        ns = t.nsecsElapsed();

        this->components().setEnergy(*this, CLJEnergy(total_cnrg,total_ljnrg));
    }
    else
    {
        double delta_cnrg(0);
        double delta_ljnrg(0);

        for (QHash<MolNum,ChangedMolecule>::const_iterator
                                      it = this->changed_mols[0].constBegin();
             it != this->changed_mols[0].constEnd();
             ++it)
        {
            if (it->nothingChanged())
            {
                continue;
            }
            else if (it->changedAll())
            {
                const CLJMolecule &cljmol = it->newMolecule();

                //loop through each CutGroup of this molecule
                const int ngroups = cljmol.coordinates().count();

                const CoordGroup *groups_array = cljmol.coordinates().constData();

                const CLJParameters::Array *params_array
                                    = cljmol.parameters().atomicParameters().constData();

                double cnrg(0);
                double ljnrg(0);

                for (int igroup=0; igroup<ngroups; ++igroup)
                {
                    const CoordGroup &group = groups_array[igroup];

                    if (not gridbox.contains(group.aaBox()))
                    {
                        //this group lies outside the grid - we need to recalculate
                        //the grid
                        this->mustNowRecalculateFromScratch();
                        this->recalculateEnergy();
                        return;
                    }

                    double icnrg, iljnrg;
                    calculateEnergy(group, params_array[igroup],
                                    icnrg, iljnrg);

                    cnrg += icnrg;
                    ljnrg += iljnrg;

                }

                CLJEnergy old_nrg = oldnrgs[cljmol.number()];
                oldnrgs[cljmol.number()] = CLJEnergy(cnrg,ljnrg);

                delta_cnrg += (cnrg - old_nrg.coulomb());
                delta_ljnrg += (ljnrg - old_nrg.lj());
            }
            else
            {
                //only calculate the change in energy for the part of the
                //molecule that has changed
                const CLJMolecule &oldmol = it->oldParts();
                const CLJMolecule &newmol = it->newParts();

                //calculate the energy of the old parts
                double old_cnrg(0);
                double old_ljnrg(0);
                {
                    const int ngroups = oldmol.coordinates().count();

                    const CoordGroup *groups_array = oldmol.coordinates().constData();

                    const CLJParameters::Array *params_array
                                    = oldmol.parameters().atomicParameters().constData();

                    for (int igroup=0; igroup<ngroups; ++igroup)
                    {
                        const CoordGroup &group = groups_array[igroup];

                        if (not gridbox.contains(group.aaBox()))
                        {
                            //this group lies outside the grid - we need to recalculate
                            //the grid
                            this->mustNowRecalculateFromScratch();
                            this->recalculateEnergy();
                            return;
                        }

                        double icnrg, iljnrg;
                        calculateEnergy(group, params_array[igroup],
                                        icnrg, iljnrg);

                        old_cnrg += icnrg;
                        old_ljnrg += iljnrg;
                    }
                }

                //calculate the energy of the new parts
                double new_cnrg(0);
                double new_ljnrg(0);
                {
                    const int ngroups = newmol.coordinates().count();

                    const CoordGroup *groups_array = newmol.coordinates().constData();

                    const CLJParameters::Array *params_array
                                    = newmol.parameters().atomicParameters().constData();

                    for (int igroup=0; igroup<ngroups; ++igroup)
                    {
                        const CoordGroup &group = groups_array[igroup];

                        if (not gridbox.contains(group.aaBox()))
                        {
                            //this group lies outside the grid - we need to recalculate
                            //the grid
                            this->mustNowRecalculateFromScratch();
                            this->recalculateEnergy();
                            return;
                        }

                        double icnrg, iljnrg;
                        calculateEnergy(group, params_array[igroup],
                                        icnrg, iljnrg);

                        new_cnrg += icnrg;
                        new_ljnrg += iljnrg;
                    }
                }

                CLJEnergy old_nrg = oldnrgs[oldmol.number()];
                CLJEnergy new_nrg = CLJEnergy(old_nrg.coulomb() + new_cnrg - old_cnrg,
                                              old_nrg.lj() + new_ljnrg - old_ljnrg);

                oldnrgs[oldmol.number()] = new_nrg;

                delta_cnrg += (new_nrg.coulomb() - old_nrg.coulomb());
                delta_ljnrg += (new_nrg.lj() - old_nrg.lj());
            }
        }

        //change the energy
        this->components().changeEnergy(*this, CLJEnergy(delta_cnrg,delta_ljnrg));

        //clear the changed molecules
        this->changed_mols[0].clear();
    }

    CLJPotential::finishedEvaluation();
    this->setClean();
}
