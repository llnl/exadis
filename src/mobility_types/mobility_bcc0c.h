/*---------------------------------------------------------------------------
 *
 *	ExaDiS
 *
 *	Nicolas Bertin
 *	bertin1@llnl.gov
 *
 *-------------------------------------------------------------------------*/

#pragma once
#ifndef EXADIS_MOBILITY_BCC0C_H
#define EXADIS_MOBILITY_BCC0C_H

#include "mobility.h"

namespace ExaDiS {

/*---------------------------------------------------------------------------
 *
 *    Class:        MobilityBCC0c
 *                  WARNING: This is an experimental non-linear mobility
 *                  law for BCC materials that implements thermally-activated
 *                  screw mobility and non-linear edge mobility.
 *                  This mobility has NOT been thouroughly tested and may
 *                  result in unphysical behaviors. It is provided as an
 *                  example only; do not use in production runs unless
 *                  fully tested and validated first.
 *
 *-------------------------------------------------------------------------*/
struct MobilityBCC0c
{
    bool non_linear = true;
    double Beclimbj, Beclimb2;
    double Bline, BlmBecl;
    double invBedge2;
    double Bedge, Bscrew0, Beclimb;
    double vscale;
    
    struct Params {
        double Medge = -1.0, Mscrew = -1.0, Mclimb = -1.0;
        double FscrewTATasym = 0.0;
        double vmax = -1.0;
        double tempK = 300.0;
        double tau0 = 0.25e9; // Peierls stress, Pa
        double p = 0.75, q = 1.1; // Enthalpy exponents, unitless
        double H0 = 1.06; // Enthalpy pre-factor, eV
        double T0 = -1.0; // Transition temperature, K
        double tauv = 0.0; // Shift stress, Pa
        double B1screw = 0.0; // Screw temperature dependant drag
        double B1edge = 0.0; // Edge temperature dependant drag
        double Dedge = 0.0; // Edge non-linear drag
        double D1edge = 0.0; // Edge temperature dependant non-linear drag
        double v0edge = -1.0; // Edge non-linear velocity, m/s
        Params() = default;
        Params(double _Medge, double _Mscrew, double _Mclimb, double _vmax=-1.0) {
            Medge = _Medge;
            Mscrew = _Mscrew;
            Mclimb = _Mclimb;
            vmax = _vmax;
        }
        Params(Dict paramslist) {
            for (auto const& [key, val] : paramslist) {
                std::string name = dict::get_key(key);
                if      (name == "Medge") Medge = dict::get_val<double>(val);
                else if (name == "Mscrew") Mscrew = dict::get_val<double>(val);
                else if (name == "Mclimb") Mclimb = dict::get_val<double>(val);
                else if (name == "FscrewTATasym") FscrewTATasym = dict::get_val<double>(val);
                else if (name == "vmax") vmax = dict::get_val<double>(val);
                else if (name == "tempK") tempK = dict::get_val<double>(val);
                else if (name == "tau0") tau0 = dict::get_val<double>(val);
                else if (name == "p") p = dict::get_val<double>(val);
                else if (name == "q") q = dict::get_val<double>(val);
                else if (name == "H0") H0 = dict::get_val<double>(val);
                else if (name == "T0") T0 = dict::get_val<double>(val);
                else if (name == "tauv") tauv = dict::get_val<double>(val);
                else if (name == "B1screw") B1screw = dict::get_val<double>(val);
                else if (name == "B1edge") B1edge = dict::get_val<double>(val);
                else if (name == "Dedge") Dedge = dict::get_val<double>(val);
                else if (name == "D1edge") D1edge = dict::get_val<double>(val);
                else if (name == "v0edge") v0edge = dict::get_val<double>(val);
                else ExaDiS_fatal("Error: unknown MobilityBCC0c input parameter '%s'\n", name.c_str());
            }
        }
    };
    Params params;
    
    MobilityBCC0c(System* system, Params& _params)
    {
        if (system->crystal.type != BCC_CRYSTAL)
            ExaDiS_fatal("Error: MobilityBCC0c() must be used with BCC crystal type\n");
        
        params = _params;
        if (params.Medge < 0 || params.Mscrew < 0.0 || params.Mclimb < 0.0)
            ExaDiS_fatal("Error: invalid or missing MobilityBCC0c input parameter values\n");
        
        Bedge   = 1.0 / params.Medge + params.B1edge * params.tempK;
        Bscrew0 = 1.0 / params.Mscrew + params.B1screw * params.tempK;
        Beclimb = 1.0 / params.Mclimb;
        vscale  = system->params.burgmag; //vscale (convert factor from m/s)
        
        // Initialization
        Beclimbj   = Beclimb;
        Beclimb2   = Beclimb * Beclimb;
        Bline      = 1.0e-2 * MIN(Bscrew0, Bedge);
        BlmBecl    = Bline - Beclimbj;
        invBedge2  = 1.0 / (Bedge*Bedge);
    }

    KOKKOS_INLINE_FUNCTION
    double chi_angle(System* system, const Vec3& burg, const Vec3& linedir, const Vec3& fi)
    {
        Vec3 bcryst = system->crystal.Rinv * burg;
        Vec3 n0 = Vec3(2*bcryst.x, -bcryst.y, -bcryst.z).normalized();
        n0 = system->crystal.R * n0; // back to lab frame
        Vec3 t0 = cross(n0, linedir).normalized();
        double chi = atan2(dot(fi, n0), dot(fi, t0));
        chi = -M_PI/6.0*cos(3.0*chi); // fold to T/AT zones
        return chi;
    }
    
    template<class N>
    KOKKOS_INLINE_FUNCTION
    Vec3 node_velocity(System* system, N* net, const int& i, const Vec3& fi)
    {
        auto nodes = net->get_nodes();
        auto segs = net->get_segs();
        auto conn = net->get_conn();
        auto cell = net->cell;
        
        Vec3 vi(0.0);
        
        if (conn[i].num >= 2 && nodes[i].constraint != PINNED_NODE) {
            
            int linejunc = 0;
            Vec3 tjunc(0.0);
            
            if (system->params.split3node) {
                
                int binaryjunc = 0;
                int planarjunc = 0;
                binaryjunc = BCC_binary_junction_node(system, net, i, tjunc, &planarjunc);                
                binaryjunc = (binaryjunc > -1);

                if (binaryjunc) {
                    if (planarjunc) {
                        linejunc = 1;
                    } else {
                        int unzipping = (dot(fi, tjunc) > 0.0);
                        // If the node is unzipping the junction and we are
                        // not treating unzipping as a purely topological
                        // operation, then orthogonalize the climb direction
                        if (unzipping) linejunc = 1;
                    }
                }
            }

            double eps = 1e-12;

            double totblen = 0.0;
            Vec3 vdir(0.0);
            int numNonZeroLenSegs = 0;
            Vec3 r1 = nodes[i].pos;

            for (int ipass = 0; ipass < 2; ipass++) {

                // Build drag matrix
                Mat33 Btotal = Mat33().zero();
                double FricForce = 0.0;
                for (int j = 0; j < conn[i].num; j++) {

                    int k = conn[i].node[j];
                    int s = conn[i].seg[j];
                    int order = conn[i].order[j];

                    Vec3 burg = order*segs[s].burg;
                    double bMag = burg.norm();
                    double bMag2 = bMag*bMag;
                    double invbMag2 = 1.0 / bMag2;

                    Vec3 r2 = cell.pbc_position(r1, nodes[k].pos);
                    
                    Vec3 dr = r2-r1;
                    double mag = dr.norm();
                    if (mag < eps) continue;
                    if (ipass == 0) {
                        numNonZeroLenSegs++;
                        totblen += mag * bMag;
                    }

                    double halfMag = 0.5 * mag;
                    double invMag  = 1.0 / mag;
                    dr = invMag * dr;

                    double costheta = dot(dr, burg);
                    double costheta2 = (costheta*costheta) * invbMag2;
                    
                    if (bMag > 1.0+eps) {
                        // [0 0 1] arms don't move as readily as other arms, so must be
                        // handled specially.
                        
                        // The junction node move along the junction line freely
                        // No drag on the junction line
                        if (linejunc == 1) {
                            Btotal[0][0] += halfMag * Bline;
                            Btotal[1][1] += halfMag * Bline;
                            Btotal[2][2] += halfMag * Bline;
                            continue;
                        }
                        
                        if (conn[i].num == 2) {
                            Btotal[0][0] += halfMag * Beclimbj;
                            Btotal[1][1] += halfMag * Beclimbj;
                            Btotal[2][2] += halfMag * Beclimbj;
                        } else {
                            Btotal[0][0] += halfMag * (dr.x*dr.x * BlmBecl + Beclimbj);
                            Btotal[0][1] += halfMag * (dr.x*dr.y * BlmBecl);
                            Btotal[0][2] += halfMag * (dr.x*dr.z * BlmBecl);
                            Btotal[1][1] += halfMag * (dr.y*dr.y * BlmBecl + Beclimbj);
                            Btotal[1][2] += halfMag * (dr.y*dr.z * BlmBecl);
                            Btotal[2][2] += halfMag * (dr.z*dr.z * BlmBecl + Beclimbj);
                        }
                    } else  {
                        // Arm is not [0 0 1], so build the drag matrix assuming the
                        // dislocation is screw type

                        double Bscrew, tau, fricStress;
                        if (ipass == 0) {
                            // Use linear Bscrew coefficient for first pass
                            // This is to estimate the velocity direction
                            Bscrew = Bscrew0;
                            tau = 0.0;
                            fricStress = 0.0;
                        } else {
                            //tau = 2.0 * fabs(dot(fi, vdir)) / totblen; // Pa
                            tau = 2.0 * fi.norm() / totblen; // Pa
                            
                            double tdotb = dot(dr, 1.0 / bMag * burg);
                            Vec3 linedir = SIGN(tdotb) / bMag * burg;
                            double chi = chi_angle(system, burg, linedir, fi);
                            double taup = params.tau0 + params.FscrewTATasym * params.tau0 * (chi+M_PI/6.0)/M_PI*3.0;

                            double dangle = 1.0 / bMag * fabs(costheta);
                            double tangle = dangle + eps;

                            double dangle_thres1 = cos(0.0 * M_PI/180.0);
                            double dangle_thres2 = cos(5.0 * M_PI/180.0);
                            if (dangle >= dangle_thres1) {
                                tangle = dangle;
                            } else if (dangle >= dangle_thres2) {
                                tangle = dangle*(dangle-dangle_thres2)/(dangle_thres1-dangle_thres2);
                            } else {
                                tangle = eps;
                            }

                            taup = tangle * taup;

                        #if 1
                            double PT = 1.0-pow(tau/taup, params.p);
                            double Tratio = (params.T0 > 0.0) ? params.tempK/params.T0 : 0.0;
                            double dGT = params.H0 * (SIGN(PT) * pow(fabs(PT), params.q) - Tratio);

                            double a = 2.0/sqrt(3.0);
                            double h = a * sqrt(2.0/3.0);
                            double w = 25.0; // reference line length
                            double Bkink = Bscrew0 / 10.0;

                            double kbT = 8.6173e-5 * params.tempK; // eV
                            double dexpGT = exp(dGT/2.0/kbT);
                            double B_kink = a * (2.0*a*dexpGT+w) / (2.0*h*w) * Bkink * dexpGT;

                            double tauvn = (params.tauv < 0.0) ? -params.tauv : 0.0;
                            double tauvd = (params.tauv > 0.0) ? params.tauv : 0.0;
                            double B_lin = (tau > 0.0) ? (tau+tauvn*dangle) / (tau+tauvd*dangle) * Bscrew0 : Bscrew0;

                            Bscrew = 0.2*B_kink + B_lin;

                        #else
                            // Test friction stress
                            double v = fmax(tau / Bscrew0 / 1e5, (tau-taup) / Bscrew0);
                            if (v > 0.0)
                                Bscrew = tau / v;
                            else
                                Bscrew = Bscrew0;
                        #endif

                            double vs = tau / Bscrew;
                            fricStress = tau - Bscrew0*vs;
                            //FricForce += fricStress * bMag * mag;
                            Bscrew = Bscrew0;
                        }

                        double Bscrew2 = Bscrew * Bscrew;
                        double invBscrew2 = 1.0 / Bscrew2;
                        double BlmBsc = Bline - Bscrew;

                        Btotal[0][0] += halfMag * (dr.x*dr.x * BlmBsc + Bscrew);
                        Btotal[0][1] += halfMag * (dr.x*dr.y * BlmBsc);
                        Btotal[0][2] += halfMag * (dr.x*dr.z * BlmBsc);
                        Btotal[1][1] += halfMag * (dr.y*dr.y * BlmBsc + Bscrew);
                        Btotal[1][2] += halfMag * (dr.y*dr.z * BlmBsc);
                        Btotal[2][2] += halfMag * (dr.z*dr.z * BlmBsc + Bscrew);

                        // Now correct the drag matrix for dislocations that are
                        // not screw
                        if ((1.0 - costheta2) > eps) {

                            double invsqrt1mcostheta2 = 1.0 / sqrt((1.0 - costheta2) * bMag2);
                            Vec3 nr = cross(burg, dr);
                            nr = invsqrt1mcostheta2 * nr;

                            Vec3 mr = cross(nr, dr);
                            
                            // Orthogonalize climb direction wrt junction line direction
                            // to avoid numerical issues with binary junction nodes
                            if (linejunc == 1) {
                                nr -= dot(nr, tjunc) * tjunc;
                                nr = nr.normalized();
                            }

                            if (ipass == 1 && params.v0edge > 0.0 && params.Dedge > 0.0) {
                                double v0 = params.v0edge / vscale;
                                double AT = Bedge;
                                double D = MAX(params.Dedge - params.D1edge * params.tempK, 0.0);
                                double ve = tau / AT;
                                if (ve > v0) {
                                    double v = v0 + (tau - AT*v0) / AT;
                                    for (int iter = 0; iter < 10; iter++) {
                                        ve = v - (AT*v + D*pow(v-v0, 1.5) - tau) / (AT + 1.5*D*pow(v-v0, 0.5));
                                        if (fabs(ve - v) / v < 1e-5) break;
                                        v = ve;
                                    }
                                    double fricStressEdge = tau - Bedge*ve;
                                    if (fabs(fricStressEdge) > fabs(fricStress))
                                        fricStress = fricStressEdge;
                                    //FricForce += fricStress * bMag * mag;
                                }
                            }
                            
                            double Bglide = sqrt(invBedge2+(invBscrew2-invBedge2)*costheta2);
                            Bglide = 1.0 / Bglide;
                            double Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) * costheta2);
                            double BclmBsc = Bclimb - Bscrew;
                            double BglmBsc = Bglide - Bscrew;

                            Btotal[0][0] += halfMag * (nr.x*nr.x * BclmBsc +
                                                    mr.x*mr.x * BglmBsc);
                            Btotal[0][1] += halfMag * (nr.x*nr.y * BclmBsc +
                                                    mr.x*mr.y * BglmBsc);
                            Btotal[0][2] += halfMag * (nr.x*nr.z * BclmBsc +
                                                    mr.x*mr.z * BglmBsc);
                            Btotal[1][1] += halfMag * (nr.y*nr.y * BclmBsc +
                                                    mr.y*mr.y * BglmBsc);
                            Btotal[1][2] += halfMag * (nr.y*nr.z * BclmBsc +
                                                    mr.y*mr.z * BglmBsc);
                            Btotal[2][2] += halfMag * (nr.z*nr.z * BclmBsc +
                                                    mr.z*mr.z * BglmBsc);
                        }

                        FricForce += fricStress * bMag * mag;

                    }  // End non-[0 0 1] arm
                }  // End loop over arms

                Btotal[1][0] = Btotal[0][1];
                Btotal[2][0] = Btotal[0][2];
                Btotal[2][1] = Btotal[1][2];
                FricForce /= 2.0;

                if (numNonZeroLenSegs == 0 || fabs(Btotal.det()) < eps) {
                    vi = Vec3(0.0);
                } else {
                    Mat33 invDragMatrix = Btotal.inverse();

                    Vec3 f = fi;
                    if (fabs(FricForce) > eps) {
                        double fmag = f.norm();
                        if (fmag > FricForce) {
                            f -= FricForce/fmag * f;
                        } else {
                            f = Vec3(0.0);
                        }
                    }

                    vi = invDragMatrix * f;
                    if (params.vmax > 0.0)
                        apply_velocity_cap(params.vmax, vscale, vi);
                }

                double vmag = vi.norm();
                if (vmag > 0.0) vdir = 1.0/vmag * vi;
            }
        }
        
        return vi;
    }
    
    static constexpr const char* name = "MobilityBCC0c";
};

namespace MobilityType {
    typedef MobilityLocal<MobilityBCC0c> BCC_0C;
}

} // namespace ExaDiS

EXADIS_MOBILITY(MobilityBCC0c, BCC_0C)

#endif
