/*---------------------------------------------------------------------------
 *
 *	ExaDiS
 *
 *	Nicolas Bertin
 *	bertin1@llnl.gov
 *
 *-------------------------------------------------------------------------*/

#pragma once
#ifndef EXADIS_FORCE_GLOBAL_H
#define EXADIS_FORCE_GLOBAL_H

#include "force.h"
#include "force_global_list.h"

namespace ExaDiS {

/*---------------------------------------------------------------------------
 *
 *    Class:        ForceGlobal
 *
 *-------------------------------------------------------------------------*/
class ForceGlobal : public Force {
public:
    // Define types and enum
    #define X(TYPE, ALIAS) typedef TYPE ALIAS;
        EXADIS_FORCE_GLOBAL_LIST
    #undef X
    
    struct ForceList {
        enum ForceTypes {
        #define X(TYPE, ALIAS) ALIAS,
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
            FORCE_END
        };
    };

private:
    // Storage, one member per registered force type.
    #define X(TYPE, ALIAS) ALIAS* force_##ALIAS;
        EXADIS_FORCE_GLOBAL_LIST
    #undef X

    // Enabled flags per force, indexed by ForceList::ForceTypes.
    bool use_force[ForceList::FORCE_END] = {0};

    inline void add_force(int f) {
        if (use_force[f])
            ExaDiS_fatal("Error: force contribution already defined\n");
        use_force[f] = true;
    }

    // Map a force alias type to its enum id.
    template<typename T> struct ForceEnum;
    
public:
    // Accessor by type.
    template<typename T> inline T*& get();
    
    // Generic, type-safe add. Prefer this over the per-force wrappers.
    template<typename T>
    ForceGlobal* add(System* system, typename T::Params params = typename T::Params()) {
        add_force(ForceEnum<T>::value);
        get<T>() = exadis_new<T>(system, params);
        return this;
    }
    
    // Composite convenience function for FFT-based model.
    inline ForceGlobal* add_LONG_FFT_SHORT_ISO(System* system, FORCE_FFT::Params params);
    
    // Hooks
    void pre_compute(System* system) {
        #define X(TYPE, ALIAS) \
            if (use_force[ForceList::ALIAS]) force_##ALIAS->pre_compute(system);
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
    }

    void compute(System* system, bool zero=true) {
        DeviceDisNet* net = system->get_device_network();
        if (zero) zero_force(net);
        #define X(TYPE, ALIAS) \
            if (use_force[ForceList::ALIAS]) force_##ALIAS->compute(system, false);
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
    }

    Vec3 node_force(System* system, const int& i) {
        Vec3 f(0.0);
        #define X(TYPE, ALIAS) \
            if (use_force[ForceList::ALIAS]) f += force_##ALIAS->node_force(system, i);
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
        return f;
    }

    template<class N>
    KOKKOS_INLINE_FUNCTION
    Vec3 node_force(System* system, N* net, const int& i, const team_handle& team)
    {
        Vec3 f(0.0);
        #define X(TYPE, ALIAS) \
            if (use_force[ForceList::ALIAS]) f += force_##ALIAS->node_force(system, net, i, team);
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
        return f;
    }
    
    ~ForceGlobal() {
        #define X(TYPE, ALIAS) exadis_delete(force_##ALIAS);
            EXADIS_FORCE_GLOBAL_LIST
        #undef X
    }
    
    const char* name() { return "ForceGlobal"; }
};

// Enum specialization by type.
#define X(TYPE, ALIAS) \
    template<> struct ForceGlobal::ForceEnum<ForceGlobal::ALIAS> \
    { static constexpr int value = ForceGlobal::ForceList::ALIAS; };
    EXADIS_FORCE_GLOBAL_LIST
#undef X

// Accessor specialization by type.
#define X(TYPE, ALIAS) \
    template<> inline ForceGlobal::ALIAS*& ForceGlobal::get<ForceGlobal::ALIAS>() \
    { return force_##ALIAS; }
    EXADIS_FORCE_GLOBAL_LIST
#undef X

// Composite convenience function for FFT-based model
inline ForceGlobal* ForceGlobal::add_LONG_FFT_SHORT_ISO(System* system, ForceGlobal::FORCE_FFT::Params params) {
    add<FORCE_FFT>(system, params);
    // Build short-range correction from long-range FFT parameters
    add<FORCE_SEGSEG_ISO_FFT>(system, FORCE_SEGSEG_ISO_FFT::Params(
        get<FORCE_FFT>()->get_neighbor_cutoff(),
        SegSegIsoFFT::Params(get<FORCE_FFT>()->get_rcgrid())
    ));
    return this;
}

namespace ForceType {
    typedef ForceGlobal GLOBAL_MODEL;
}

} // namespace ExaDiS

#endif
