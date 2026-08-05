/*---------------------------------------------------------------------------
 *
 *	ExaDiS python binding registry mechanisms
 *
 *	Nicolas Bertin
 *	bertin1@llnl.gov
 *
 *-------------------------------------------------------------------------*/

#ifdef EXADIS_PYBIND

#pragma once
#ifndef EXADIS_PYBIND_REGISTRY_H
#define EXADIS_PYBIND_REGISTRY_H

#include <string>
#include <unordered_map>
#include <pybind11/pybind11.h>

namespace py = pybind11;

// Helper to make unique names
#define EXADIS_CONCAT_INNER(a, b) a##b
#define EXADIS_CONCAT(a, b) EXADIS_CONCAT_INNER(a, b)

/*---------------------------------------------------------------------------
 *
 *    Submodules registry
 *
 *-------------------------------------------------------------------------*/
namespace ExaDiS { namespace pybind {

using PybindInitModule = void(*)(py::module_&);
struct ModuleBindingRegistry {
    static std::vector<PybindInitModule>& registry() {
        static std::vector<PybindInitModule> r;
        return r;
    }

    static void add(PybindInitModule s) {
        registry().push_back(s);
    }

    static void bind_all(py::module_& m) {
        for (auto s : registry()) {
            s(m);
        }
    }
};

} } // namespace ExaDiS::pybind

#define EXADIS_MODULE_BINDING(BinderName)                                               \
    namespace ExaDiS { namespace pybind {                                               \
        static inline void EXADIS_CONCAT(bind_, BinderName)(py::module& m);             \
        struct EXADIS_CONCAT(ModuleRegistrar_, BinderName) {                            \
            EXADIS_CONCAT(ModuleRegistrar_, BinderName)() {                             \
                ModuleBindingRegistry::add(&EXADIS_CONCAT(bind_, BinderName));          \
            }                                                                           \
        } EXADIS_CONCAT(_module_registrar_, BinderName);                                \
    } }                                                                                 \
    static inline void ExaDiS::pybind::EXADIS_CONCAT(bind_, BinderName)(py::module& m)


/*---------------------------------------------------------------------------
 *
 *    Mobility registry
 *
 *-------------------------------------------------------------------------*/
namespace ExaDiS { namespace pybind {

// Registry that deduplicates by a string key
// for safety if headers are included in multiple TUs.
using PybindMobFn = void (*)(py::module&);
struct MobilityBindingRegistry {
    static std::unordered_map<std::string, PybindMobFn>& map() {
        static std::unordered_map<std::string, PybindMobFn> m;
        return m;
    }

    static void add(const std::string& name, PybindMobFn fn) {
        auto& m = map();
        if (m.find(name) == m.end()) m.emplace(name, fn);
    }

    static void bind_all(py::module& m) {
        for (auto& kv : map()) kv.second(m);
    }
};

// Registry for factories keyed by python mobility name
using FactoryFn = MobilityBind (*)(Params&, py::dict);
struct MobilityFactoryRegistry {
    static std::unordered_map<std::string, FactoryFn>& map() {
        static std::unordered_map<std::string, FactoryFn> m;
        return m;
    }
    static void add(const std::string& name, FactoryFn fn) {
        auto& m = map();
        if (m.find(name) == m.end()) m.emplace(name, fn);
    }
    static FactoryFn* get(const std::string& name) {
        auto& m = map();
        auto it = m.find(name);
        return it == m.end() ? nullptr : &it->second;
    }
};

} } // namespace ExaDiS::pybind

#define EXADIS_MOBILITY(BinderName, MobilityTypeT)                                            \
    namespace ExaDiS { namespace pybind {                                                     \
        /* Generic factory for name-based creation */                                         \
        static MobilityBind EXADIS_CONCAT(factory_, BinderName)(Params& params,               \
                                                                py::dict mobparams) {         \
            return make_mobility<MobilityType::MobilityTypeT>(params, mobparams);             \
        }                                                                                     \
        static inline void EXADIS_CONCAT(bind_, BinderName)(py::module& m);                   \
        struct EXADIS_CONCAT(MobilityRegistrar_, BinderName) {                                \
            EXADIS_CONCAT(MobilityRegistrar_, BinderName)() {                                 \
                MobilityBindingRegistry::add(#BinderName, &EXADIS_CONCAT(bind_, BinderName)); \
                MobilityFactoryRegistry::add(#MobilityTypeT,                                  \
                                             &EXADIS_CONCAT(factory_, BinderName));           \
            }                                                                                 \
        };                                                                                    \
        static EXADIS_CONCAT(MobilityRegistrar_, BinderName)                                  \
            EXADIS_CONCAT(_mobility_registrar_instance_, BinderName);                         \
    } }                                                                                       \
    static inline void ExaDiS::pybind::EXADIS_CONCAT(bind_, BinderName)(py::module& m){}


#endif // EXADIS_PYBIND_REGISTRY_H
#endif // EXADIS_PYBIND
